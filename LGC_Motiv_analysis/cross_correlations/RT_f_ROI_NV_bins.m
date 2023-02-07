%% look whether within each input parameter bin, the activity of the selected ROI 
% determines which choices should be made.
%

%% general parameters
% display figures?
dispFig = true;

%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% selection of participants
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% extract ROI activity for all subjects
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);
% define which ROI, and which time period is of interest to you
% define ROI
ROI_names = fieldnames(ROI_trial_b_trial);
ROI_subList = ROI_trial_b_trial.subject_id;
ROI_names(strcmp(ROI_names,'subject_id')) = [];
if length(ROI_names) > 1
    error(['There should be only 1 ROI selected, not ',num2str(length(ROI_names))])
else
    ROI_nm = ROI_names;
end
ROI_short_nm = inputdlg('ROI short name?');
ROI_short_nm = ROI_short_nm{1};
% define task
task_names = {'Ep','Em','EpEmPool'};
which_ROI_task = listdlg('PromptString','Which task for ROI?','ListString',task_names);
ROI_task_to_look = task_names{which_ROI_task};
% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};

%% select parameters of interest
potential_input_prm = {'NV_hE','deltaNV','uncertainty','E_level','money_level',...
    'deltaMoney_level'};
which_input = listdlg('PromptString','please select input parameter',...
    'ListString',potential_input_prm);
input_prm_nm = potential_input_prm{which_input};
which_bhv_task = listdlg('PromptString','Which task for behavior?','ListString',task_names);
behavioral_task_to_look = task_names{which_bhv_task};

%% number of bins
switch input_prm_nm
    case {'E_level','money_level','deltaMoney_level'}
        nBins = 3;
    case {'NV_hE','deltaNV','uncertainty'}
        nBins = 6;
    otherwise
        nBins = 6;
end

%% select which model to use (if relevant)
if ismember(input_prm_nm,{'NV_hE','deltaNV','uncertainty'})
    needModeling = true;
    [mdlType, mdlN] = behavioral_model_selection;
else
    needModeling = false;
end

%% extract the parameters of interest
switch behavioral_task_to_look
    case {'Ep','Em'}
        nRunsPerTask = 2;
        nTrialsPerRun = 54;
        nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
        [input_prm, ROI_mediator, RT] = deal(NaN(nTrialsPerTask, NS));
    case 'EpEmPool'
        nRuns = 4;
        nTrialsPerRun = 54;
        nTrials = nTrialsPerRun*nRuns;
        % input/mediator/output
        [input_prm, ROI_mediator, RT] = deal(NaN(nTrials, NS));
end

% bin variables
[input_f_input_bin, ROI_f_input_bin, RT_f_input_bin,...
    input_f_input_low_ROI, ROI_f_input_low_ROI, RT_f_input_low_ROI,...
    input_f_input_high_ROI, ROI_f_input_high_ROI, RT_f_input_high_ROI] = deal(NaN(nBins, NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];

    % extract confidence based on the model
    if needModeling == true
        switch mdlType
            case 'simple'
                [~, dataInferred] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                    0, 'levels', 6, 6);
            case 'bayesian'
                gitResultsFolder = [fullfile('C:','Users','clairis','Desktop',...
                    'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];
        end
    end

    % extract runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    jRun = 0;
    for iRun = 1:length(okRuns)
        kRun = okRuns(iRun);
        run_nm = num2str(kRun);
        jRun = jRun + 1;
        task_nm_tmp = taskNames{jRun};
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        % define which task session it is
        switch kRun
            case {1,2}
                taskRun_idx = 1;
            case {3,4}
                taskRun_idx = 2;
        end
        run_nm_bis = ['run',num2str(taskRun_idx)];
        % define trial index for relevant variable to extract
        switch behavioral_task_to_look
            case {'Ep','Em'}
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(taskRun_idx-1);
            case 'EpEmPool'
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(kRun-1);
        end
        
        % filter task based on what was selected in the inputs
        if strcmp(behavioral_task_to_look,'EpEmPool') ||...
                strcmp(behavioral_task_to_look, task_nm_tmp)
            %% load the data
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
                '_task.mat']);
            choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
            switch task_nm_tmp
                case 'Em'
                    choiceAndPerf_tmp = behaviorStruct_tmp.mentalE_perf;
                case 'Ep'
                    choiceAndPerf_tmp = behaviorStruct_tmp.physicalPerf;
            end
            
            %% default side
            defaultSide_tmp = choiceOptions_tmp.default_LR;
            %% extract R or P
            RP_var_tmp = strcmp(choiceOptions_tmp.R_or_P,'R');
            
            %% effort level
            E_highE_tmp = (choiceOptions_tmp.E.left).*(defaultSide_tmp == 1) +...
                (choiceOptions_tmp.E.right).*(defaultSide_tmp == -1);
            
            %% high-effort money amount
            money_hE_tmp = ((choiceOptions_tmp.monetary_amount.left).*(defaultSide_tmp == 1) +...
                (choiceOptions_tmp.monetary_amount.right).*(defaultSide_tmp == -1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
            money_lE_tmp = ((choiceOptions_tmp.monetary_amount.left).*(defaultSide_tmp == -1) +...
                (choiceOptions_tmp.monetary_amount.right).*(defaultSide_tmp == 1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
            
            %% delta between high and low effort options
            deltaMoney_tmp = money_hE_tmp - money_lE_tmp;
            
            %% RT
            onsets_tmp = behaviorStruct_tmp.onsets;
            switch task_nm_tmp
                case 'Ep'
                    onsets_tmp = behaviorStruct_tmp.physicalPerf.onsets;
                    choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                case 'Em'
                    onsets_tmp = behaviorStruct_tmp.mentalE_perf.onsets;
                    choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
            end
            RT_tmp = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
            
            %% confidence rating
            uncertaintyRtg_tmp = NaN(1,length(choice_LR_tmp));
            uncertaintyRtg_tmp(abs(choice_LR_tmp) == 2) = 0; % high confidence = low uncertainty
            uncertaintyRtg_tmp(abs(choice_LR_tmp) == 1) = 1; % low confidence = high uncertainty
            
            %% net value and confidence inferred by the model
            if needModeling == 1
                switch mdlType
                    case 'simple'
                        NV_hE_tmp = dataInferred.NV_varOption.(task_nm_tmp).(['mdl_',mdlN]).(run_nm_bis);
                        trial_idx = (1:nTrialsPerRun) + nTrialsPerRun*(kRun >= 3);
                        deltaNV_tmp = dataInferred.deltaNV.(['mdl_',mdlN]).(task_nm_tmp)(trial_idx);
                        uncertainty_tmp = - dataInferred.confidenceFitted.(['mdl_',mdlN]).(run_nm_bis); % revert sign to transform confidence into uncertainty
                    case 'bayesian'
                        [~, NV_hE_tmp, confidence_tmp] = extract_bayesian_mdl(gitResultsFolder, subBehaviorFolder,...
                            sub_nm, run_nm, task_fullName, ['mdl_',mdlN]);
                        uncertainty_tmp = -confidence_tmp;
                end
            end
            
            %% extract input behavioral variable
            switch input_prm_nm
                case 'uncertainty'
                    input_prm(runTrials_idx, iS) = uncertainty_tmp;
                case 'E_level'
                    input_prm(runTrials_idx, iS) = E_highE_tmp;
                case 'money_level'
                    input_prm(runTrials_idx, iS) = money_hE_tmp;
                case 'deltaMoney_level'
                    input_prm(runTrials_idx, iS) = deltaMoney_tmp;
                case 'deltaNV'
                    input_prm(runTrials_idx, iS) = deltaNV_tmp;
                case 'NV_hE'
                    input_prm(runTrials_idx, iS) = NV_hE_tmp;
                otherwise
                    error(['input = ',input_prm_nm,' not ready yet']);
            end
            
            %% extract RT = output behavioral variable
            RT(runTrials_idx, iS) = RT_tmp;
        end % task filter
        %% extract fMRI ROI mediator
        if strcmp(ROI_task_to_look,'EpEmPool') ||...
                (strcmp(ROI_task_to_look, task_nm_tmp))
            ROI_mediator(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
        end
        
    end % run loop
    

    %% extract the bins
    % 1) split the data according to input parameter bins
    [ROI_f_input_bin(:,iS), input_f_input_bin(:,iS), bin_idx_tmp] = do_bin2(ROI_mediator(:, iS), input_prm(:,iS), nBins, 0);
    RT_f_input_bin(:,iS) = do_bin2(RT(:, iS), input_prm(:,iS), nBins, 0);

    % 2) perform a median split according to the ROI activity for each bin
    for iBin = 1:nBins
        curr_bin_idx_tmp = bin_idx_tmp == iBin;
        med_tmp = median(ROI_mediator(curr_bin_idx_tmp, iS), 1,'omitnan');
        % extract data for the current bin of interest
        input_bin_tmp   = input_prm(curr_bin_idx_tmp, iS);
        ROI_bin_tmp     = ROI_mediator(curr_bin_idx_tmp, iS);
        RT_bin_tmp  = RT(curr_bin_idx_tmp, iS);
        % extract data for low ROI activity
        curr_bin_low_ROI_idx = ROI_bin_tmp < med_tmp;
        input_f_input_low_ROI(iBin, iS) = mean(input_bin_tmp(curr_bin_low_ROI_idx),1,'omitnan');
        ROI_f_input_low_ROI(iBin, iS) = mean(ROI_bin_tmp(curr_bin_low_ROI_idx),1,'omitnan');
        RT_f_input_low_ROI(iBin, iS) =  mean(RT_bin_tmp(curr_bin_low_ROI_idx),1,'omitnan');
        % extract data for high ROI activity
        curr_bin_high_ROI_idx = ROI_bin_tmp > med_tmp;
        input_f_input_high_ROI(iBin, iS) = mean(input_bin_tmp(curr_bin_high_ROI_idx),1,'omitnan');
        ROI_f_input_high_ROI(iBin, iS) = mean(ROI_bin_tmp(curr_bin_high_ROI_idx),1,'omitnan');
        RT_f_input_high_ROI(iBin, iS) =  mean(RT_bin_tmp(curr_bin_high_ROI_idx),1,'omitnan');
    end
    
end % subject loop

%% average data
% basic data
[m_input_f_input,...
    sem_input_f_input] = mean_sem_sd(input_f_input_bin, 2);
[m_ROI_f_input,...
    sem_ROI_f_input] = mean_sem_sd(ROI_f_input_bin, 2);
[m_RT_f_input,...
    sem_RT_f_input] = mean_sem_sd(RT_f_input_bin, 2);
% low ROI median split
[m_input_f_input_low_ROI,...
    sem_input_f_input_low_ROI] = mean_sem_sd(input_f_input_low_ROI, 2);
[m_ROI_f_input_low_ROI,...
    sem_ROI_f_input_low_ROI] = mean_sem_sd(ROI_f_input_low_ROI, 2);
[m_RT_f_input_low_ROI,...
    sem_RT_f_input_low_ROI] = mean_sem_sd(RT_f_input_low_ROI, 2);
% high ROI median split
[m_input_f_input_high_ROI,...
    sem_input_f_input_high_ROI] = mean_sem_sd(input_f_input_high_ROI, 2);
[m_ROI_f_input_high_ROI,...
    sem_ROI_f_input_high_ROI] = mean_sem_sd(ROI_f_input_high_ROI, 2);
[m_RT_f_input_high_ROI,...
    sem_RT_f_input_high_ROI] = mean_sem_sd(RT_f_input_high_ROI, 2);

%% figures
if dispFig == true
    lWidth = 3;
    black = [0 0 0];
    purple = [153 142 195]./255;
    orange = [241 163 64]./255;
    pSize = 50;
    input_prm_nm = strrep(input_prm_nm,'_',' ');
    switch behavioral_task_to_look
        case 'Ep'
            full_bhv_taskName = 'physical task';
        case 'Em'
            full_bhv_taskName = 'mental task';
        case 'EpEmPool'
            full_bhv_taskName = 'both tasks';
    end
    switch ROI_task_to_look
        case 'Ep'
            full_ROI_taskName = 'physical task';
        case 'Em'
            full_ROI_taskName = 'mental task';
        case 'EpEmPool'
            full_ROI_taskName = 'both tasks';
    end
    
    % look at the general figure (choice = f(inputs), ROI=f(inputs)
    fig;
    % RT = f(inputs)
    subplot(1,2,1);
    hold on;
    gal_data_hdl = errorbar(m_input_f_input,...
        m_RT_f_input,...
        sem_RT_f_input);
    gal_data_hdl.LineStyle = 'none';
    gal_data_hdl.LineWidth = lWidth;
    gal_data_hdl.MarkerEdgeColor = black;
    xlabel({input_prm_nm; full_bhv_taskName});
    ylabel({'RT (s)';...
        full_bhv_taskName});
    legend_size(pSize);

    subplot(1,2,2);
    % ROI = f(inputs)
    hold on;
    ROI_f_input_hdl = errorbar(m_input_f_input,...
        m_ROI_f_input,...
        sem_ROI_f_input);
    ROI_f_input_hdl.LineStyle = 'none';
    ROI_f_input_hdl.LineWidth = lWidth;
    ROI_f_input_hdl.MarkerEdgeColor = black;
    xlabel({input_prm_nm; full_bhv_taskName});
    ylabel({[ROI_short_nm,' BOLD '];full_ROI_taskName});
    legend_size(pSize);

    %% sanity check: ROI = f(inputs) + ROI median split (did it work?)
    fig;
    hold on;
    ROI_f_input_low_ROI_hdl = errorbar(m_input_f_input_low_ROI,...
        m_ROI_f_input_low_ROI,...
        sem_ROI_f_input_low_ROI);
    ROI_f_input_high_ROI_hdl = errorbar(m_input_f_input_high_ROI,...
        m_ROI_f_input_high_ROI,...
        sem_ROI_f_input_high_ROI);
    ROI_f_input_low_ROI_hdl.LineStyle = 'none';
    ROI_f_input_low_ROI_hdl.LineWidth = lWidth;
    ROI_f_input_low_ROI_hdl.MarkerEdgeColor = purple;
    ROI_f_input_high_ROI_hdl.LineStyle = 'none';
    ROI_f_input_high_ROI_hdl.LineWidth = lWidth;
    ROI_f_input_high_ROI_hdl.MarkerEdgeColor = orange;
    xlabel({input_prm_nm; full_bhv_taskName});
    ylabel({[ROI_short_nm,' BOLD '];full_ROI_taskName});
    legend([ROI_f_input_low_ROI_hdl, ROI_f_input_high_ROI_hdl],...
        {['low ',ROI_short_nm],['high ',ROI_short_nm]});
    legend('boxoff');
    legend_size(pSize);

    %% main figure: RT = f(inputs) + ROI median split
    fig;
    hold on;
    RT_f_inputs_low_ROI_data_hdl = errorbar(m_input_f_input_low_ROI,...
        m_RT_f_input_low_ROI,...
        sem_RT_f_input_low_ROI);
    RT_f_inputs_high_ROI_data_hdl = errorbar(m_input_f_input_high_ROI,...
        m_RT_f_input_high_ROI,...
        sem_RT_f_input_high_ROI);
    RT_f_inputs_low_ROI_data_hdl.LineWidth = lWidth;
    RT_f_inputs_low_ROI_data_hdl.MarkerEdgeColor = purple;
    RT_f_inputs_high_ROI_data_hdl.LineWidth = lWidth;
    RT_f_inputs_high_ROI_data_hdl.MarkerEdgeColor = orange;
    xlabel({input_prm_nm; full_bhv_taskName});
    ylabel({'RT (s)';...
        full_bhv_taskName});
    legend([RT_f_inputs_low_ROI_data_hdl, RT_f_inputs_high_ROI_data_hdl],...
        {['low ',ROI_short_nm],['high ',ROI_short_nm]});
    legend('boxoff');
    legend_size(pSize);
end

%% compare slopes

[slope_low, slope_high,...
    intercept_low, intercept_high] = deal(NaN(1,NS));
for iS = 1:NS
    b1 = glmfit(1:nBins, RT_f_input_low_ROI(:,iS),'normal');
    intercept_low(iS) = b1(1);
    slope_low(iS) = b1(2);
    b2 = glmfit(1:nBins, RT_f_input_high_ROI(:,iS),'normal');
    intercept_high(iS) = b2(1);
    slope_high(iS) = b2(2);
end
[~,pval.intercept]=ttest(intercept_low,intercept_high);
[~,pval.slope]=ttest(slope_low,slope_high);