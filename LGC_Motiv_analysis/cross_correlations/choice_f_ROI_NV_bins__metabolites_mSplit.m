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
which_task = listdlg('PromptString','Which task?','ListString',task_names);
task_to_look = task_names{which_task};
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
nRuns = 4;
nTrialsPerRun = 54;
nTrials = nTrialsPerRun*nRuns;
% input/mediator/output
[input_prm, ROI_mediator, choice_hE] = deal(NaN(nTrials, NS));

% bin variables
[input_f_input_bin, ROI_f_input_bin, choice_f_input_bin,...
    input_f_input_low_ROI, ROI_f_input_low_ROI, choice_f_input_low_ROI,...
    input_f_input_high_ROI, ROI_f_input_high_ROI, choice_f_input_high_ROI] = deal(NaN(nBins, NS));

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
            error('please update for bayesian model');
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
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(kRun-1);
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        % define which task session it is
        switch kRun
            case {1,2}
                run_nm_bis = ['run',num2str(1)];
            case {3,4}
                run_nm_bis = ['run',num2str(2)];
        end
        
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
                otherwise
                    error('not ready yet');
            end
        end
        
        %% choice
        choice_LR_tmp = choiceAndPerf_tmp.choice;
        % remove confidence info from choice:
        choice_LR_tmp(choice_LR_tmp == 2) = 1;
        choice_LR_tmp(choice_LR_tmp == -2) = -1;
        % extract high effort choice
        choice_highE_tmp = NaN(1,length(choice_LR_tmp));
        choice_highE_tmp(choice_LR_tmp == -defaultSide_tmp) = 1;
        choice_highE_tmp(choice_LR_tmp == defaultSide_tmp) = 0;
        
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
        
        %% extract fMRI ROI mediator
        if strcmp(task_to_look,'EpEmPool') ||...
                (strcmp(task_to_look, task_nm))
            ROI_mediator(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
        end
        %% extract output behavioral variable
        choice_hE(runTrials_idx, iS) = choice_highE_tmp;
    end % run loop
    

    %% extract the bins
    % 1) split the data according to input parameter bins
    [ROI_f_input_bin(:,iS), input_f_input_bin(:,iS), bin_idx_tmp] = do_bin2(ROI_mediator(:, iS), input_prm(:,iS), nBins, 0);
    choice_f_input_bin(:,iS) = do_bin2(choice_hE(:, iS), input_prm(:,iS), nBins, 0);

    % 2) perform a median split according to the ROI activity for each bin
    for iBin = 1:nBins
        curr_bin_idx_tmp = bin_idx_tmp == iBin;
        med_tmp = median(ROI_mediator(curr_bin_idx_tmp, iS), 1,'omitnan');
        % extract data for the current bin of interest
        input_bin_tmp   = input_prm(curr_bin_idx_tmp, iS);
        ROI_bin_tmp     = ROI_mediator(curr_bin_idx_tmp, iS);
        choice_bin_tmp  = choice_hE(curr_bin_idx_tmp, iS);
        % extract data for low ROI activity
        curr_bin_low_ROI_idx = ROI_bin_tmp < med_tmp;
        input_f_input_low_ROI(iBin, iS) = mean(input_bin_tmp(curr_bin_low_ROI_idx),1,'omitnan');
        ROI_f_input_low_ROI(iBin, iS) = mean(ROI_bin_tmp(curr_bin_low_ROI_idx),1,'omitnan');
        choice_f_input_low_ROI(iBin, iS) =  mean(choice_bin_tmp(curr_bin_low_ROI_idx),1,'omitnan');
        % extract data for high ROI activity
        curr_bin_high_ROI_idx = ROI_bin_tmp > med_tmp;
        input_f_input_high_ROI(iBin, iS) = mean(input_bin_tmp(curr_bin_high_ROI_idx),1,'omitnan');
        ROI_f_input_high_ROI(iBin, iS) = mean(ROI_bin_tmp(curr_bin_high_ROI_idx),1,'omitnan');
        choice_f_input_high_ROI(iBin, iS) =  mean(choice_bin_tmp(curr_bin_high_ROI_idx),1,'omitnan');
    end
    
end % subject loop

%% load group indexes depending on the level of metabolites
[met_subs.low, met_subs.high, metabolite_nm, ROI_nm] = medSplit_metabolites(study_nm, subject_id);
met_groups = {'low','high'};
n_met_grps = length(met_groups);
%% average data according to which group (low/high metabolite)
% basic data
for iGrp = 1:n_met_grps
    grp_nm = met_groups{iGrp};
    grp_nm_bis = [grp_nm,'_',metabolite_nm];
    % extract data
    [m_input_f_input.(grp_nm_bis),...
        sem_input_f_input.(grp_nm_bis)] = mean_sem_sd(input_f_input_bin(:,met_subs.(grp_nm)), 2);
    [m_ROI_f_input.(grp_nm_bis),...
        sem_ROI_f_input.(grp_nm_bis)] = mean_sem_sd(ROI_f_input_bin(:,met_subs.(grp_nm)), 2);
    [m_choice_f_input.(grp_nm_bis),...
        sem_choice_f_input.(grp_nm_bis)] = mean_sem_sd(choice_f_input_bin(:,met_subs.(grp_nm)), 2);
    % low ROI median split
    [m_input_f_input_low_ROI.(grp_nm_bis),...
        sem_input_f_input_low_ROI.(grp_nm_bis)] = mean_sem_sd(input_f_input_low_ROI(:,met_subs.(grp_nm)), 2);
    [m_ROI_f_input_low_ROI.(grp_nm_bis),...
        sem_ROI_f_input_low_ROI.(grp_nm_bis)] = mean_sem_sd(ROI_f_input_low_ROI(:,met_subs.(grp_nm)), 2);
    [m_choice_f_input_low_ROI.(grp_nm_bis),...
        sem_choice_f_input_low_ROI.(grp_nm_bis)] = mean_sem_sd(choice_f_input_low_ROI(:,met_subs.(grp_nm)), 2);
    % high ROI median split
    [m_input_f_input_high_ROI.(grp_nm_bis),...
        sem_input_f_input_high_ROI.(grp_nm_bis)] = mean_sem_sd(input_f_input_high_ROI(:,met_subs.(grp_nm)), 2);
    [m_ROI_f_input_high_ROI.(grp_nm_bis),...
        sem_ROI_f_input_high_ROI.(grp_nm_bis)] = mean_sem_sd(ROI_f_input_high_ROI(:,met_subs.(grp_nm)), 2);
    [m_choice_f_input_high_ROI.(grp_nm_bis),...
        sem_choice_f_input_high_ROI.(grp_nm_bis)] = mean_sem_sd(choice_f_input_high_ROI(:,met_subs.(grp_nm)), 2);
end % loop median split

%% figures
if dispFig == true
    lWidth = 3;
    black = [0 0 0];
    grey = [143 143 143]./255;
    purple = [153 142 195]./255;
    orange = [241 163 64]./255;
    pSize = 50;
    input_prm_nm = strrep(input_prm_nm,'_',' ');
    low_met_nm = ['low_',metabolite_nm];
    high_met_nm = ['high_',metabolite_nm];
    low_met_nm_bis = ['low ',metabolite_nm];
    high_met_nm_bis = ['high ',metabolite_nm];
    
    % look at the general figure (choice = f(inputs), ROI=f(inputs)
    fig;
    % choices = f(inputs)
    subplot(1,2,1);
    hold on;
    gal_data_hdl_low_met = errorbar(m_input_f_input.(low_met_nm),...
        m_choice_f_input.(low_met_nm),...
        sem_choice_f_input.(low_met_nm));
    gal_data_hdl_high_met = errorbar(m_input_f_input.(high_met_nm),...
        m_choice_f_input.(high_met_nm),...
        sem_choice_f_input.(high_met_nm));
    gal_data_hdl_low_met.LineStyle = '--';
    gal_data_hdl_low_met.LineWidth = lWidth;
    gal_data_hdl_low_met.MarkerEdgeColor = grey;
    gal_data_hdl_high_met.LineStyle = '-';
    gal_data_hdl_high_met.LineWidth = lWidth;
    gal_data_hdl_high_met.MarkerEdgeColor = black;
    legend([gal_data_hdl_high_met, gal_data_hdl_low_met],...
        {high_met_nm_bis, low_met_nm_bis});
    legend('boxoff');
    xlabel(input_prm_nm);
    ylabel('Choice = high effort (%)');
    legend_size(pSize);

    subplot(1,2,2);
    % ROI = f(inputs)
    hold on;
    ROI_f_input_hdl_low_met = errorbar(m_input_f_input.(low_met_nm),...
        m_ROI_f_input.(low_met_nm),...
        sem_ROI_f_input.(low_met_nm));
    ROI_f_input_hdl_high_met = errorbar(m_input_f_input.(high_met_nm),...
        m_ROI_f_input.(high_met_nm),...
        sem_ROI_f_input.(high_met_nm));
    ROI_f_input_hdl_low_met.LineStyle = '--';
    ROI_f_input_hdl_low_met.LineWidth = lWidth;
    ROI_f_input_hdl_low_met.MarkerEdgeColor = grey;
    ROI_f_input_hdl_high_met.LineStyle = '-';
    ROI_f_input_hdl_high_met.LineWidth = lWidth;
    ROI_f_input_hdl_high_met.MarkerEdgeColor = black;
    legend([ROI_f_input_hdl_high_met, ROI_f_input_hdl_low_met],...
        {high_met_nm_bis, low_met_nm_bis});
    legend('boxoff');
    xlabel(input_prm_nm);
    ylabel([ROI_short_nm,' BOLD']);
    legend_size(pSize);

    %% sanity check: ROI = f(inputs) + ROI median split (did it work?)
    fig;
    hold on;
    ROI_f_input_low_ROI_low_met_hdl = errorbar(m_input_f_input_low_ROI.(low_met_nm),...
        m_ROI_f_input_low_ROI.(low_met_nm),...
        sem_ROI_f_input_low_ROI.(low_met_nm));
    ROI_f_input_high_ROI_low_met_hdl = errorbar(m_input_f_input_high_ROI.(low_met_nm),...
        m_ROI_f_input_high_ROI.(low_met_nm),...
        sem_ROI_f_input_high_ROI.(low_met_nm));
    ROI_f_input_low_ROI_high_met_hdl = errorbar(m_input_f_input_low_ROI.(high_met_nm),...
        m_ROI_f_input_low_ROI.(high_met_nm),...
        sem_ROI_f_input_low_ROI.(high_met_nm));
    ROI_f_input_high_ROI_high_met_hdl = errorbar(m_input_f_input_high_ROI.(high_met_nm),...
        m_ROI_f_input_high_ROI.(high_met_nm),...
        sem_ROI_f_input_high_ROI.(high_met_nm));
    ROI_f_input_low_ROI_low_met_hdl.LineStyle = '--';
    ROI_f_input_low_ROI_low_met_hdl.LineWidth = lWidth;
    ROI_f_input_low_ROI_low_met_hdl.MarkerEdgeColor = purple;
    ROI_f_input_high_ROI_low_met_hdl.LineStyle = '--';
    ROI_f_input_high_ROI_low_met_hdl.LineWidth = lWidth;
    ROI_f_input_high_ROI_low_met_hdl.MarkerEdgeColor = orange;
    ROI_f_input_low_ROI_high_met_hdl.LineStyle = '-';
    ROI_f_input_low_ROI_high_met_hdl.LineWidth = lWidth;
    ROI_f_input_low_ROI_high_met_hdl.MarkerEdgeColor = purple;
    ROI_f_input_high_ROI_high_met_hdl.LineStyle = '-';
    ROI_f_input_high_ROI_high_met_hdl.LineWidth = lWidth;
    ROI_f_input_high_ROI_high_met_hdl.MarkerEdgeColor = orange;
    xlabel(input_prm_nm);
    ylabel([ROI_short_nm,' BOLD']);
    legend([ROI_f_input_low_ROI_low_met_hdl, ROI_f_input_high_ROI_low_met_hdl,...
        ROI_f_input_low_ROI_high_met_hdl, ROI_f_input_high_ROI_high_met_hdl],...
        {['low ',ROI_short_nm,' - low ',metabolite_nm],...
        ['high ',ROI_short_nm,' - low ',metabolite_nm],...
        ['low ',ROI_short_nm,' - high ',metabolite_nm],...
        ['high ',ROI_short_nm,' - high ',metabolite_nm]});
    legend('boxoff');
    legend_size(pSize);

    %% main figure: choices = f(inputs) + ROI median split
    fig;
    hold on;
    choices_f_inputs_low_ROI_data_low_met_hdl = errorbar(m_input_f_input_low_ROI.(low_met_nm),...
        m_choice_f_input_low_ROI.(low_met_nm),...
        sem_choice_f_input_low_ROI.(low_met_nm));
    choices_f_inputs_high_ROI_data_low_met_hdl = errorbar(m_input_f_input_high_ROI.(low_met_nm),...
        m_choice_f_input_high_ROI.(low_met_nm),...
        sem_choice_f_input_high_ROI.(low_met_nm));
    choices_f_inputs_low_ROI_data_high_met_hdl = errorbar(m_input_f_input_low_ROI.(high_met_nm),...
        m_choice_f_input_low_ROI.(high_met_nm),...
        sem_choice_f_input_low_ROI.(high_met_nm));
    choices_f_inputs_high_ROI_data_high_met_hdl = errorbar(m_input_f_input_high_ROI.(high_met_nm),...
        m_choice_f_input_high_ROI.(high_met_nm),...
        sem_choice_f_input_high_ROI.(high_met_nm));
    choices_f_inputs_low_ROI_data_low_met_hdl.LineStyle = '--';
    choices_f_inputs_low_ROI_data_low_met_hdl.LineWidth = lWidth;
    choices_f_inputs_low_ROI_data_low_met_hdl.Color = purple;
    choices_f_inputs_high_ROI_data_low_met_hdl.LineStyle = '--';
    choices_f_inputs_high_ROI_data_low_met_hdl.LineWidth = lWidth;
    choices_f_inputs_high_ROI_data_low_met_hdl.Color = orange;
    choices_f_inputs_low_ROI_data_high_met_hdl.LineStyle = '-';
    choices_f_inputs_low_ROI_data_high_met_hdl.LineWidth = lWidth;
    choices_f_inputs_low_ROI_data_high_met_hdl.Color = purple;
    choices_f_inputs_high_ROI_data_high_met_hdl.LineStyle = '-';
    choices_f_inputs_high_ROI_data_high_met_hdl.LineWidth = lWidth;
    choices_f_inputs_high_ROI_data_high_met_hdl.Color = orange;
    legend([choices_f_inputs_low_ROI_data_low_met_hdl,...
        choices_f_inputs_high_ROI_data_low_met_hdl,...
        choices_f_inputs_low_ROI_data_high_met_hdl,...
        choices_f_inputs_high_ROI_data_high_met_hdl],...
        {['low ',ROI_short_nm,' - low ',metabolite_nm],...
        ['high ',ROI_short_nm,' - low ',metabolite_nm],...
        ['low ',ROI_short_nm,' - high ',metabolite_nm],...
        ['high ',ROI_short_nm,' - high ',metabolite_nm]});
    legend('boxoff');
    xlabel(input_prm_nm);
    ylabel('Choice = high effort (%)');
    legend_size(pSize);

    
end

%% compare slopes: should do an ANOVA or GLM with dmPFC and GSH as factors and look at the interaction(s)

% [slope_low, slope_high,...
%     intercept_low, intercept_high] = deal(NaN(1,NS));
% for iS = 1:NS
%     b1 = glmfit(1:nBins, choice_f_input_low_ROI(:,iS),'normal');
%     intercept_low(iS) = b1(1);
%     slope_low(iS) = b1(2);
%     b2 = glmfit(1:nBins, choice_f_input_high_ROI(:,iS),'normal');
%     intercept_high(iS) = b2(1);
%     slope_high(iS) = b2(2);
% end
% [~,pval.intercept]=ttest(intercept_low,intercept_high);
% [~,pval.slope]=ttest(slope_low,slope_high);