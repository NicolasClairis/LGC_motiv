%% script looking at whether the brain activity between two brain areas is 
% correlated across trials. You may also want to look at DCM eventually

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
if length(ROI_names) ~= 2
    if length(ROI_names) < 2
        error('number of selected ROIs is too low, should be = 2');
    else
        error('number of selected ROIs is too high, should be = 2');
    end
else
    ROI_nm = ROI_names;
end
ROI_short_nm1 = inputdlg(['short name for ',ROI_names{1},'?']);
ROI_short_nm2 = inputdlg(['short name for ',ROI_names{2},'?']);
ROI_short_nm1 = ROI_short_nm1{1};
ROI_short_nm2 = ROI_short_nm2{1};
% define task
task_names = {'Ep','Em','EpEmPool'};
which_task = listdlg('PromptString','Which task?','ListString',task_names);
task_to_look = task_names{which_task};
% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(ROI_nm{1}).Ep.run1);
which_timePeriod1 = listdlg('PromptString',...
    ['trial period for ',ROI_short_nm1,'?'],...
    'listString',timePeriods);
timePeriod_nm1 = timePeriods{which_timePeriod1};
which_timePeriod2 = listdlg('PromptString',...
    ['trial period for ',ROI_short_nm2,'?'],...
    'listString',timePeriods);
timePeriod_nm2 = timePeriods{which_timePeriod2};

%% select parameters of interest
potential_ROI1_prm = {'uncertainty','E_level','money_level',...
    'deltaMoney_level','choice_hE','E_chosen'};
potential_ROI2_prm = {'RT','uncertainty_rtg','choice_hE','E_chosen'};

which_prm_ROI1 = listdlg('PromptString',['select prm for ',ROI_short_nm1],...
    'ListString',potential_ROI1_prm);
ROI1_prm_nm = potential_ROI1_prm{which_prm_ROI1};
which_prm_ROI2 = listdlg('PromptString',['select prm for ',ROI_short_nm2],...
    'ListString',potential_ROI2_prm);
ROI2_prm_nm = potential_ROI2_prm{which_prm_ROI2};

%% if a parameter is based on modeling, you need to decide which model to use
% bayesian across tasks or simple model each task separately?
if ismember(ROI1_prm_nm,{'uncertainty'})
    needModeling = 1;
    if ~exist('mdlType','var') || isempty(mdlType)
        listPossibleModels = {'bayesian','simple'};
        mdlType_idx = listdlg('promptstring','Which model type?',...
            'ListString',listPossibleModels);
        mdlType = listPossibleModels{mdlType_idx};
    end
    % warning until this gets fixed
    if strcmp(mdlType,'bayesian')
        error('uncertainty extraction from bayesian model not ready yet.');
    end
    
    % which model number to use?
    if ~exist('mdlN','var') || isempty(mdlN)
        switch mdlType
            case 'bayesian'
                listPossibleModelNumbers = {'1','3'};
            case 'simple'
                listPossibleModelNumbers = {'1','2','3','4'};
        end
        mdlN_idx = listdlg('promptstring','Which model number?',...
            'ListString',listPossibleModelNumbers);
        mdlN = listPossibleModelNumbers{mdlN_idx};
    end
else
    needModeling = 0;
end

%% extract the parameters of interest
nRuns = 4;
nTrialsPerRun = 54;
nTrials = nTrialsPerRun*nRuns;
% input/mediator/output
[ROI1_prm, ROI1_activity, ROI2_activity, ROI2_prm,...
    ROI1_activity_fit, ROI2_activity_fit] = deal(NaN(nTrials, NS));
% bin variables
nBins = 6;
[ROI1_prm_f_ROI1_prm_bin, ROI1_f_ROI1_prm_bin,...
    ROI1_fit_f_ROI1_prm_bin,...
    ROI2_prm_f_ROI2_prm_bin, ROI2_f_ROI2_prm_bin,...
    ROI2_fit_f_ROI2_prm_bin,...
    ROI2_f_ROI1_bin, ROI1_f_ROI1_bin] = deal(NaN(nBins, NS));
[betas_ROI1, betas_ROI2,...
    betas_ROI1_vs_ROI2_slopes] = deal(NaN(2,NS));

for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
    % extract confidence based on the model
    if needModeling == 1
        if strcmp(mdlType,'simple')
            [~, dataInferred] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                0, 'levels', 6, 6);
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
        
        %% choice
        choice_LR_tmp = choiceAndPerf_tmp.choice;
        % remove confidence info from choice:
        choice_LR_tmp(choice_LR_tmp == 2) = 1;
        choice_LR_tmp(choice_LR_tmp == -2) = -1;
        % extract high effort choice
        choice_highE_tmp = NaN(1,length(choice_LR_tmp));
        choice_highE_tmp(choice_LR_tmp == -defaultSide_tmp) = 1;
        choice_highE_tmp(choice_LR_tmp == defaultSide_tmp) = 0;
        
        %% effort level
        E_highE_tmp = (choiceOptions_tmp.E.left).*(defaultSide_tmp == 1) +...
            (choiceOptions_tmp.E.right).*(defaultSide_tmp == -1);
        E_chosen_tmp = (choiceOptions_tmp.E.left).*(choice_LR_tmp == -1) +...
            (choiceOptions_tmp.E.right).*(choice_LR_tmp == 1);
        
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
        
        %% confidence inferred by the model
        if needModeling == 1 
            switch mdlType
                case 'simple'
                    uncertainty_tmp = - dataInferred.confidenceFitted.(['mdl_',mdlN]).(run_nm_bis); % revert sign to transform confidence into uncertainty
                otherwise
                    error('not ready yet');
            end
        end
        
        %% extract input behavioral variable
        switch ROI1_prm_nm
            case 'uncertainty'
                ROI1_prm(runTrials_idx, iS) = uncertainty_tmp;
            case 'E_level'
                ROI1_prm(runTrials_idx, iS) = E_highE_tmp;
            case 'money_level'
                ROI1_prm(runTrials_idx, iS) = money_hE_tmp;
            case 'deltaMoney_level'
                ROI1_prm(runTrials_idx, iS) = deltaMoney_tmp;
            case 'choice_hE'
                ROI1_prm(runTrials_idx, iS) = choice_highE_tmp;
            case 'E_chosen'
                ROI1_prm(runTrials_idx, iS) = E_chosen_tmp;
            otherwise
                error(['input = ',ROI1_prm_nm,' not ready yet']);
        end
        
        %% extract fMRI ROI mediator
        if strcmp(task_to_look,'EpEmPool') ||...
                (strcmp(task_to_look, task_nm))
            ROI1_activity(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm1)(:, iS);
            ROI2_activity(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{2}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm2)(:, iS);
        end
        %% extract output behavioral variable
        switch ROI2_prm_nm
            case 'RT'
                ROI2_prm(runTrials_idx, iS) = RT_tmp;
            case 'uncertainty_rtg'
                ROI2_prm(runTrials_idx, iS) = uncertaintyRtg_tmp;
            case 'choice_hE'
                ROI2_prm(runTrials_idx, iS) = choice_highE_tmp;
            case 'E_chosen'
                ROI2_prm(runTrials_idx, iS) = E_chosen_tmp;
            otherwise
                error(['output = ',ROI2_prm_nm,' not ready yet']);
        end
    end % run loop
    
    %% do the bins per subject
    % ROI 1 = f(prm 1)
    [ROI1_f_ROI1_prm_bin(:,iS), ROI1_prm_f_ROI1_prm_bin(:,iS)] = do_bin2(...
        ROI1_activity(:,iS), ROI1_prm(:,iS), nBins, 0);
    betas_ROI1(:,iS) = glmfit(ROI1_prm(:,iS), ROI1_activity(:,iS),'normal');
    ROI1_activity_fit(:,iS) = glmval(betas_ROI1(:,iS), ROI1_prm(:,iS), 'identity');
    ROI1_fit_f_ROI1_prm_bin(:,iS) = do_bin2(...
        ROI1_activity_fit(:,iS), ROI1_prm(:,iS), nBins, 0);

    % ROI 2 = f(prm 2)
    [ROI2_f_ROI2_prm_bin(:,iS), ROI2_prm_f_ROI2_prm_bin(:,iS)] = do_bin2(...
        ROI2_activity(:,iS), ROI2_prm(:,iS), nBins, 0);
    betas_ROI2(:,iS) = glmfit(ROI2_prm(:,iS), ROI2_activity(:,iS),'normal');
    ROI2_activity_fit(:,iS) = glmval(betas_ROI2(:,iS), ROI2_prm(:,iS), 'identity');
    ROI2_fit_f_ROI2_prm_bin(:,iS) = do_bin2(...
        ROI2_activity_fit(:,iS), ROI2_prm(:,iS), nBins, 0);

    % slope ROI 1=f(prm1) vs slope ROI2 = f(prm 2)
    [ROI2_f_ROI1_bin(:,iS), ROI1_f_ROI1_bin(:,iS)] = do_bin2(...
        ROI2_activity(:,iS), ROI1_activity(:,iS), nBins, 0);
    betas_ROI1_vs_ROI2_slopes(:,iS) = glmfit(betas_ROI1(:,iS), betas_ROI2(:,iS), 'normal');
    

    % ideally one should also be testing the serial mediation from dmPFC
    % activity to effort chosen to SMA activity during effort performance
end % subject loop

%% average
% ROI 1 = f(prm 1)
[m_ROI1_prm_f_ROI1_prm_bin,...
    sem_ROI1_prm_f_ROI1_prm_bin,...
    sd_ROI1_prm_f_ROI1_prm_bin] = mean_sem_sd(ROI1_prm_f_ROI1_prm_bin, 2);
[m_ROI1_f_ROI1_prm_bin,...
    sem_ROI1_f_ROI1_prm_bin,...
    sd_ROI1_f_ROI1_prm_bin] = mean_sem_sd(ROI1_f_ROI1_prm_bin, 2);
[m_betas.ROI1,...
    sem_betas.ROI1,...
    sd_betas.ROI1] = mean_sem_sd(betas_ROI1, 2);
% ROI 2 = f(prm 2)
[m_ROI2_prm_f_ROI2_prm_bin,...
    sem_ROI2_prm_f_ROI2_prm_bin,...
    sd_ROI2_prm_f_ROI2_prm_bin] = mean_sem_sd(ROI2_prm_f_ROI2_prm_bin, 2);
[m_ROI2_f_ROI2_prm_bin,...
    sem_ROI2_f_ROI2_prm_bin,...
    sd_ROI2_f_ROI2_prm_bin] = mean_sem_sd(ROI2_f_ROI2_prm_bin, 2);
[m_betas.ROI2,...
    sem_betas.ROI2,...
    sd_betas.ROI2] = mean_sem_sd(betas_ROI2, 2);
% ROI 2 vs ROI 1
[m_ROI2_f_ROI1_bin,...
    sem_ROI2_f_ROI1_bin,...
    sd_ROI2_f_ROI1_bin] = mean_sem_sd(ROI2_f_ROI1_bin, 2);
[m_ROI1_f_ROI1_bin,...
    sem_ROI1_f_ROI1_bin,...
    sd_ROI1_f_ROI1_bin] = mean_sem_sd(ROI1_f_ROI1_bin, 2);
[m_betas.ROI2_vs_ROI1,...
    sem_betas.ROI2_vs_ROI1,...
    sd_betas.ROI2_vs_ROI1] = mean_sem_sd(betas_ROI1_vs_ROI2_slopes, 2);

[m_betas.ROI2_vs_ROI1_bis,~,stats_tmp] = glmfit(betas_ROI1(2,:),...
    betas_ROI2(2,:),...
    'normal');
pval.([ROI_short_nm2,'_f_',ROI_short_nm1]) = stats_tmp.p;
betas_ROI2_vs_ROI1_bis_fit = glmval(m_betas.ROI2_vs_ROI1_bis, betas_ROI1(2,:), 'identity');

%% figure showing different bins
pSize = 50;
lWidth = 3;

%% start figure 1 where each ROI is compared to a parameter
fig1 = fig;

% ROI1 = f(prm 1)
subplot(1,3,1);
hdl = errorbar(m_ROI1_prm_f_ROI1_prm_bin,...
    m_ROI1_f_ROI1_prm_bin,...
    sem_ROI1_f_ROI1_prm_bin);
hdl.LineWidth = lWidth;
xlabel(ROI1_prm_nm);
ylabel(ROI_short_nm1);
legend_size(pSize);

% ROI2 = f(prm 2)
subplot(1,3,2);
hdl = errorbar(m_ROI2_prm_f_ROI2_prm_bin,...
    m_ROI2_f_ROI2_prm_bin,...
    sem_ROI2_f_ROI2_prm_bin);
hdl.LineWidth = lWidth;
legend_size(pSize);
xlabel(ROI2_prm_nm);
ylabel(ROI_short_nm2);

% ROI2 = f(ROI1) independent of behavioral parameters
subplot(1,3,3);
hdl = errorbar(m_ROI2_f_ROI1_bin,...
    m_ROI2_f_ROI1_bin,...
    sem_ROI2_f_ROI1_bin);
hdl.LineWidth = lWidth;
legend_size(pSize);
xlabel(ROI_short_nm1);
ylabel(ROI_short_nm2);

%% figure 2 check at the slope
fig2 = fig;
hold on;
hdl = scatter(betas_ROI1(2,:),...
    betas_ROI2(2,:));
hdl.MarkerEdgeColor = 'k';
hdl.LineWidth = lWidth;
[betas_ROI1_sorted, idx_b_ROI1] = sort(betas_ROI1(2,:));
betas_ROI2_vs_ROI1_bis_fit_sorted = betas_ROI2_vs_ROI1_bis_fit(idx_b_ROI1);
hdl_fit = plot(betas_ROI1_sorted,...
    betas_ROI2_vs_ROI1_bis_fit_sorted,...
    'LineStyle','--',...
    'LineWidth',lWidth,...
    'Color','k');
hdl_fit.LineWidth = lWidth;
xlabel([ROI_short_nm1,' b']);
ylabel([ROI_short_nm2,' b']);
legend_size(pSize);

