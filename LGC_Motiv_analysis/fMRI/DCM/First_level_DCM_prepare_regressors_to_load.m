function[] = First_level_DCM_prepare_regressors_to_load(study_nm, sub_nm, sub_idx, runs, n_runs,...
        subBehaviorFolder, computerRoot, n_scansPerRun, TR, DCM_mode)
%
% First_level_DCM_prepare_regressors_to_load will pool the regressors
% together depending on DCM_mode value so that all data are pooled either
% per run, per task or across everything.
%
% INPUTS
% DCM_mode:
% (1) all sessions modeled independently like in a classic univariate GLM
% => hard to manipulate for DCM but could be useful for testing
% session-specific effects or comparing sessions
% (2) sessions pooled within each task (ex: session 1 and 3 of physical
% effort will be concatenated into one single regressor) but each task will
% be modeled separately
% (3) all sessions pooled together
% (4) all trial periods are pooled together across sessions except for
% choice and effort which are modeled independently for each task (but
% pooled across sessions of the same task)
% (5) all trial periods are pooled together across sessions except for
% the effort period which is modeled independently for each task (but
% pooled across sessions of the same task)
%
% OUTPUTS
% all onsets, durations and regressors required for the 1st level

for iRun = 1:n_runs
    jRun = runs.runsToKeep(iRun);
    run_nm = num2str(jRun);
    task_id = runs.tasks{iRun};
    task_fullName = task_fullName_extraction(task_short_nm);
    %% determine task to work on
    switch task_fullName
        case 'physical'
            task_behavioral_id = 'physicalPerf';
        case 'mental'
            task_behavioral_id = 'mentalE_perf';
    end
    
    %% compute delay to add to all onsets to compensate for the fact that all sessions are now grouped
    if iRun == 1
        DCM_delay_onset = 0;
    else
        DCM_delay_onset = sum(n_scansPerRun(1:(iRun - 1))).*TR;
    end
    
    %% load data
    % warning('script needs to be updated to be able to model failed trials (choice or perf) separately');
    currRunBehaviorFileName = ls([subBehaviorFolder,'*_session',run_nm,'_*_task.mat']);
    if size(currRunBehaviorFileName,1) > 1
        error(['problem file identification: too many files popping out with run number',run_nm]);
    end
    if strcmp(currRunBehaviorFileName(16:23),'physical') ||...
            strcmp(currRunBehaviorFileName(17:24),'physical')
        task_nm = 'physical';
    elseif strcmp(currRunBehaviorFileName(16:21),'mental') ||...
            strcmp(currRunBehaviorFileName(17:22),'mental')
        task_nm = 'mental';
    else
        error('problem in identifying task type because file name doesn''t match');
    end
    if ~strcmp(task_id,task_nm)
        error('problem in identifying task type because file names don''t match between fMRI and behavior');
    end
    behavioralDataStruct = load([subBehaviorFolder, currRunBehaviorFileName]);
    
    %% extract all onsets of interest
    T0 = behavioralDataStruct.onsets.T0;
    if strcmp(study_nm,'fMRI_pilots') && ismember(subBehaviorFolder((end-30):end),...
            {[filesep,'fMRI_pilots',filesep,'pilot_s1',filesep,'behavior',filesep],...
            [filesep,'fMRI_pilots',filesep,'pilot_s2',filesep,'behavior',filesep]})
        % for pilots s1 & s2 where there was only 1 fixation cross before
        % choice, no cross before effort
        preChoiceCrossOnsets    = behavioralDataStruct.(task_behavioral_id).onsets.cross - T0 + DCM_delay_onset; % only cross appearing BEFORE the start of each trial
        preEffortCrossOnsets = [];
    else% for all other subjects
        preChoiceCrossOnsets    = behavioralDataStruct.(task_behavioral_id).onsets.preChoiceCross - T0 + DCM_delay_onset; % only cross appearing BEFORE the start of each trial
        preEffortCrossOnsets    = behavioralDataStruct.(task_behavioral_id).onsets.preEffortCross - T0 + DCM_delay_onset;
    end
    [allCrossesOnsets, allCrossesIdx] = sort([preChoiceCrossOnsets, preEffortCrossOnsets]);
    whiteCrossOnsets        = [preChoiceCrossOnsets, behavioralDataStruct.onsets.finalCross - T0 + DCM_delay_onset]; % add final cross (at the end of the experiment)
    dispChoiceOptionOnsets  = behavioralDataStruct.(task_behavioral_id).onsets.dispChoiceOptions - T0 + DCM_delay_onset;
    choiceOnsets            = behavioralDataStruct.(task_behavioral_id).onsets.choice - T0 + DCM_delay_onset;
    dispChosenOnsets        = behavioralDataStruct.(task_behavioral_id).onsets.dispChoice - T0 + DCM_delay_onset;
    n_trials = length(dispChosenOnsets);
    
    % extract trials where no choice was made
    choiceMissedTrials = isnan(choiceOnsets);
    
    % extract onsets when the effort was performed
    EperfOnsets = NaN(1,n_trials);
    for iTrial = 1:n_trials
        switch task_fullName
            case 'physical'
                EperfOnsets(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.effort_phase - T0 + DCM_delay_onset;
            case 'mental'
                EperfOnsets(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.nb_1 - T0 + DCM_delay_onset;
        end
    end
    % feedback onsets
    fbkOnsets = behavioralDataStruct.(task_behavioral_id).onsets.fbk - T0 + DCM_delay_onset;
    
    %% durations of each event
    if strcmp(study_nm,'fMRI_pilots') &&...
            ismember(subBehaviorFolder((end-30):end),...
            {[filesep,'fMRI_pilots',filesep,'pilot_s1',filesep,'behavior',filesep],...
            [filesep,'fMRI_pilots',filesep,'pilot_s2',filesep,'behavior',filesep],...
            [filesep,'fMRI_pilots',filesep,'pilot_s3',filesep,'behavior',filesep]})
        % duration not recorded yet in those pilots
        preChoiceCrossDur = dispChoiceOptionOnsets - preChoiceCrossOnsets;
        dispChoiceOptionsDur = dispChosenOnsets - dispChoiceOptionOnsets;
        if strcmp(study_nm,'fMRI_pilots') &&...
                ismember(subBehaviorFolder((end-30):end),...
                {[filesep,'fMRI_pilots',filesep,'pilot_s1',filesep,'behavior',filesep],...
                [filesep,'fMRI_pilots',filesep,'pilot_s2',filesep,'behavior',filesep]})
            % pilot s1 & s2 no fixation cross before effort
            dispChosenDur = EperfOnsets - dispChosenOnsets;
        else
            dispChosenDur = preEffortCrossOnsets - dispChosenOnsets;
            preEffortCrossDur = EperfOnsets - preEffortCrossOnsets;
        end
        EperfDur = fbkOnsets - EperfOnsets;
        fbkDur = whiteCrossOnsets(2:end) - fbkOnsets;
    else
        preChoiceCrossDur = behavioralDataStruct.(task_behavioral_id).durations.preChoiceCross;
        dispChoiceOptionsDur = behavioralDataStruct.(task_behavioral_id).durations.dispChoiceOptions;
        dispChosenDur = behavioralDataStruct.(task_behavioral_id).durations.dispChoice;
        preEffortCrossDur = behavioralDataStruct.(task_behavioral_id).durations.preEffortCross;
        EperfDur = behavioralDataStruct.(task_behavioral_id).durations.effortPeriod;
        fbkDur = behavioralDataStruct.(task_behavioral_id).durations.fbk;
    end
    allCrossesDur = [preChoiceCrossDur, preEffortCrossDur];
    allCrossesDur = allCrossesDur(allCrossesIdx);
    choice_RT = choiceOnsets - dispChoiceOptionOnsets;
    % choice_RT(isnan(choice_RT))=5;
    
    %% extract regressors of interest
    % loading R vs P trial
    RP_var_binary = strcmp( behavioralDataStruct.(task_behavioral_id).choiceOptions.R_or_P, 'R');
    RP_var = NaN(1, length(RP_var_binary));
    RP_var(RP_var_binary == 1) = 1;
    RP_var(RP_var_binary == 0) = -1;
    R_trials = RP_var == 1;
    P_trials = RP_var == -1;
    % loading default option side
    defaultSide = behavioralDataStruct.(task_behavioral_id).choiceOptions.default_LR;
    highE_on_the_right = defaultSide == -1;
    highE_on_the_left = defaultSide == 1;
    % loading money choice
    money_amount_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.left.*RP_var;
    money_amount_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.right.*RP_var;
    money_amount_sum = money_amount_left + money_amount_right;
    money_amount_varOption = money_amount_left.*highE_on_the_left + money_amount_right.*highE_on_the_right;
    abs_money_amount_varOption = abs(money_amount_varOption);
    money_level_left_v0 = behavioralDataStruct.(task_behavioral_id).choiceOptions.R.left;
    money_level_right_v0 = behavioralDataStruct.(task_behavioral_id).choiceOptions.R.right;
    money_amount_fixedOption = money_amount_left.*highE_on_the_right + money_amount_right.*highE_on_the_left;
    R_amount_varOption = money_amount_varOption.*R_trials;
    P_amount_varOption = money_amount_varOption.*P_trials;
    
    % fix money levels
    [money_level_left, money_level_right] = First_level_fix_money_levels(study_nm, sub_nm, task_fullName,...
        money_level_left_v0, money_level_right_v0, RP_var);
    
    % extract other relevant variables
    money_level_varOption = (money_level_left.*highE_on_the_left +...
        money_level_right.*highE_on_the_right) + 4*R_trials; % add 4 to distinguish rewards and punishments
    money_level_fixedOption = money_level_left.*highE_on_the_right +...
        money_level_right.*highE_on_the_left + 4*R_trials; % add 4 to distinguish rewards and punishments
    
    % consider that levels of reward (R1 to R4) are equivalent to levels of punishments (P1 to P4) and pool them together
    abs_money_level_varOption = (money_level_left.*highE_on_the_left +...
        money_level_right.*highE_on_the_right);
    
    R_level_varOption = (money_level_left.*highE_on_the_left +...
        money_level_right.*highE_on_the_right).*R_trials;
    P_level_varOption = (money_level_left.*highE_on_the_left +...
        money_level_right.*highE_on_the_right).*P_trials;
    % loading effort choice
    E_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.left;
    E_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.right;
    E_sum = E_left + E_right;
    E_varOption = E_left.*highE_on_the_left + E_right.*highE_on_the_right;
    E_fixedOption = E_left.*highE_on_the_right + E_right.*highE_on_the_left;
    % loading chosen option
    choice_LRandConf = behavioralDataStruct.(task_behavioral_id).choice;
    % transform into one variable equal to -1 when choice = left and 1 when
    % choice = right
    choice_LR = choice_LRandConf;
    choice_LR(choice_LRandConf == -2) = -1;
    choice_LR(choice_LRandConf == 2) = 1;
    choice_left = choice_LR == -1;
    choice_right = choice_LR == 1;
    % choice = high effort
    run_nm = num2str(jRun); % careful: use jRun for the run name
    [choice_hE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName); % (0) for low E chosen and (+1) for high E chosen
    choice_hE_bis = (choice_hE == 1) - (choice_hE == 0); % (-1) for low E chosen and (+1) for high E chosen
    % loading effort variables
    E_chosen = behavioralDataStruct.(task_behavioral_id).E_chosen; % 0 for low E chosen and 1/2/3 for high effort chosen
    E_chosen_bis = E_varOption.*choice_hE_bis; % -1/-2/-3 for low E chosen and +1/+2/+3 for high effort chosen
    E_unchosen = E_left.*choice_right + E_right.*choice_left;
    E_chosen_min_E_unchosen = E_chosen - E_unchosen;
    Ech_min_Efixed = E_chosen - E_fixedOption;
    
    % money levels
    money_level_chosen = money_level_left.*choice_left + money_level_right.*choice_right + 4*R_trials;% add 4 to distinguish rewards and punishments
    money_level_unchosen = money_level_left.*choice_right + money_level_right.*choice_left + 4*R_trials;% add 4 to distinguish rewards and punishments
    % consider that levels of reward (R1 to R4) are equivalent to levels of punishments (P1 to P4) and pool them together
    abs_money_level_chosen = money_level_left.*choice_left + money_level_right.*choice_right;
    abs_money_level_unchosen = money_level_left.*choice_right + money_level_right.*choice_left;
    moneyChosen_min_moneyUnchosen_level = money_level_chosen - money_level_unchosen;
    absMoneyChosen_min_moneyUnchosen_level = abs_money_level_chosen - abs_money_level_unchosen;
    % money chosen - money fixed option
    moneyChosen_min_moneyFixed_level = money_level_chosen - money_level_fixedOption;
    R_level_chosen = (money_level_left.*choice_left + money_level_right.*choice_right).*R_trials;
    P_level_chosen = (money_level_left.*choice_left + money_level_right.*choice_right).*P_trials;
    R_level_unchosen = (money_level_left.*choice_right + money_level_right.*choice_left).*R_trials;
    P_level_unchosen = (money_level_left.*choice_right + money_level_right.*choice_left).*P_trials;
    
    % amount variables
    money_amount_chosen = behavioralDataStruct.(task_behavioral_id).R_chosen.*RP_var;
    money_amount_unchosen = money_amount_left.*choice_right + money_amount_right.*choice_left;
    % money chosen - money unchosen option
    moneyChosen_min_moneyFixed_amount = money_amount_chosen - money_amount_fixedOption;
    moneyChosen_min_moneyUnchosen_amount = money_amount_chosen - money_amount_unchosen;
    absMoneyChosen_min_moneyUnchosen_amount = abs(money_amount_chosen - money_amount_unchosen);
    abs_money_amount_chosen = abs(money_amount_chosen);
    abs_money_amount_unchosen = abs(money_amount_unchosen);
    R_amount_chosen = money_amount_chosen.*R_trials;
    P_amount_chosen = money_amount_chosen.*P_trials;
    
    % interaction terms
    money_level_x_E_varOption = money_level_varOption.*E_varOption;
    money_level_x_E_chosen = money_level_chosen.*E_varOption;
    R_level_x_E_varOption = R_level_varOption.*E_varOption;
    R_level_x_E_chosen = R_level_chosen.*E_chosen;
    P_level_x_E_varOption = P_level_varOption.*E_varOption;
    P_level_x_E_chosen = P_level_chosen.*E_chosen;
    
    switch task_fullName
        case 'physical'
            [latency_tmp,...
                AUC_tmp, forcePeak_tmp, AUC_overshoot_tmp,...
                AUC_N_tmp, forcePeak_N_tmp, AUC_overshoot_N_tmp] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
            % time to start squeezing
            latency = latency_tmp.allTrials;
            % force variables in Volts
            AUC = AUC_tmp.allTrials;
            forcePeak = forcePeak_tmp.allTrials;
            AUC_overshoot = AUC_overshoot_tmp.allTrials;
            % force variables converted in Newtons
            AUC_N = AUC_N_tmp.allTrials;
            forcePeak_N = forcePeak_N_tmp.allTrials;
            AUC_overshoot_N = AUC_overshoot_N_tmp.allTrials;
            % extract fatigue
            fatigue = NaN(1,n_trials);
            for iTrial = 1:n_trials
                if iTrial == 1
                    fatigue(iTrial) = 0;
                else
                    fatigue(iTrial) = fatigue(iTrial-1) + AUC(iTrial-1);
                end
            end
        case 'mental'
            [~, n_correct_tmp, n_errors_tmp, RT_avg_tmp,...
                ~,...
                ~,...
                efficacy_with2first_tmp,...
                efficacy_pureNback_tmp,...
                efficacy_bis_with2first_tmp,...
                efficacy_bis_pureNback_tmp,...
                latency_tmp] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
            n_correct = n_correct_tmp.allTrials;
            n_errors = n_errors_tmp.allTrials;
            RT_avg = RT_avg_tmp.allTrials;
            efficacy_with2first = efficacy_with2first_tmp.allTrials;
            efficacy_pureNback = efficacy_pureNback_tmp.allTrials;
            efficacy_bis_with2first = efficacy_bis_with2first_tmp.allTrials;
            efficacy_bis_pureNback = efficacy_bis_pureNback_tmp.allTrials;
            [prevEfficacy_with2first,...
                prevEfficacy_pureNback,...
                prevEfficacy_bis_with2first,...
                prevEfficacy_bis_pureNback] = deal(NaN(1,n_trials));
            for iTrial = 1:n_trials
                if iTrial == 1
                    prevEfficacy_with2first(iTrial) = 0;
                    prevEfficacy_pureNback(iTrial) = 0;
                    prevEfficacy_bis_with2first(iTrial) = 0;
                    prevEfficacy_bis_pureNback(iTrial) = 0;
                else
                    prevEfficacy_with2first(iTrial) = efficacy_with2first(iTrial-1);
                    prevEfficacy_pureNback(iTrial) = efficacy_pureNback(iTrial-1);
                    prevEfficacy_bis_with2first(iTrial) = efficacy_bis_with2first(iTrial-1);
                    prevEfficacy_bis_pureNback(iTrial) = efficacy_bis_pureNback(iTrial-1);
                end
            end
            latency = latency_tmp.allTrials;
    end
    
    % loading feedback
    money_amount_obtained = behavioralDataStruct.(task_behavioral_id).gain; % could be different from reward chosen (in case of failure) but mostly similar
    win_vs_loss_fbk = money_amount_obtained > 0;
    % load run name
    switch jRun
        case {1,2}
            run_nm_bis = 'run1';
        case {3,4}
            run_nm_bis = 'run2';
        otherwise
            error('case not ready yet: net value in the model and more than 4 runs.');
    end
    % load net value and confidence
    RPconds = {'R','P','RP'};
    EsplitConditions = {'E','E1','E2','E3',...
        'Ech0','Ech1','Ech2','Ech3',...
        'lEch','hEch'};
    for iRP = 1:length(RPconds)
        RP_nm = RPconds{iRP};
        
        for iE = 1:length(EsplitConditions)
            Esplit_nm = EsplitConditions{iE};
            
            % load net value
            if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                    GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0 ||...
                    GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                    GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0 ||...
                    GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                    GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0
                
                % extract NV model name
                if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                        GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0
                    NV_mdl_nm = GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_mdl;
                elseif GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                        GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0
                    NV_mdl_nm = GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_mdl;
                elseif GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                        GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0
                    NV_mdl_nm = GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_mdl;
                end
                % extract net value
                if strcmp(NV_mdl_nm(1:4),'mdl_') % classic model
                    % load net value
                    [~, modelledDataStruct] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                        0, 'levels', 6, 6);
                    NV_ch_min_unch = modelledDataStruct.NV_chosen.(task_id).(NV_mdl_nm).(run_nm_bis);
                    NV_varOption = modelledDataStruct.NV_varOption.(task_id).(NV_mdl_nm).(run_nm_bis);
                    pChoice_hE = modelledDataStruct.pChoice_hE.(task_id).(NV_mdl_nm).(run_nm_bis);
                elseif strcmp(NV_mdl_nm(1:14),'bayesianModel_') % bayesian model
                    bayesianMdl_nm = strrep(NV_mdl_nm,'bayesianModel','mdl');
                    gitResultsFolder = [fullfile('C:','Users','clairis','Desktop',...
                        'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];
                    [NV_ch_min_unch, NV_varOption,~,pChoice_hE,...
                        NV_varOption_plus_bias, NV_ch_min_unch_with_bias] = extract_bayesian_mdl(gitResultsFolder, subBehaviorFolder,...
                        sub_nm, run_nm, task_fullName, bayesianMdl_nm);
                else
                    error(['model with ',NV_mdl_nm,' not ready yet']);
                end
                
                % probability of chosen option
                pChosen = pChoice_hE.*(choice_hE == 1)' + (1 - pChoice_hE).*(choice_hE == 0)';
            end % net value
            
            % load confidence
            if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence > 0 ||...
                    GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence > 0 ||...
                    GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence > 0 ||...
                    GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence > 0
                
                % use confidence based on ratings
                if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence == 1 ||...
                        GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence == 1 ||...
                        GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence == 1 ||...
                        GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence == 1
                    confidence = abs(choice_LRandConf) == 2; % 0 when low confidence and 1 when high confidence
                else
                    % extract NV model name
                    if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence > 1
                        conf_mdl_nm = GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).conf_mdl;
                    elseif GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence > 1
                        conf_mdl_nm = GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).conf_mdl;
                    elseif GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence > 1
                        conf_mdl_nm = GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).conf_mdl;
                    elseif GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence > 1
                        conf_mdl_nm = GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).conf_mdl;
                    end
                    % extract confidence
                    if strcmp(conf_mdl_nm(1:4),'mdl_') % classic model
                        % load inferred confidence
                        [~, modelledDataStruct] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                            0, 'levels', 6, 6);
                        confidence_highE = modelledDataStruct.confidenceFitted_highE.(conf_mdl_nm).(task_id).(run_nm_bis);
                        pChoice = modelledDataStruct.choicesFitted.(conf_mdl_nm).(task_id).(run_nm_bis);
                    elseif strcmp(conf_mdl_nm(1:14),'bayesianModel_') % bayesian model
                        bayesianMdl_nm = strrep(conf_mdl_nm,'bayesianModel','mdl');
                        gitResultsFolder = [fullfile('C:','Users','clairis','Desktop',...
                            'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];
                        [~, ~, confidence_highE, pChoice] = extract_bayesian_mdl(gitResultsFolder, subBehaviorFolder,...
                            sub_nm, run_nm, task_fullName, bayesianMdl_nm);
                    else
                        error(['model with ',conf_mdl_nm,' not ready yet']);
                    end
                    
                    % define confidence proxy to use
                    if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence == 2 ||...
                            GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence == 2 ||...
                            GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence == 2 ||...
                            GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence == 2
                        confidence = confidence_highE;
                    elseif GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence == 3 ||...
                            GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence == 3 ||...
                            GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence == 3 ||...
                            GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence == 3
                        confidence_left = (pChoice.*highE_on_the_left + (1 - pChoice).*highE_on_the_right - 0.5).^2;
                        confidence = confidence_left./0.25; % normalise to be between 0 and 1
                    elseif GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence == 4 ||...
                            GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence == 4 ||...
                            GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence == 4 ||...
                            GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence == 4
                        confidence = ((pChoice - 0.5).^2)./0.25;
                    end
                    
                end % which confidence to use
            end % confidence
        end % Effort conditions
    end % R/P/RP loop
    
    % extract fatigue
    n_theoreticalTrialsPerRun = 54;
    trialN = (1:n_theoreticalTrialsPerRun) - 1; % trial number - 1 to reflect fatigue level
    trialN_dEch = trialN.*E_chosen_min_E_unchosen;
    trialN_dEnonDef_min_Edef = trialN.*(E_varOption - E_fixedOption);
    trialN_dEnonDef = trialN.*E_varOption;
    
    %% remove trials where no choice was performed
    if sum(choiceMissedTrials) > 0
        % note: allCross not cleaned from missed choices on purpose as it would
        % make no sense
        
        % extract onsets and durations of missed trials
        preChoiceCrossOnsets_missedOnsets = preChoiceCrossOnsets(choiceMissedTrials);
        dispChoiceOption_missedOnsets = dispChoiceOptionOnsets(choiceMissedTrials);
        choice_missedOnsets = choiceOnsets(choiceMissedTrials);
        dispChosen_missedOnsets = dispChosenOnsets(choiceMissedTrials);
        preEffortCross_missedOnsets = preEffortCrossOnsets(choiceMissedTrials);
        Eperf_missedOnsets = EperfOnsets(choiceMissedTrials);
        fbk_missedOnsets = fbkOnsets(choiceMissedTrials);
        % durations
        preChoiceCross_missedDur = preChoiceCrossDur(choiceMissedTrials);
        dispChoiceOptions_missedDur = dispChoiceOptionsDur(choiceMissedTrials);
        dispChosen_missedDur = dispChosenDur(choiceMissedTrials);
        preEffortCross_missedDur = preEffortCrossDur(choiceMissedTrials);
        Eperf_missedDur = EperfDur(choiceMissedTrials);
        fbk_missedDur = fbkDur(choiceMissedTrials);
        
        % onsets
        preChoiceCrossOnsets(choiceMissedTrials) = [];
        dispChoiceOptionOnsets(choiceMissedTrials) = [];
        choiceOnsets(choiceMissedTrials) = [];
        dispChosenOnsets(choiceMissedTrials) = [];
        preEffortCrossOnsets(choiceMissedTrials) = [];
        EperfOnsets(choiceMissedTrials) = [];
        fbkOnsets(choiceMissedTrials) = [];
        % durations
        preChoiceCrossDur(choiceMissedTrials) = [];
        dispChoiceOptionsDur(choiceMissedTrials) = [];
        dispChosenDur(choiceMissedTrials) = [];
        preEffortCrossDur(choiceMissedTrials) = [];
        EperfDur(choiceMissedTrials) = [];
        fbkDur(choiceMissedTrials) = [];
        % regressors
        R_trials(choiceMissedTrials) = [];
        P_trials(choiceMissedTrials) = [];
        defaultSide(choiceMissedTrials) = [];
        RP_var_binary(choiceMissedTrials) = [];
        choice_hE(choiceMissedTrials) = [];
        choice_hE_bis(choiceMissedTrials) = [];
        money_level_chosen(choiceMissedTrials) = [];
        R_level_chosen(choiceMissedTrials) = [];
        P_level_chosen(choiceMissedTrials) = [];
        R_level_unchosen(choiceMissedTrials) = [];
        P_level_unchosen(choiceMissedTrials) = [];
        money_amount_left(choiceMissedTrials) = [];
        money_amount_right(choiceMissedTrials) = [];
        money_amount_sum(choiceMissedTrials) = [];
        money_amount_varOption(choiceMissedTrials) = [];
        R_amount_varOption(choiceMissedTrials) = [];
        P_amount_varOption(choiceMissedTrials) = [];
        money_amount_fixedOption(choiceMissedTrials) = [];
        money_level_left(choiceMissedTrials) = [];
        money_level_right(choiceMissedTrials) = [];
        money_level_varOption(choiceMissedTrials) = [];
        R_level_varOption(choiceMissedTrials) = [];
        P_level_varOption(choiceMissedTrials) = [];
        abs_money_amount_varOption(choiceMissedTrials) = [];
        abs_money_level_varOption(choiceMissedTrials) = [];
        E_left(choiceMissedTrials) = [];
        E_right(choiceMissedTrials) = [];
        E_sum(choiceMissedTrials) = [];
        E_varOption(choiceMissedTrials) = [];
        money_level_x_E_varOption(choiceMissedTrials) = [];
        money_level_x_E_chosen(choiceMissedTrials) = [];
        R_level_x_E_varOption(choiceMissedTrials) = [];
        R_level_x_E_chosen(choiceMissedTrials) = [];
        P_level_x_E_varOption(choiceMissedTrials) = [];
        P_level_x_E_chosen(choiceMissedTrials) = [];
        choice_LRandConf(choiceMissedTrials) = [];
        choice_LR(choiceMissedTrials) = [];
        if exist('confidence','var') && ~isempty(confidence)
            confidence(choiceMissedTrials) = [];
        end
        E_chosen(choiceMissedTrials) = [];
        E_chosen_bis(choiceMissedTrials) = [];
        E_unchosen(choiceMissedTrials) = [];
        E_chosen_min_E_unchosen(choiceMissedTrials) = [];
        Ech_min_Efixed(choiceMissedTrials) = [];
        money_amount_chosen(choiceMissedTrials) = [];
        R_amount_chosen(choiceMissedTrials) = [];
        P_amount_chosen(choiceMissedTrials) = [];
        abs_money_amount_chosen(choiceMissedTrials) = [];
        abs_money_level_chosen(choiceMissedTrials) = [];
        abs_money_amount_unchosen(choiceMissedTrials) = [];
        abs_money_level_unchosen(choiceMissedTrials) = [];
        money_amount_unchosen(choiceMissedTrials) = [];
        money_level_unchosen(choiceMissedTrials) = [];
        moneyChosen_min_moneyUnchosen_amount(choiceMissedTrials) = [];
        absMoneyChosen_min_moneyUnchosen_amount(choiceMissedTrials) = [];
        moneyChosen_min_moneyUnchosen_level(choiceMissedTrials) = [];
        absMoneyChosen_min_moneyUnchosen_level(choiceMissedTrials) = [];
        moneyChosen_min_moneyFixed_amount(choiceMissedTrials) = [];
        moneyChosen_min_moneyFixed_level(choiceMissedTrials) = [];
        money_amount_obtained(choiceMissedTrials) = [];
        win_vs_loss_fbk(choiceMissedTrials) = [];
        choice_RT(choiceMissedTrials) = [];
        if exist('NV_mdl_nm','var') &&...
                (strcmp(NV_mdl_nm(1:4),'mdl_') ||...
                strcmp(NV_mdl_nm(1:14),'bayesianModel_'))
            NV_ch_min_unch(choiceMissedTrials) = [];
            NV_ch_min_unch_with_bias(choiceMissedTrials) = [];
            NV_varOption(choiceMissedTrials) = [];
            NV_varOption_plus_bias(choiceMissedTrials) = [];
            pChoice_hE(choiceMissedTrials) = [];
            pChosen(choiceMissedTrials) = [];
        end
        trialN(choiceMissedTrials) = [];
        trialN_dEch(choiceMissedTrials) = [];
        trialN_dEnonDef_min_Edef(choiceMissedTrials) = [];
        trialN_dEnonDef(choiceMissedTrials) = [];
        switch task_fullName
            case 'physical'
                AUC(choiceMissedTrials) = [];
                forcePeak(choiceMissedTrials) = [];
                AUC_overshoot(choiceMissedTrials) = [];
                AUC_N(choiceMissedTrials) = [];
                forcePeak_N(choiceMissedTrials) = [];
                AUC_overshoot_N(choiceMissedTrials) = [];
                fatigue(choiceMissedTrials) = [];
            case 'mental'
                n_correct(choiceMissedTrials) = [];
                n_errors(choiceMissedTrials) = [];
                RT_avg(choiceMissedTrials) = [];
                efficacy_with2first(choiceMissedTrials) = [];
                efficacy_pureNback(choiceMissedTrials) = [];
                efficacy_bis_with2first(choiceMissedTrials) = [];
                efficacy_bis_pureNback(choiceMissedTrials) = [];
                prevEfficacy_with2first(choiceMissedTrials) = [];
                prevEfficacy_pureNback(choiceMissedTrials) = [];
                prevEfficacy_bis_with2first(choiceMissedTrials) = [];
                prevEfficacy_bis_pureNback(choiceMissedTrials) = [];
        end
        latency(choiceMissedTrials) = [];
    end
    
    % zscored R/P ignoring P/R trials (respectively)
    [z_R_level_chosen, z_P_level_chosen,...
        z_R_level_unchosen, z_P_level_unchosen,...
        z_R_amount_chosen, z_P_amount_chosen,...
        z_R_level_varOption, z_P_level_varOption,...
        z_R_amount_varOption, z_P_amount_varOption] = deal(zeros(size(R_level_chosen)));
    z_R_level_chosen(R_trials) = zscore(R_level_chosen(R_trials));
    z_P_level_chosen(P_trials) = zscore(P_level_chosen(P_trials));
    z_R_level_unchosen(R_trials) = zscore(R_level_unchosen(R_trials));
    z_P_level_unchosen(P_trials) = zscore(P_level_unchosen(P_trials));
    z_R_amount_chosen(R_trials) = zscore(R_amount_chosen(R_trials));
    z_P_amount_chosen(P_trials) = zscore(P_amount_chosen(P_trials));
    z_R_level_varOption(R_trials) = zscore(R_level_varOption(R_trials));
    z_P_level_varOption(P_trials) = zscore(P_level_varOption(P_trials));
    z_R_amount_varOption(R_trials) = zscore(R_amount_varOption(R_trials));
    z_P_amount_varOption(P_trials) = zscore(P_amount_varOption(P_trials));
    
    % extract trials where no performance was achieved
    perfOkTrials = ~isnan(latency);

    %% pool the data depending on DCM_mode
    switch DCM_mode
        case 1 % no change to be done

        case 2 % pooling sessions/task but each task is modeled independently

        case 3 % pooling all sessions

        case 4 % pooling all sessions except for effort and choice

        case 5 % pooling all sessions except for effort period

    end % DCM_mode

end % run loop to prepare onsets and regressors of interest

end % function