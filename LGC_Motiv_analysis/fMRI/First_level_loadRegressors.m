function[matlabbatch] = First_level_loadRegressors(matlabbatch, GLMprm, study_nm, sub_idx, iRun,...
    subj_behavior_folder, currRunBehaviorFileName, task_nm)
% [matlabbatch] = First_level_loadRegressors(matlabbatch, GLMprm, study_nm, sub_idx, iRun,...
%     subj_behavior_folder, currRunBehaviorFileName, task_nm)
%
% First_level_loadRegressors will load the regressors of interest for each
% task
%
% INPUTS
% matlabbatch: structure with the First level batch
%
% GLMprm: structure with GLM parameters
%
% study_nm: study name
%
% sub_idx: index for matlabbatch
%
% iRun: run number
%
% subj_behavior_folder: subject behavioral folder name
%
% currRunBehaviorFilename: name of the file containing the behavioral data
%
% task_nm: task name 'mental' or 'physical'
%
% OUTPUTS
% matlabbatch: structure updated with the regressors of interest depending
% on GLMprm values

%% determine task to work on
switch task_nm
    case 'physical'
        task_id = 'Ep';
        task_behavioral_id = 'physicalPerf';
    case 'mental'
        task_id = 'Em';
        task_behavioral_id = 'mentalE_perf';
end

%% load data
warning('script needs to be updated to take failed trials into account and modelled separately');
behavioralDataStruct = load([subj_behavior_folder, currRunBehaviorFileName]);

%% extract all onsets of interest
T0 = behavioralDataStruct.onsets.T0;
if strcmp(study_nm,'fMRI_pilots') && ismember(subj_behavior_folder((end-30):end),...
        {[filesep,'fMRI_pilots',filesep,'pilot_s1',filesep,'behavior',filesep],...
        [filesep,'fMRI_pilots',filesep,'pilot_s2',filesep,'behavior',filesep]})
    % for pilots s1 & s2 where there was only 1 fixation cross before
    % choice, no cross before effort
    preChoiceCrossOnsets    = behavioralDataStruct.(task_behavioral_id).onsets.cross - T0; % only cross appearing BEFORE the start of each trial
else% for all other subjects
    preChoiceCrossOnsets    = behavioralDataStruct.(task_behavioral_id).onsets.preChoiceCross - T0; % only cross appearing BEFORE the start of each trial
    preEffortCrossOnsets    = behavioralDataStruct.(task_behavioral_id).onsets.preEffortCross - T0;
end
whiteCrossOnsets        = [preChoiceCrossOnsets, behavioralDataStruct.onsets.finalCross - T0]; % add final cross (at the end of the experiment)
dispChoiceOptionOnsets  = behavioralDataStruct.(task_behavioral_id).onsets.dispChoiceOptions - T0;
choiceOnsets            = behavioralDataStruct.(task_behavioral_id).onsets.choice - T0;
dispChosenOnsets        = behavioralDataStruct.(task_behavioral_id).onsets.dispChoice - T0;
n_trials = length(dispChosenOnsets);

% extract trials where no choice was made
choiceMissedTrials = isnan(choiceOnsets);

% extract onsets when the effort was performed
EperfOnsets = NaN(1,n_trials);
for iTrial = 1:n_trials
    switch task_nm
        case 'physical'
            EperfOnsets(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.effort_phase - T0;
        case 'mental'
            EperfOnsets(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.nb_1 - T0;
    end
end
% feedback onsets
fbkOnsets = behavioralDataStruct.(task_behavioral_id).onsets.fbk - T0;

%% durations of each event
if strcmp(study_nm,'fMRI_pilots') &&...
        ismember(subj_behavior_folder((end-30):end),...
        {[filesep,'fMRI_pilots',filesep,'pilot_s1',filesep,'behavior',filesep],...
        [filesep,'fMRI_pilots',filesep,'pilot_s2',filesep,'behavior',filesep],...
        [filesep,'fMRI_pilots',filesep,'pilot_s3',filesep,'behavior',filesep]})
    % duration not recorded yet in those pilots
    preChoiceCrossDur = dispChoiceOptionOnsets - preChoiceCrossOnsets;
    dispChoiceOptionsDur = dispChosenOnsets - dispChoiceOptionOnsets;
    if strcmp(study_nm,'fMRI_pilots') &&...
            ismember(subj_behavior_folder((end-30):end),...
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
choice_RT = choiceOnsets - dispChoiceOptionOnsets;

%% extract regressors of interest
% loading R vs P trial
RP_var_binary = strcmp( behavioralDataStruct.(task_behavioral_id).choiceOptions.R_or_P, 'R');
RP_var = RP_var_binary;
RP_var(RP_var_binary == 0) = -1;
% loading default option side
defaultSide = behavioralDataStruct.(task_behavioral_id).choiceOptions.default_LR;
% loading money choice
money_amount_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.left.*RP_var;
money_amount_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.right.*RP_var;
money_amount_sum = money_amount_left + money_amount_right;
money_amount_varOption = money_amount_left.*(defaultSide == 1) + money_amount_right.*(defaultSide == -1);
abs_money_amount_varOption = abs(money_amount_varOption);
money_level_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.R.left.*RP_var;
money_level_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.R.right.*RP_var;
money_level_varOption = money_level_left.*(defaultSide == 1) + money_level_right.*(defaultSide == -1);
abs_money_level_varOption = abs(money_level_varOption);
% loading effort choice
E_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.left.*RP_var;
E_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.left.*RP_var;
E_sum = E_left + E_right;
E_varOption = E_left.*(defaultSide == 1) + E_right.*(defaultSide == -1);
% loading chosen option
choice_LRandConf = behavioralDataStruct.(task_behavioral_id).choice;
% transform into one variable equal to -1 when choice = left and 1 when
% choice = right
choice_LR = choice_LRandConf;
choice_LR(choice_LRandConf == -2) = -1;
choice_LR(choice_LRandConf == 2) = 1;
confidence = abs(choice_LRandConf) == 2; % 0 when low confidence and 1 when high confidence
% loading chosen variables
E_chosen = behavioralDataStruct.(task_behavioral_id).E_chosen;
E_unchosen = E_left.*(choice_LR == 1) + E_right.*(choice_LR == -1);
E_chosen_min_E_unchosen = E_chosen - E_unchosen;
money_amount_chosen = behavioralDataStruct.(task_behavioral_id).R_chosen;
money_level_chosen = money_level_left.*(choice_LR == -1) + money_level_right.*(choice_LR == 1);
money_level_unchosen = money_level_left.*(choice_LR == 1) + money_level_right.*(choice_LR == -1);
abs_money_amount_chosen = abs(money_amount_chosen);
abs_money_level_chosen = abs(money_level_chosen);
money_amount_unchosen = money_amount_left.*(choice_LR == 1) + money_amount_right.*(choice_LR == -1);
abs_money_amount_unchosen = abs(money_amount_unchosen);
abs_money_level_unchosen = abs(money_level_unchosen);
moneyChosen_min_moneyUnchosen_amount = money_amount_chosen - money_amount_unchosen;
absMoneyChosen_min_moneyUnchosen_amount = abs(money_amount_chosen - money_amount_unchosen);
moneyChosen_min_moneyUnchosen_level = money_level_chosen - money_level_unchosen;
absMoneyChosen_min_moneyUnchosen_level = abs(money_level_chosen - money_level_unchosen);
% loading feedback
money_amount_obtained = behavioralDataStruct.(task_behavioral_id).gain; % could be different from reward chosen (in case of failure) but mostly similar
win_vs_loss_fbk = money_amount_obtained > 0;

%% remove trials where no choice was performed
if sum(choiceMissedTrials) > 0
    % extract onsets and durations of missed trials
    dispChoiceOption_missedOnsets = dispChoiceOptionOnsets(choiceMissedTrials);
    choice_missedOnsets = choiceOnsets(choiceMissedTrials);
    dispChosen_missedOnsets = dispChosenOnsets(choiceMissedTrials);
    Eperf_missedOnsets = EperfOnsets(choiceMissedTrials);
    fbk_missedOnsets = fbkOnsets(choiceMissedTrials);
    % durations
    dispChoiceOptions_missedDur = dispChoiceOptionsDur(choiceMissedTrials);
    dispChosen_missedDur = dispChosenDur(choiceMissedTrials);
    Eperf_missedDur = EperfDur(choiceMissedTrials);
    fbk_missedDur = fbkDur(choiceMissedTrials);
    
    % onsets
    dispChoiceOptionOnsets(choiceMissedTrials) = [];
    choiceOnsets(choiceMissedTrials) = [];
    dispChosenOnsets(choiceMissedTrials) = [];
    EperfOnsets(choiceMissedTrials) = [];
    fbkOnsets(choiceMissedTrials) = [];
    % durations
    dispChoiceOptionsDur(choiceMissedTrials) = [];
    dispChosenDur(choiceMissedTrials) = [];
    EperfDur(choiceMissedTrials) = [];
    fbkDur(choiceMissedTrials) = [];
    % regressors
    RP_var_binary(choiceMissedTrials) = [];
    money_amount_left(choiceMissedTrials) = [];
    money_amount_right(choiceMissedTrials) = [];
    money_amount_sum(choiceMissedTrials) = [];
    money_amount_varOption(choiceMissedTrials) = [];
    money_level_left(choiceMissedTrials) = [];
    money_level_right(choiceMissedTrials) = [];
    money_level_varOption(choiceMissedTrials) = [];
    abs_money_amount_varOption(choiceMissedTrials) = [];
    abs_money_level_varOption(choiceMissedTrials) = [];
    E_left(choiceMissedTrials) = [];
    E_right(choiceMissedTrials) = [];
    E_sum(choiceMissedTrials) = [];
    E_varOption(choiceMissedTrials) = [];
    choice_LRandConf(choiceMissedTrials) = [];
    choice_LR(choiceMissedTrials) = [];
    confidence(choiceMissedTrials) = [];
    E_chosen(choiceMissedTrials) = [];
    E_unchosen(choiceMissedTrials) = [];
    E_chosen_min_E_unchosen(choiceMissedTrials) = [];
    money_amount_chosen(choiceMissedTrials) = [];
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
    money_amount_obtained(choiceMissedTrials) = [];
    win_vs_loss_fbk(choiceMissedTrials) = [];
    choice_RT(choiceMissedTrials) = [];
end

%% load the batch according to GLMprm variables
% load general GLM parameters
orth_vars = GLMprm.gal.orth_vars;
zPerRun = GLMprm.gal.zPerRun;
switch zPerRun
    case 0 % raw data = no change
        raw_or_z = @(x) x;
    case 1 % zscore all variables per run
        raw_or_z = @(x) zscore(x);
end

% GLM parameters for choice
choice_RPpool = GLMprm.choice.(task_id).RPpool;
switch choice_RPpool
    case 0
        RPchoiceCond = {'R','P'};
    case 1
        RPchoiceCond = {'RP'};
end
chosen_RPpool = GLMprm.chosen.(task_id).RPpool;
switch chosen_RPpool
    case 0
        RPchosenCond = {'R','P'};
    case 1
        RPchosenCond = {'RP'};
end
Eperf_RPpool = GLMprm.Eperf.(task_id).RPpool;
switch Eperf_RPpool
    case 0
        RPperfCond = {'R','P'};
    case 1
        RPperfCond = {'RP'};
end
fbk_RPpool = GLMprm.fbk.(task_id).RPpool;
switch fbk_RPpool
    case 0
        RPfbkCond = {'R','P'};
    case 1
        RPfbkCond = {'RP'};
end

%% initialize conditions
iCond = 0;

%% fixation cross before choice
preChoiceCrossModel = GLMprm.model_onset.(task_id).preChoiceCross;
if ismember(preChoiceCrossModel,{'stick','boxcar'})
    iCond = iCond + 1;
    switch preChoiceCrossModel
        case 'stick'
            modelPreChoiceCrossdur = 0;
        case 'boxcar'
            modelPreChoiceCrossdur = preChoiceCrossDur;
    end
    [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
        'preChoice fixation cross', whiteCrossOnsets, modelPreChoiceCrossdur, 0, '', [], orth_vars);
end


%% choice period
choiceModel = GLMprm.model_onset.(task_id).choice;
if ismember(choiceModel,{'stick','boxcar'})

    for iRP_choice = 1:length(RPchoiceCond)
        RP_choice_nm                                = RPchoiceCond{iRP_choice};
        choiceModel_RP                              = GLMprm.choice.(task_id).(RP_choice_nm).R_vs_P;
        choiceModel_moneyLeft                       = GLMprm.choice.(task_id).(RP_choice_nm).money_left;
        choiceModel_moneyRight                      = GLMprm.choice.(task_id).(RP_choice_nm).money_left;
        choiceModel_moneySum                        = GLMprm.choice.(task_id).(RP_choice_nm).money_sum;
        choiceModel_E_sum                           = GLMprm.choice.(task_id).(RP_choice_nm).E_sum;
        choiceModel_moneyChosen                     = GLMprm.choice.(task_id).(RP_choice_nm).money_chosen;
        choiceModel_moneyUnchosen                   = GLMprm.choice.(task_id).(RP_choice_nm).money_unchosen;
        choiceModel_moneyNonDefault                 = GLMprm.choice.(task_id).(RP_choice_nm).money_varOption;
        choiceModel_money_chosen_min_money_unchosen = GLMprm.choice.(task_id).(RP_choice_nm).money_ch_min_unch;
        choiceModel_E_chosen                        = GLMprm.choice.(task_id).(RP_choice_nm).E_chosen;
        choiceModel_E_unchosen                      = GLMprm.choice.(task_id).(RP_choice_nm).E_unchosen;
        choiceModel_E_left                          = GLMprm.choice.(task_id).(RP_choice_nm).E_left;
        choiceModel_E_right                         = GLMprm.choice.(task_id).(RP_choice_nm).E_right;
        choiceModel_E_nonDefault                    = GLMprm.choice.(task_id).(RP_choice_nm).E_varOption;
        choiceModel_E_chosen_min_E_unchosen         = GLMprm.choice.(task_id).(RP_choice_nm).E_ch_min_unch;
        choiceModel_RT                              = GLMprm.choice.(task_id).(RP_choice_nm).RT;
        choiceModel_conf                            = GLMprm.choice.(task_id).(RP_choice_nm).confidence;

        % extract trial index for the current loop
        switch RP_choice_nm
            case 'RP'
                choice_trial_idx = 1:length(RP_var_binary);
            case 'R'
                choice_trial_idx = RP_var_binary == 1;
            case 'P'
                choice_trial_idx = RP_var_binary == 0;
        end
        
        % onset
        iCond = iCond + 1;
        modelChoiceOnset = dispChoiceOptionOnsets(choice_trial_idx);
        % duration
        switch choiceModel
            case 'stick'
                modelChoiceDur = 0;
            case 'boxcar'
                modelChoiceDur = dispChoiceOptionsDur(choice_trial_idx);
        end

        % modulators
        n_choiceMods = 0;
        choice_modNames = cell(1,1);
        choice_modVals = [];

        % reward vs punishments
        if choiceModel_RP == 1
            if strcmp(RP_choice_nm,'RP')
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'R vs P';
                choice_modVals(n_choiceMods,:) = RP_var_binary; % binary variable => no zscore
            else
                error('cannot split R and P trials and add a variable representing R/P trials') ;
            end
        end

        % money left
        if choiceModel_moneyLeft > 0
            error('not ready yet');
        end

        % money right
        if choiceModel_moneyRight > 0
            error('not ready yet');
        end

        % money chosen
        if choiceModel_moneyChosen > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'money chosen';
            switch choiceModel_moneyChosen
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_chosen(choice_trial_idx));
                case 2
                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_chosen(choice_trial_idx));
                case 3
                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_chosen(choice_trial_idx));
                case 4
                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_chosen(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % money unchosen
        if choiceModel_moneyUnchosen > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'money unchosen';
            switch choiceModel_moneyUnchosen
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_unchosen(choice_trial_idx));
                case 2
                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_unchosen(choice_trial_idx));
                case 3
                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_unchosen(choice_trial_idx));
                case 4
                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_unchosen(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % money chosen - money unchosen
        if choiceModel_money_chosen_min_money_unchosen > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'money chosen - unchosen';
            switch choiceModel_money_chosen_min_money_unchosen
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen_amount(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % sum of money
        if choiceModel_moneySum > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'money sum';
            switch choiceModel_moneySum
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_sum(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % money associated to the non-default option
        if choiceModel_moneyNonDefault > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'money non-default';
            switch choiceModel_moneyNonDefault
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_varOption(choice_trial_idx));
                case 2
                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_varOption(choice_trial_idx));
                case 3
                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_varOption(choice_trial_idx));
                case 4
                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_varOption(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % effort left
        if choiceModel_E_left > 0
            error('not ready yet');
        end

        % effort right
        if choiceModel_E_right > 0
            error('not ready yet');
        end

        % effort chosen
        if choiceModel_E_chosen > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'effort chosen';
            switch choiceModel_E_chosen
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % effort unchosen
        if choiceModel_E_unchosen > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'effort unchosen';
            switch choiceModel_E_unchosen
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(E_unchosen(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % effort non-default option
        if choiceModel_E_nonDefault > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'effort non-default option';
            switch choiceModel_E_nonDefault
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(E_varOption(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % effort chosen - unchosen
        if choiceModel_E_chosen_min_E_unchosen > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'effort chosen - unchosen';
            switch choiceModel_E_chosen_min_E_unchosen
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_min_E_unchosen(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % sum of efforts
        if choiceModel_E_sum > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'effort sum';
            switch choiceModel_E_sum
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(E_sum(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % choice confidence
        if choiceModel_conf > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'confidence';
            switch choiceModel_conf % binary variable => no zscore
                case 1
                    choice_modVals(n_choiceMods,:) = confidence(choice_trial_idx);
                otherwise
                    error('not ready yet');
            end
        end

        % RT
        if choiceModel_RT > 0
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'choice RT';
            switch choiceModel_RT
                case 1
                    choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT(choice_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
            ['choice_',RP_choice_nm], modelChoiceOnset, modelChoiceDur,...
            n_choiceMods, choice_modNames, choice_modVals,...
            orth_vars);
    end % RP
end % model choice

%% chosen period
chosenModel = GLMprm.model_onset.(task_id).chosen;
if ismember(chosenModel,{'stick','boxcar'})

    for iRP_chosen = 1:length(RPchosenCond)
        RP_chosen_nm = RPchosenCond{iRP_chosen};

        chosenModel_moneyChosen = GLMprm.chosen.(task_id).(RP_chosen_nm).money_chosen;
        chosenModel_moneyNonDefault = GLMprm.chosen.(task_id).(RP_chosen_nm).money_varOption;
        chosenModel_Echosen = GLMprm.chosen.(task_id).(RP_chosen_nm).E_chosen;
        chosenModel_EnonDefault = GLMprm.chosen.(task_id).(RP_chosen_nm).E_varOption;
        chosenModel_moneyUnchosen = GLMprm.chosen.(task_id).(RP_chosen_nm).money_unchosen;
        chosenModel_Eunchosen = GLMprm.chosen.(task_id).(RP_chosen_nm).E_unchosen;
        chosenModel_confidence = GLMprm.chosen.(task_id).(RP_chosen_nm).confidence;

        % extract trial index for the current loop
        switch RP_chosen_nm
            case 'RP'
                chosen_trial_idx = 1:length(RP_var_binary);
            case 'R'
                chosen_trial_idx = RP_var_binary == 1;
            case 'P'
                chosen_trial_idx = RP_var_binary == 0;
        end
        
        % onset
        iCond = iCond + 1;
        modelChosenOnset = dispChosenOnsets(chosen_trial_idx);
        % duration
        switch chosenModel
            case 'stick'
                modelChosenDur = 0;
            case 'boxcar'
                modelChosenDur = dispChosenDur(chosen_trial_idx);
        end

        % modulators
        n_chosenMods = 0;
        chosen_modNames = cell(1,1);
        chosen_modVals = [];

        % money chosen
        if chosenModel_moneyChosen > 0
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'money chosen';
            switch chosenModel_moneyChosen
                case 1
                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_chosen(chosen_trial_idx));
                case 2
                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_chosen(chosen_trial_idx));
                case 3
                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_chosen(chosen_trial_idx));
                case 4
                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_chosen(chosen_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % money non-default option
        if chosenModel_moneyNonDefault > 0
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'money non-default option';
            switch chosenModel_moneyNonDefault
                case 1
                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_varOption(chosen_trial_idx));
                case 2
                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_varOption(chosen_trial_idx));
                case 3
                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_varOption(chosen_trial_idx));
                case 4
                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_varOption(chosen_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % effort chosen
        if chosenModel_Echosen > 0
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'effort chosen';
            switch chosenModel_Echosen
                case 1
                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen(chosen_trial_idx));
                otherwise
                    error('ready yet');
            end
        end

        % money unchosen
        if chosenModel_moneyUnchosen > 0
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'money unchosen';
            switch chosenModel_moneyUnchosen
                case 1
                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_unchosen(chosen_trial_idx));
                case 2
                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_unchosen(chosen_trial_idx));
                case 3
                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_unchosen(chosen_trial_idx));
                case 4
                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_unchosen(chosen_trial_idx));
                otherwise
                    error('ready yet');
            end
        end

        % effort unchosen
        if chosenModel_Eunchosen > 0
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'effort unchosen';
            switch chosenModel_Eunchosen
                case 1
                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_unchosen(chosen_trial_idx));
                otherwise
                    error('ready yet');
            end
        end

        % effort non-default option
        if chosenModel_EnonDefault > 0
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'effort non-default option';
            switch chosenModel_EnonDefault
                case 1
                    chosen_modVals(n_chosenMods,:) = E_varOption(chosen_trial_idx); % binary variable => no zscore
                otherwise
                    error('ready yet');
            end
        end

        % confidence
        if chosenModel_confidence > 0
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'confidence';
            switch chosenModel_confidence
                case 1
                    chosen_modVals(n_chosenMods,:) = confidence(chosen_trial_idx); % binary variable => no zscore
                otherwise
                    error('ready yet');
            end
        end

        [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
            ['dispChosen_',RP_chosen_nm], modelChosenOnset, modelChosenDur,...
            n_chosenMods, chosen_modNames, chosen_modVals,...
            orth_vars);
    end % RP
end % model chosen period

%% fixation cross before effort
preEffortCrossModel = GLMprm.model_onset.(task_id).preEffortCross;
if ismember(preEffortCrossModel,{'stick','boxcar'})
    iCond = iCond + 1;
    switch preEffortCrossModel
        case 'stick'
            modelPreEffortCrossdur = 0;
        case 'boxcar'
            modelPreEffortCrossdur = preEffortCrossDur;
    end
    [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
        'preEffort fixation cross', preEffortCrossOnsets, modelPreEffortCrossdur, 0, '', [], orth_vars);
end

%% effort performance
EperfModel = GLMprm.model_onset.(task_id).Eperf;
if ismember(EperfModel,{'stick','boxcar'})

    for iRP_Eperf = 1:length(RPperfCond)
        RP_Eperf_nm = RPperfCond{iRP_Eperf};
        %
        EperfModel_money_chosen = GLMprm.Eperf.(task_id).(RP_Eperf_nm).money_chosen;
        EperfModel_effort_chosen = GLMprm.Eperf.(task_id).(RP_Eperf_nm).E_chosen;
        EperfModel_RT1stAnswer = GLMprm.Eperf.(task_id).(RP_Eperf_nm).RT_1stAnswer;

        % extract trial index for the current loop
        switch RP_Eperf_nm
            case 'RP'
                Eperf_trial_idx = 1:length(RP_var_binary);
            case 'R'
                Eperf_trial_idx = RP_var_binary == 1;
            case 'P'
                Eperf_trial_idx = RP_var_binary == 0;
        end
        
        % onset
        iCond = iCond + 1;
        modelEperfOnset = EperfOnsets(Eperf_trial_idx);
        % duration
        switch EperfModel
            case 'stick'
                modelEperfDur = 0;
            case 'boxcar'
                modelEperfDur = EperfDur(Eperf_trial_idx);
        end

        % modulators
        n_EperfMods = 0;
        Eperf_modNames = cell(1,1);
        Eperf_modVals = [];

        % money chosen
        if EperfModel_money_chosen > 0
            n_EperfMods = n_EperfMods + 1;
            Eperf_modNames{n_EperfMods} = 'money chosen';
            switch EperfModel_money_chosen
                case 1
                    Eperf_modVals(n_EperfMods,:) = raw_or_z(money_amount_chosen(Eperf_trial_idx));
                case 2
                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_amount_chosen(Eperf_trial_idx));
                case 3
                    Eperf_modVals(n_EperfMods,:) = raw_or_z(money_level_unchosen(Eperf_trial_idx));
                case 4
                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_level_unchosen(Eperf_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % effort chosen
        if EperfModel_effort_chosen > 0
            n_EperfMods = n_EperfMods + 1;
            Eperf_modNames{n_EperfMods} = 'effort chosen';
            switch EperfModel_effort_chosen
                case 1
                    Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen(Eperf_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end
        
        % RT 1st answer
        if EperfModel_RT1stAnswer > 0
            n_EperfMods = n_EperfMods + 1;
            Eperf_modNames{n_EperfMods} = 'RT 1st answer';
            switch EperfModel_RT1stAnswer
                case 1
                    error('not ready yet');
                    %             Eperf_modVals(n_EperfMods,:) = ;
                otherwise
                    error('not ready yet');
            end
        end

        [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
            ['Eperf_',RP_Eperf_nm], modelEperfOnset, modelEperfDur,...
            n_EperfMods, Eperf_modNames, Eperf_modVals,...
            orth_vars);
    end % RP
end % model chosen period

%% feedback
fbkModel = GLMprm.model_onset.(task_id).fbk;
if ismember(fbkModel,{'stick','boxcar'})

    for iRP_fbk = 1:length(RPfbkCond)
        RP_fbk_nm = RPfbkCond{iRP_fbk};
        %
        fbkModel_moneyObtained = GLMprm.fbk.(task_id).(RP_fbk_nm).money_obtained;
        fbkModel_winVSloss = GLMprm.fbk.(task_id).(RP_fbk_nm).win_vs_loss;
        fbkModel_Emade = GLMprm.fbk.(task_id).(RP_fbk_nm).E_made;
        fbkModel_confidence = GLMprm.fbk.(task_id).(RP_fbk_nm).confidence;

        % extract trial index for the current loop
        switch RP_fbk_nm
            case 'RP'
                fbk_trial_idx = 1:length(RP_var_binary);
            case 'R'
                fbk_trial_idx = RP_var_binary == 1;
            case 'P'
                fbk_trial_idx = RP_var_binary == 0;
        end
        
        % onset
        iCond = iCond + 1;
        modelFbkOnset = fbkOnsets(fbk_trial_idx);
        % duration
        switch fbkModel
            case 'stick'
                modelFbkDur = 0;
            case 'boxcar'
                modelFbkDur = fbkDur(fbk_trial_idx);
        end

        % modulators
        n_fbkMods = 0;
        fbk_modNames = cell(1,1);
        fbk_modVals = [];

        % win vs loss
        if fbkModel_winVSloss > 0
            if strcmp(RP_fbk_nm,'RP')
                n_fbkMods = n_fbkMods + 1;
                fbk_modNames{n_fbkMods} = 'win vs loss';
                switch fbkModel_winVSloss % binary variable => no zscore
                    case 1
                        fbk_modVals(n_fbkMods,:) = win_vs_loss_fbk;
                    otherwise
                        error('not ready yet');
                end
            else
                error('cannot split R and P trials and add a variable representing R/P trials') ;
            end
        end

        % money chosen
        if fbkModel_moneyObtained > 0
            n_fbkMods = n_fbkMods + 1;
            fbk_modNames{n_fbkMods} = 'money obtained';
            switch fbkModel_moneyObtained
                case 1
                    fbk_modVals(n_fbkMods,:) = raw_or_z(money_amount_obtained(fbk_trial_idx));
                otherwise
                    error('not ready yet');
            end
        end

        % effort performed
        if fbkModel_Emade > 0
            n_fbkMods = n_fbkMods + 1;
            fbk_modNames{n_fbkMods} = 'effort made';
            switch fbkModel_Emade
                case 1
                    error('not ready yet');
                    %             fbk_modVals(n_fbkMods,:) = ;
                otherwise
                    error('not ready yet');
            end
        end

        % confidence
        if fbkModel_confidence > 0
            n_fbkMods = n_fbkMods + 1;
            fbk_modNames{n_fbkMods} = 'confidence';
            switch fbkModel_confidence
                case 1
                    error('not ready yet');
                    %             fbk_modVals(n_fbkMods,:) = ;
                otherwise
                    error('not ready yet');
            end
        end

        [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
            ['fbk_',RP_fbk_nm], modelFbkOnset, modelFbkDur,...
            n_fbkMods, fbk_modNames, fbk_modVals,...
            orth_vars);
    end % RP
end % model chosen period

end % function