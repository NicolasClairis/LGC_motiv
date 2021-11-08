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
EperfOnset = NaN(1,n_trials);
for iTrial = 1:n_trials
    if choiceMissedTrials(iTrial) == 0 % if no choice was made during this trial, no effort was made either
        switch task_nm
            case 'physical'
                EperfOnset(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.effort_phase - T0;
            case 'mental'
                EperfOnset(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.nb_1 - T0;
        end
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
        dispChosenDur = EperfOnset - dispChosenOnsets;
    else
        dispChosenDur = preEffortCrossOnsets - dispChosenOnsets;
        preEffortCrossDur = EperfOnset - preEffortCrossOnsets;
    end
    EperfDur = fbkOnsets - EperfOnset;
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
% loading money choice
money_amount_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.left.*RP_var;
money_amount_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.right.*RP_var;
monetary_amount_sum = money_amount_left + money_amount_right;
% loading effort choice
E_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.left.*RP_var;
E_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.left.*RP_var;
E_sum = E_left + E_right;
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
money_chosen = behavioralDataStruct.(task_behavioral_id).R_chosen;
abs_money_chosen = abs(money_chosen);
money_unchosen = money_amount_left.*(choice_LR == 1) + money_amount_right.*(choice_LR == -1);
moneyChosen_min_moneyUnchosen = money_chosen - money_unchosen;
absMoneyChosen_min_moneyUnchosen = abs(money_chosen - money_unchosen);
% loading feedback
money_obtained = behavioralDataStruct.(task_behavioral_id).gain; % could be different from reward chosen (in case of failure) but mostly similar
win_vs_loss_fbk = money_obtained > 0;

%% remove trials where no choice was performed
choiceMiss_onsets = [];
if sum(choiceMissedTrials) > 0
    % onsets
    dispChoiceOptionOnsets(choiceMissedTrials) = [];
    choiceOnsets(choiceMissedTrials) = [];
    dispChosenOnsets(choiceMissedTrials) = [];
    EperfOnset(choiceMissedTrials) = [];
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
    monetary_amount_sum(choiceMissedTrials) = [];
    E_left(choiceMissedTrials) = [];
    E_right(choiceMissedTrials) = [];
    E_sum(choiceMissedTrials) = [];
    choice_LRandConf(choiceMissedTrials) = [];
    choice_LR(choiceMissedTrials) = [];
    confidence(choiceMissedTrials) = [];
    E_chosen(choiceMissedTrials) = [];
    E_unchosen(choiceMissedTrials) = [];
    E_chosen_min_E_unchosen(choiceMissedTrials) = [];
    money_chosen(choiceMissedTrials) = [];
    abs_money_chosen(choiceMissedTrials) = [];
    money_unchosen(choiceMissedTrials) = [];
    moneyChosen_min_moneyUnchosen(choiceMissedTrials) = [];
    absMoneyChosen_min_moneyUnchosen(choiceMissedTrials) = [];
    money_obtained(choiceMissedTrials) = [];
    win_vs_loss_fbk(choiceMissedTrials) = [];
    choice_RT(choiceMissedTrials) = [];
end

% remove trials where effort was not successfull
warning('trials where effort had to be repeated because of subject failure, should also be removed in the future');

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
chosen_RPpool = GLMprm.choice.(task_id).RPpool;
switch chosen_RPpool
    case 0
        RPchosenCond = {'R','P'};
    case 1
        RPchosenCond = {'RP'};
end
Eperf_RPpool = GLMprm.choice.(task_id).RPpool;
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
        choiceModel_moneySum                        = GLMprm.choice.(task_id).(RP_choice_nm).money_sum;
        choiceModel_E_sum                           = GLMprm.choice.(task_id).(RP_choice_nm).E_sum;
        choiceModel_moneyChosen                     = GLMprm.choice.(task_id).(RP_choice_nm).money_chosen;
        choiceModel_moneyUnchosen                   = GLMprm.choice.(task_id).(RP_choice_nm).money_unchosen;
        choiceModel_money_chosen_min_money_unchosen = GLMprm.choice.(task_id).(RP_choice_nm).money_ch_min_unch;
        choiceModel_E_chosen                        = GLMprm.choice.(task_id).(RP_choice_nm).E_chosen;
        choiceModel_E_unchosen                      = GLMprm.choice.(task_id).(RP_choice_nm).E_unchosen;
        choiceModel_E_chosen_min_E_unchosen         = GLMprm.choice.(task_id).(RP_choice_nm).E_ch_min_unch;
        choiceModel_RT                              = GLMprm.choice.(task_id).(RP_choice_nm).RT;
        choiceModel_conf                            = GLMprm.choice.(task_id).(RP_choice_nm).confidence;
        
        if ~strcmp(RP_choice_nm,'RP')
            error('just split onsets and regressors depending on reward/punishment trial');
        end
        iCond = iCond + 1;
        switch choiceModel
            case 'stick'
                modelChoiceDur = 0;
            case 'boxcar'
                modelChoiceDur = dispChoiceOptionsDur;
        end
        
        % modulators
        n_choiceMods = 0;
        choice_modNames = cell(1,1);
        choice_modVals = [];
        
        % reward vs punishments
        if choiceModel_RP == 1
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'R vs P';
            choice_modVals(n_choiceMods,:) = RP_var_binary; % binary variable => no zscore
        end
        
        % money left
        % money right
        
        % money chosen
        switch choiceModel_moneyChosen
            case 1
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'money chosen';
                choice_modVals(n_choiceMods,:) = raw_or_z(money_chosen);
            case 2
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'money chosen';
                choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_chosen);
        end
        
        % money unchosen
        switch choiceModel_moneyUnchosen
            case 1
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'money unchosen';
                choice_modVals(n_choiceMods,:) = raw_or_z(money_unchosen);
        end
        
        % money chosen - money unchosen
        switch choiceModel_money_chosen_min_money_unchosen
            case 1
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'money chosen - unchosen';
                choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen);
        end
        
        % sum of money
        if choiceModel_moneySum == 1
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'money sum';
            choice_modVals(n_choiceMods,:) = raw_or_z(monetary_amount_sum);
        end
        
        % effort left
        
        % effort right
        
        % effort chosen
        if choiceModel_E_chosen == 1
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'effort chosen';
            choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen);
        end
        
        % effort unchosen
        switch choiceModel_E_unchosen
            case 1
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'effort unchosen';
                choice_modVals(n_choiceMods,:) = raw_or_z(E_unchosen);
        end
        
        % effort chosen - unchosen
        switch choiceModel_E_chosen_min_E_unchosen
            case 1
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'effort chosen - unchosen';
                choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_min_E_unchosen);
        end
        
        % sum of efforts
        if choiceModel_E_sum == 1
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'effort sum';
            choice_modVals(n_choiceMods,:) = raw_or_z(E_sum);
        end
        
        % choice confidence
        switch choiceModel_conf
            case 1
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'confidence';
                choice_modVals(n_choiceMods,:) = confidence; % binary variable => no zscore
        end
        
        % RT
        if choiceModel_RT == 1
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'choice RT';
            choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT);
        end
        
        [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
            ['choice_',RP_choice_nm], dispChoiceOptionOnsets, modelChoiceDur,...
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
        chosenModel_Echosen = GLMprm.chosen.(task_id).(RP_chosen_nm).E_chosen;
        chosenModel_moneyUnchosen = GLMprm.chosen.(task_id).(RP_chosen_nm).money_unchosen;
        chosenModel_Eunchosen = GLMprm.chosen.(task_id).(RP_chosen_nm).E_unchosen;
        chosenModel_confidence = GLMprm.chosen.(task_id).(RP_chosen_nm).confidence;
        
        if ~strcmp(RP_chosen_nm,'RP')
            error('just split onsets and regressors depending on reward/punishment trial');
        end
        iCond = iCond + 1;
        switch chosenModel
            case 'stick'
                modelChosenDur = 0;
            case 'boxcar'
                modelChosenDur = dispChosenDur;
        end
        
        % modulators
        n_chosenMods = 0;
        chosen_modNames = cell(1,1);
        chosen_modVals = [];
        % money chosen
        switch chosenModel_moneyChosen
            case 1
                n_chosenMods = n_chosenMods + 1;
                chosen_modNames{n_chosenMods} = 'money chosen';
                chosen_modVals(n_chosenMods,:) = raw_or_z(money_chosen);
            case 2
                n_chosenMods = n_chosenMods + 1;
                chosen_modNames{n_chosenMods} = 'money chosen';
                chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_chosen);
        end
        % effort chosen
        if chosenModel_Echosen == 1
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'effort chosen';
            chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen);
        end
        % money unchosen
        if chosenModel_moneyUnchosen == 1
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'money unchosen';
            chosen_modVals(n_chosenMods,:) = raw_or_z(money_unchosen);
        end
        % effort unchosen
        if chosenModel_Eunchosen == 1
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'effort unchosen';
            chosen_modVals(n_chosenMods,:) = raw_or_z(E_unchosen);
        end
        % confidence
        if chosenModel_confidence == 1
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'confidence';
            chosen_modVals(n_chosenMods,:) = confidence; % binary variable => no zscore
        end
            
        [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
            ['dispChosen_',RP_chosen_nm], dispChosenOnsets, modelChosenDur,...
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
        
        if ~strcmp(RP_Eperf_nm,'RP')
            error('just split onsets and regressors depending on reward/punishment trial');
        end
        iCond = iCond + 1;
        switch EperfModel
            case 'stick'
                modelEperfDur = 0;
            case 'boxcar'
                modelEperfDur = EperfDur;
        end
        
        % modulators
        n_EperfMods = 0;
        Eperf_modNames = cell(1,1);
        Eperf_modVals = [];
        % money chosen
        switch EperfModel_money_chosen
            case 1
                n_EperfMods = n_EperfMods + 1;
                Eperf_modNames{n_EperfMods} = 'money chosen';
                Eperf_modVals(n_EperfMods,:) = raw_or_z(money_chosen);
            case 2
                n_EperfMods = n_EperfMods + 1;
                Eperf_modNames{n_EperfMods} = 'money chosen';
                Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_chosen);
        end
        % effort chosen
        if EperfModel_effort_chosen == 1
            n_EperfMods = n_EperfMods + 1;
            Eperf_modNames{n_EperfMods} = 'effort chosen';
            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen);
        end
        %
        % RT 1st answer
        if EperfModel_RT1stAnswer == 1
            n_EperfMods = n_EperfMods + 1;
            Eperf_modNames{n_EperfMods} = 'RT 1st answer';
            error('not ready yet');
%             Eperf_modVals(n_EperfMods,:) = ;
        end
        
        [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
            ['Eperf_',RP_Eperf_nm], EperfOnset, modelEperfDur,...
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
        
        if ~strcmp(RP_fbk_nm,'RP')
            error('just split onsets and regressors depending on reward/punishment trial');
        end
        iCond = iCond + 1;
        switch fbkModel
            case 'stick'
                modelFbkDur = 0;
            case 'boxcar'
                modelFbkDur = fbkDur;
        end
        
        % modulators
        n_fbkMods = 0;
        fbk_modNames = cell(1,1);
        fbk_modVals = [];
        % win vs loss
        if fbkModel_winVSloss == 1
            n_fbkMods = n_fbkMods + 1;
            fbk_modNames{n_fbkMods} = 'win vs loss';
            fbk_modVals(n_fbkMods,:) = win_vs_loss_fbk; % binary variable => no zscore
        end
        % money chosen
        if fbkModel_moneyObtained == 1
            n_fbkMods = n_fbkMods + 1;
            fbk_modNames{n_fbkMods} = 'money obtained';
            fbk_modVals(n_fbkMods,:) = raw_or_z(money_obtained);
        end
        % effort performed
        if fbkModel_Emade == 1
            n_fbkMods = n_fbkMods + 1;
            fbk_modNames{n_fbkMods} = 'effort made';
            error('not ready yet');
%             fbk_modVals(n_fbkMods,:) = ;
        end
        
        [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
            ['fbk_',RP_fbk_nm], fbkOnsets, modelFbkDur,...
            n_fbkMods, fbk_modNames, fbk_modVals,...
            orth_vars);
    end % RP
end % model chosen period

end % function