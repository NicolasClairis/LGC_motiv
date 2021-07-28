function[matlabbatch] = First_level_loadRegressors(matlabbatch, GLMprm, sub_idx, iRun,...
    subj_behavior_folder, currRunBehaviorFileName, task_nm)
% [matlabbatch] = First_level_loadRegressors(matlabbatch, GLMprm, sub_idx, iRun,...
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
behavioralDataStruct = load([subj_behavior_folder, currRunBehaviorFileName]);
% extract all onsets of interest and corresponding durations
T0 = behavioralDataStruct.onsets.T0;
initialCrossOnsets      = behavioralDataStruct.(task_behavioral_id).onsets.cross - T0; % only cross appearing BEFORE the start of each trial
crossOnsets             = [initialCrossOnsets, behavioralDataStruct.onsets.finalCross - T0]; % add final cross (at the end of the experiment)
dispChoiceOptionOnsets  = behavioralDataStruct.(task_behavioral_id).onsets.dispChoiceOptions - T0;
choiceOnsets            = behavioralDataStruct.(task_behavioral_id).onsets.choice - T0;
dispChosenOnsets        = behavioralDataStruct.(task_behavioral_id).onsets.dispChoice - T0;
n_trials = length(dispChosenOnsets);

EperfOnset = NaN(1,n_trials);
for iTrial = 1:n_trials
    switch task_nm
        case 'physical'
            EperfOnset(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.effort_phase - T0;
        case 'mental'
            EperfOnset(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.nb_1 - T0;
    end
end
fbkOnsets = behavioralDataStruct.(task_behavioral_id).onsets.fbk - T0;
% durations of each event
crossDur = dispChoiceOptionOnsets - initialCrossOnsets;
dispChoiceDur = dispChosenOnsets - dispChoiceOptionOnsets;
choice_RT = choiceOnsets - dispChoiceOptionOnsets;
dispChosenDur = EperfOnset - dispChosenOnsets;
EperfDur = fbkOnsets - EperfOnset;
fbkDur = crossOnsets(2:end) - fbkOnsets;
% extract regressors of interest
RP_var_binary = strcmp( behavioralDataStruct.(task_behavioral_id).choiceOptions.R_or_P, 'R');
RP_var = RP_var_binary;
RP_var(RP_var_binary == 0) = -1;
money_amount_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.left.*RP_var;
money_amount_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.right.*RP_var;
monetary_amount_sum = money_amount_left + money_amount_right;
E_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.left.*RP_var;
E_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.left.*RP_var;
E_sum = E_left + E_right;
money_chosen = behavioralDataStruct.(task_behavioral_id).R_chosen;
% money_unchosen = ;
E_chosen = behavioralDataStruct.(task_behavioral_id).E_chosen;
% E_unchosen = ;
money_obtained = behavioralDataStruct.(task_behavioral_id).gain;
win_vs_loss_fbk = money_obtained > 0;

%% load the batch according to GLMprm variables
orth_vars = GLMprm.gal.orth_vars;
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

%% fixation cross
crossModel = GLMprm.model_onset.(task_id).cross;
if ismember(crossModel,{'stick','boxcar'})
    iCond = iCond + 1;
    switch crossModel
        case 'stick'
            modelCrossdur = 0;
        case 'boxcar'
            modelCrossdur = crossDur;
    end
    [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
        'fixation cross', crossOnsets, modelCrossdur, 0, '', [], orth_vars);
end


%% choice period
choiceModel = GLMprm.model_onset.(task_id).choice;
if ismember(choiceModel,{'stick','boxcar'})
    
    for iRP_choice = 1:length(RPchoiceCond)
        RP_choice_nm = RPchoiceCond{iRP_choice};
        choiceModel_RP = GLMprm.choice.(task_id).(RP_choice_nm).R_vs_P;
        choiceModel_moneySum = GLMprm.choice.(task_id).(RP_choice_nm).money_sum;
        choiceModel_E_sum = GLMprm.choice.(task_id).(RP_choice_nm).E_sum;
        choiceModel_RT = GLMprm.choice.(task_id).(RP_choice_nm).RT;
        
        if ~strcmp(RP_choice_nm,'RP')
            error('just split onsets and regressors depending on reward/punishment trial');
        end
        iCond = iCond + 1;
        switch choiceModel
            case 'stick'
                modelChoiceDur = 0;
            case 'boxcar'
                modelChoiceDur = dispChoiceDur;
        end
        
        % modulators
        n_choiceMods = 0;
        choice_modNames = cell(1,1);
        choice_modVals = [];
        % reward vs punishments
        if choiceModel_RP == 1
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'R vs P';
            choice_modVals(n_choiceMods,:) = RP_var_binary;
        end
        % sum of money
        if choiceModel_moneySum == 1
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'money sum';
            choice_modVals(n_choiceMods,:) = monetary_amount_sum;
        end
        % sum of efforts
        if choiceModel_E_sum == 1
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'effort sum';
            choice_modVals(n_choiceMods,:) = E_sum;
        end
        % RT
        if choiceModel_RT == 1
            n_choiceMods = n_choiceMods + 1;
            choice_modNames{n_choiceMods} = 'choice RT';
            choice_modVals(n_choiceMods,:) = choice_RT;
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
        if chosenModel_moneyChosen == 1
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'money chosen';
            chosen_modVals(n_chosenMods,:) = money_chosen;
        end
        % effort chosen
        if chosenModel_Echosen == 1
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'effort chosen';
            chosen_modVals(n_chosenMods,:) = E_chosen;
        end
        % money unchosen
        if chosenModel_moneyUnchosen == 1
            error('not ready yet');
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'money unchosen';
            chosen_modVals(n_chosenMods,:) = money_unchosen;
        end
        % effort unchosen
        if chosenModel_Eunchosen == 1
            error('not ready yet');
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'effort unchosen';
            chosen_modVals(n_chosenMods,:) = E_unchosen;
        end
        % confidence
        if chosenModel_confidence == 1
            error('not ready yet');
            n_chosenMods = n_chosenMods + 1;
            chosen_modNames{n_chosenMods} = 'confidence';
            chosen_modVals(n_chosenMods,:) = confidence;
        end
            
        [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
            ['dispChosen_',RP_chosen_nm], dispChosenOnsets, modelChosenDur,...
            n_chosenMods, chosen_modNames, chosen_modVals,...
            orth_vars);
    end % RP
end % model chosen period

%% effort performance
EperfModel = GLMprm.model_onset.(task_id).Eperf;
if ismember(EperfModel,{'stick','boxcar'})
    
    for iRP_Eperf = 1:length(RPperfCond)
        RP_Eperf_nm = RPperfCond{iRP_Eperf};
        %
        EperfModel_money = GLMprm.Eperf.(task_id).(RP_Eperf_nm).money;
        EperfModel_effort = GLMprm.Eperf.(task_id).(RP_Eperf_nm).effort;
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
        if EperfModel_money == 1
            n_EperfMods = n_EperfMods + 1;
            Eperf_modNames{n_EperfMods} = 'money chosen';
            Eperf_modVals(n_EperfMods,:) = money_chosen;
        end
        % effort chosen
        if EperfModel_effort == 1
            n_EperfMods = n_EperfMods + 1;
            Eperf_modNames{n_EperfMods} = 'effort chosen';
            Eperf_modVals(n_EperfMods,:) = E_chosen;
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
            fbk_modVals(n_fbkMods,:) = win_vs_loss_fbk;
        end
        % money chosen
        if fbkModel_moneyObtained == 1
            n_fbkMods = n_fbkMods + 1;
            fbk_modNames{n_fbkMods} = 'money obtained';
            fbk_modVals(n_fbkMods,:) = money_obtained;
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