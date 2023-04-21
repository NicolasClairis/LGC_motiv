function[money_chosen, money_level_chosen] = extract_money_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [money_chosen, money_level_chosen] = extract_money_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName)
%
% INPUTS
% subBehaviorFolder: folder where data is stored
%
% sub_nm: string with subject CID name
%
% run_nm: string with run name
%
% task_fullName: task full name 'mental'/'physical'
%
% OUTPUTS
% money_chosen: 1*nTrials vector with information about the money amount 
% chosen for the current study, subject and run.
%
% money_level_chosen: 1*nTrials vector with information about the money level 
% chosen for the current study, subject and run.

%% load data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
    '_task.mat']);
choiceOptions = behaviorStruct.choice_opt;
switch task_fullName
    case 'mental'
        choiceAndPerf = behaviorStruct.mentalE_perf;
    case 'physical'
        choiceAndPerf = behaviorStruct.physicalPerf;
end
%% choice was left or right?
choice_LR = choiceAndPerf.choice;
% remove confidence info from choice:
choice_LR(choice_LR == 2) = 1;
choice_LR(choice_LR == -2) = -1;
%% extract R/P
R_trials = strcmp(choiceOptions.R_or_P,'R');
%% extract money chosen
money_chosen = ((choiceOptions.monetary_amount.left).*(choice_LR == -1) +...
            (choiceOptions.monetary_amount.right).*(choice_LR == 1)).*...
            ((R_trials == 1) - (R_trials == 0));
money_level_chosen = (choiceOptions.R.left).*(choice_LR == -1) +...
            (choiceOptions.R.right).*(choice_LR == 1);
        
end % function