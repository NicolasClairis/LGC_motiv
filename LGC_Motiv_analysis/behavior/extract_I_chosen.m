function[I_chosen] = extract_I_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [I_chosen] = extract_I_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% I_chosen: 1*nTrials vector with information about the incentive level 
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
%% extract incentive chosen
Ich_v0 = (choiceOptions.R.left).*(choice_LR == -1) +...
    (choiceOptions.R.right).*(choice_LR == 1);
%% valence
RP_var_tmp = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
%% re-organize to account for punishments order (3/2/1):
% highest incentive is the lowest punishment level
I_chosen = Ich_v0;
I_chosen((Ich_v0 == 3).*(RP_var_tmp == -1) == 1) = 1;
I_chosen((Ich_v0 == 1).*(RP_var_tmp == -1) == 1) = 3;
   
end % function