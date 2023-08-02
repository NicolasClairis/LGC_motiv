function[P_level_chosen] = extract_P_level_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [P_level_chosen] = extract_P_level_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% P_level_chosen: 1*nTrials vector with information about the punishment level 
% chosen for the current study, subject and run. Note that the value will 
% be given in only 1/2/3/4 levels for reward and 0 for when the trial is a 
% reward trial.

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
%% extract high reward level
hI_level_v0 = (choiceOptions.R.left).*(choice_LR == 1) +...
    (choiceOptions.R.right).*(choice_LR == -1);
% adding (+1) to distinguish choice of default option (now equal +1)
% from reward trials (P_level_chosen = 0)
hI_level_v1 = hI_level_v0 + 1;
%% valence
RP_var_tmp = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
%% store and remove all punishment trials
P_level_chosen = hI_level_v1;
P_level_chosen(RP_var_tmp == 1) = 0;
P_level_chosen_forSwap = P_level_chosen;

% swap low and high punishment numbers to make it more similar to rewards
% (P2 will now be the least motivating amount (ie worst punishment, closer
% to the default option P1)
% P4 will now be the more motivating amount (ie "best" punishment))
P_level_chosen(P_level_chosen_forSwap == 4) = 2;
P_level_chosen(P_level_chosen_forSwap == 2) = 4;

end % function