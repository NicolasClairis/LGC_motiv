function[hP_level] = extract_hP_level(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [hP_level] = extract_hP_level(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% hP_level: 1*nTrials vector with information about the punishment level 
% associated to the high effort option for the current study, subject and
% run. Note that the value will be given in only 1/2/3 levels for punishment 
% and 0 for when the trial is a reward trial.

%% load data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
    '_task.mat']);
choiceOptions = behaviorStruct.choice_opt;

%% default side
defaultSide = choiceOptions.default_LR;
%% extract high effort level
hI_level_v0 = (choiceOptions.R.left).*(defaultSide == 1) +...
    (choiceOptions.R.right).*(defaultSide == -1);
%% valence
RP_var_tmp = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
%% store and remove all reward trials
hP_level = hI_level_v0;
% remove reward trials
hP_level(RP_var_tmp == 1) = 0;
% swap low and high punishment numbers to make it more similar to rewards
% (P1 will now be the least motivating amount (ie worst punishment, closer 
% to the default option)
% P3 will now be the more motivating amount (ie "best" punishment))
hP_level((hI_level_v0 == 3).*(RP_var_tmp == -1) == 1) = 1;
hP_level((hI_level_v0 == 1).*(RP_var_tmp == -1) == 1) = 3;
        
end % function