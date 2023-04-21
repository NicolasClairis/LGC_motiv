function[hR_level] = extract_hR_level(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [hR_level] = extract_hR_level(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% hR_level: 1*nTrials vector with information about the reward level 
% associated to the high effort option for the current study, subject and
% run. Note that the value will be given in only 1/2/3 levels for reward 
% and 0 for when the trial is a punishment trial.

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
%% store and remove all punishment trials
hR_level = hI_level_v0;
hR_level(RP_var_tmp == -1) = 0;
        
end % function