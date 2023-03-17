function[hI_level] = extract_hI_level(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [hI_level] = extract_hI_level(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% hI_level: 1*nTrials vector with information about the incentive level 
% associated to the high effort option for the current study, subject and
% run. Note that the value will be given in only 0/1/2/3 levels pooling
% both reward and punishment together.

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
%% re-organize to account for punishments order (3/2/1): 
% highest incentive is the lowest punishment level
hI_level = hI_level_v0;
hI_level((hI_level_v0 == 3).*(RP_var_tmp == -1) == 1) = 1;
hI_level((hI_level_v0 == 1).*(RP_var_tmp == -1) == 1) = 3;
        
end % function