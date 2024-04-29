function[RP_trial] = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [RP_trial] = extract_money_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% RP_trial: simple vector indicating whether trial was Reward (+1) or
% Punishment (-1)

%% load data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
    '_task.mat']);
choiceOptions = behaviorStruct.choice_opt;
%% extract R/P (R=+1; P=-1)
RP_trial = strcmp(choiceOptions.R_or_P,'R') - strcmp(choiceOptions.R_or_P,'P');

end % function