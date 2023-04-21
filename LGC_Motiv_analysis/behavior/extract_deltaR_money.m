function[deltaR_money] = extract_deltaR_money(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [deltaR_money] = extract_deltaR_money(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% deltaR_money: 1*nTrials vector with information about the difference
% between the high effort option reward amount and the baseline. For
% punishment trials, will be equal to 0.

%% load data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
    '_task.mat']);
choiceOptions = behaviorStruct.choice_opt;

%% default side
defaultSide = choiceOptions.default_LR;
%% extract high effort level
highE_left = defaultSide == 1;
highE_right = defaultSide == -1;
hI_level_v0 = (choiceOptions.monetary_amount.left).*highE_left +...
    (choiceOptions.monetary_amount.right).*highE_right;
lowI_level_v0 = (choiceOptions.monetary_amount.left).*highE_right +...
    (choiceOptions.monetary_amount.right).*highE_left;
%% valence
RP_var_tmp = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
%% store and remove all punishment trials
deltaR_money = hI_level_v0 - lowI_level_v0;
deltaR_money(RP_var_tmp == -1) = 0;
        
end % function