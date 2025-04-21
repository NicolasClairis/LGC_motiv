function[deltaI_money, low_I_money, high_I_money] = extract_monetary_incentive(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [deltaI_money, lowI, highI] = extract_monetary_incentive(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% extract_monetary_incentive will extract the difference between the high and the
% low effort options in monetary amount.
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
% deltaI_money: 1*nTrials vector with information about the difference
% between the high effort option monetary amount and the baseline across
% all trials (including reward and punishment together)
%
% low_I_money: 1*nTrials vector with amount associated to the low effort option
%
% high_I_money: 1*nTrials vector with amount associated to the high effort option

%% load data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
    '_task.mat']);
choiceOptions = behaviorStruct.choice_opt;

%% default side
defaultSide = choiceOptions.default_LR;
%% extract high effort level
% extract side of high effort (left or right) to know which is the low and
% which is the high effort option
highE_left = defaultSide == 1;
highE_right = defaultSide == -1;
% extract valence (reward or punishment trial?) to have signed values
RP_trial = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
% extract low and high effort option monetary incentive
low_I_money = ((choiceOptions.monetary_amount.left).*highE_right +...
    (choiceOptions.monetary_amount.right).*highE_left).*RP_trial;
high_I_money = ((choiceOptions.monetary_amount.left).*highE_left +...
    (choiceOptions.monetary_amount.right).*highE_right).*RP_trial;
% extract difference between high and low monetary incentives
deltaI_money = high_I_money - low_I_money;
        
end % function