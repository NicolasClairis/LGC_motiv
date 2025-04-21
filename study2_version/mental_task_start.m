function[mental_taskType_trialStart] = mental_task_start(n_trials)
%[mental_taskType_trialStart] = mental_task_start(n_trials)
%
% mental_task_start will randomize the type of task for the first
% question of each trial
%
% INPUTS
% n_trials: number of trials for the block
%
% OUTPUTS
% mental_taskType_trialStart: (0/1) variable indicating whether the trial
% should start with an odd/even question (0) or lower/higher than 5
% question (1)

%% check n_trials is pair
if mod(n_trials, 2) == 0
    n_trials_eachType = n_trials/2;
else
    error(['total number of trials is not pair => can''t you fixed that? ',...
        'If not, you need to define which question type ',...
        '(odd/even or lower/higher than 5) will have more trials starting with it',...
        ' or try to balance this across participants.']);
end

%% define a big vector with the trial start
mental_taskType_trialStart_bf_orderRandomization = [zeros(1,n_trials_eachType), ones(1,n_trials_eachType)];

%% randomize the order
order_rdm = randperm(n_trials);
mental_taskType_trialStart = mental_taskType_trialStart_bf_orderRandomization(order_rdm);
end % function