function[mental_taskType_trialStart] = LGCM_mental_task_start(n_trials)
%[mental_taskType_trialStart] = LGCM_mental_task_start(n_trials)
%
% INPUTS
% n_trials: number of trials
%
% OUTPUTS
% mental_taskType_trialStart: vector which tells you what is the type of
% the first trial (odd/even task or higher/lower than 5 task)
%

%% split trials in half
n_half_trials = floor(n_trials/2);
% be careful here this means that if you have an impair number of trials,
% one of the two tasks will be systematically starting for one more trial
% than the other

%% define all trials type
task_type_list = [zeros(1, n_half_trials), ones(1, n_half_trials)];

%% randomize the order
trial_idx = randperm(n_trials);
mental_taskType_trialStart = task_type_list(trial_idx);

end % function