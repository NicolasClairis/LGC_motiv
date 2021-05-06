function[task_seq] = mental_effort_task_switches(prev_task_type, n_correctAnswerToProvide, n_switch)
% [task_seq] = mental_effort_task_switches(prev_task_type, n_correctAnswerToProvide, n_switch)
%
% INPUTS
% prev_task_type: initial task type = starting type of the sequence
%
% n_correctAnswerToProvide: total correct answers to provide
%
% n_switch: number of switches to implement
%
% OUTPUTS
% task_seq: sequence indicating task type for each question based on the
% proportion of switches and on the total number of correct answers to
% provide
%

%% initialize variable of interest
task_seq = NaN(1, n_correctAnswerToProvide);
task_seq(1) = prev_task_type;

%% define the sequence
if n_switch > 0
    % randomize the moment when the switch appears
    switch_idx = randperm(n_correctAnswerToProvide - 1, n_switch) + 1;
    
    % for each question, determine if a switch
    for iQuest = 2:n_correctAnswerToProvide
        
        if ~ismember(iQuest, switch_idx) % no switch
            task_seq(iQuest) = task_seq(iQuest - 1);
            
        elseif ismember(iQuest, switch_idx) % switch
            task_seq(iQuest) = 1 - task_seq(iQuest - 1);
        end
    end
    
else % learning session with no switch: all trials are of the same type
    task_seq(:) = prev_task_type;
end

end % function