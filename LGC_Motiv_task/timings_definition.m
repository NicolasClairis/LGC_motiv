function[trainingTimes, calibTimes, taskTimes] = timings_definition(scr, trainingConditions, n_R_levels, n_E_levels, nTrials, effort_type)
%[trainingTimes, calibTimes, taskTimes] = timings_definition(scr, n_trainingTrials, taskTrainingConditions, n_R_levels, n_E_levels, nTrials, effort_type)
% timings_definition defines the time duration for each period of the
% experiment.
% All timings are expressed in seconds.
%
% INPUTS
% scr: structure with screen parameters
%
% trainingConditions: 'R' or {'R','P','RP'} (reward, punishment, reward and punishment)
% prepare jitter timings accordingly
%
% n_R_levels: number of reward conditions
%
% n_E_levels: number of effort conditions
%
% nTrials: number of trials in the main task
%
% effort_type: string indicating the nature of the current task
% 'mental'
% 'physical'
%
% OUTPUTS
% trainingTimes: structure with training timings
%
% calibTimes: structure with calibration timings
%
% taskTimes: structure with main task timings
%
%

%% calibration timings
calibTimes.instructions = 5;
switch effort_type % in case you use different numbers for each effort type
    case 'mental'
        calibTimes.effort_max = 10; % maximal time to perform the task (calibrated time should be shorter)
    case 'physical'
        calibTimes.instructions_bis = 10;
        calibTimes.effort_max = 5;% time to perform the task
        calibTimes.physicalReadWait = 0.075; % Arthur manual definition
        calibTimes.MVC_rest = 7; % rest after each MVC calibration
        calibTimes.ifi = 1/15; % manual definition to match with read frame rate
end
calibTimes.fbk = 2;

%% main task timings
jitterMin = 0.5;
jitterMax = 3.5;
jitters = linspace(jitterMin, jitterMax, nTrials);
jitterRdmPerm = randperm(nTrials);
t_cross = jitters(jitterRdmPerm);

t_finalCross = 5;
t_choice = 7;
t_dispChoice = 3;
switch effort_type
    case 'physical' % in case you use different numbers for each effort type
        t_max_effort = 5; % time to perform the task
        taskTimes.max_effort = t_max_effort;
        
        % store frame rate for physical effort task
        % query the frame duration (inter-frame interval)
%         taskTimes.ifi = Screen('GetFlipInterval', scr.window);
        taskTimes.ifi = 1/15; % manual definition to match with read frame rate
        
        % define pause duration after read to make it work without losing
        % too much in the display
        taskTimes.physicalReadWait = 0.075; % Arthur manual definition
    case 'mental'
        taskTimes.t_min_scalingFactor = 150/100; % multiply calibrated minimal time by this value
end
t_fbk = 1; % feedback display
taskTimes.cross = t_cross;
taskTimes.choice = t_choice;
taskTimes.dispChoice = t_dispChoice;
taskTimes.feedback = t_fbk;
taskTimes.finalCross = t_finalCross;

%% training timings
trainingTimes.trainingInstructions = 5;
trainingTimes.trainingEnd   = 5;
% jitters for fixation cross during training
n_trainingCond = length(trainingConditions);
for iTraining = 1:n_trainingCond
    trainingCond = trainingConditions{iTraining};
    [~, n_trainingTrials, ~] = training_options(trainingCond, n_R_levels, n_E_levels);
    jittersTraining = linspace(jitterMin, jitterMax, n_trainingTrials);
    jitterTrainingRdmPerm = randperm(n_trainingTrials);
    t_trainingCross = jittersTraining(jitterTrainingRdmPerm);
    trainingTimes.cross.(trainingCond) = t_trainingCross;
end
% other times
switch effort_type
    case 'physical' % in case you use different numbers for each effort type
        trainingTimes.max_effort = t_max_effort;
        
        % store frame rate for physical effort task
        % query the frame duration (inter-frame interval)
%         trainingTimes.ifi = Screen('GetFlipInterval', scr.window);
        trainingTimes.ifi = 1/15; % manual definition to match with read frame rate
        
        % define pause duration after read to make it work without losing
        % too much in the display
        trainingTimes.physicalReadWait = 0.075; % Arthur manual definition
    case 'mental'
        trainingTimes.t_min_scalingFactor = 150/100; % multiply calibrated minimal time by this value
end
trainingTimes.choice        = t_choice;
trainingTimes.dispChoice    = t_dispChoice;
trainingTimes.fbk           = t_fbk;

end % function