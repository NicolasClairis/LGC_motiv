function[trainingTimes, calibTimes, taskTimes] = timings_definition(n_trainingTrials, nTrials)
%[trainingTimes, calibTimes, taskTimes] = timings_definition(n_trainingTrials, nTrials)
% timings_definition defines the time duration for each period of the
% experiment.
% All timings are expressed in seconds.
%
% INPUTS
% n_trainingTrials: number of training trials
%
% nTrials: number of trials in the main task
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
        %         calibTimes.
        calibTimes.t_MVC_rest = 7; % rest after each MVC calibration
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
jittersTraining = linspace(jitterMin, jitterMax, n_trainingTrials);
jitterTrainingRdmPerm = randperm(n_trainingTrials);
t_trainingCross = jittersTraining(jitterTrainingRdmPerm);
trainingTimes.t_cross = t_trainingCross;
% other times
switch effort_type
    case 'physical' % in case you use different numbers for each effort type
        trainingTimes.max_effort = t_max_effort;
    case 'mental'
        trainingTimes.t_min_scalingFactor = 150/100; % multiply calibrated minimal time by this value
end
trainingTimes.choice        = t_choice;
trainingTimes.dispChoice    = t_dispChoice;
trainingTimes.fbk           = t_fbk;

end % function