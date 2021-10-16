function[trainingTimes, calibTimes, learningTimes, taskTimes, mainTimes] = timings_definition(trainingConditions, n_R_levels, n_E_levels, nTrials, effort_type)
%[trainingTimes, calibTimes, learningTimes, taskTimes, mainTimes] = timings_definition(scr, n_trainingTrials, taskTrainingConditions, n_R_levels, n_E_levels, nTrials, effort_type)
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
% learningTimes: structure with learning timings
%
% taskTimes: structure with main task timings
%
% mainTimes: structure with times useful for both tasks

%% manual calibration values
if strcmp(effort_type, 'physical')
    t_ifi = 1/15;
    t_readWait = 0.075;
end

%% calibration timings
calibTimes.instructions = 5;
switch effort_type % in case you use different numbers for each effort type
    case 'mental'
        calibTimes.effort_max = 12; % maximal time to perform the task (calibrated time should be shorter)
    case 'physical'
        calibTimes.instructions_bis = 6;
        calibTimes.effort_max = 5;% time to perform the task
        calibTimes.physicalReadWait = t_readWait; % Arthur manual definition
        calibTimes.MVC_rest = 7; % rest after each MVC calibration
        calibTimes.ifi = t_ifi; % manual definition to match with read frame rate
end
calibTimes.fbk = 2;
calibTimes.fail_and_repeat_fbk = 5;

%% learning timings
switch effort_type % in case you use different numbers for each effort type
    case 'mental'
        learningTimes.learning_rest = 3;
    case 'physical'
        learningTimes.ifi = t_ifi;
        learningTimes.max_effort = [];
        learningTimes.physicalReadWait = t_readWait;
        learningTimes.learning_rest = 3;
end
learningTimes.fail_and_repeat_fbk = 5;

%% main task timings

% initial cross
jitterMin_choiceCross = 0.5;
jitterMax_choiceCross = 3.5;
jitters_choiceCross = linspace(jitterMin_choiceCross, jitterMax_choiceCross, nTrials);
jitterRdmPerm_choiceCross = randperm(nTrials);
t_choiceCross = jitters_choiceCross(jitterRdmPerm_choiceCross);

% cross between choice and effort
jitterMin_effortCross = 0.5;
jitterMax_effortCross = 1.5;
jitters_effortCross = linspace(jitterMin_effortCross, jitterMax_effortCross, nTrials);
jitterRdmPerm_effortCross = randperm(nTrials);
t_effortCross = jitters_effortCross(jitterRdmPerm_effortCross);

t_finalCross = 10;
t_choice = 15;
t_dispChoice = 2;
switch effort_type
    case 'physical' % in case you use different numbers for each effort type
        t_max_effort = 6; % time to perform the task
        taskTimes.max_effort = t_max_effort;
        
        % store frame rate for physical effort task
        % query the frame duration (inter-frame interval)
        %         taskTimes.ifi = Screen('GetFlipInterval', scr.window);
        taskTimes.ifi = t_ifi; % manual definition to match with read frame rate
        
        % define pause duration after read to make it work without losing
        % too much in the display
        taskTimes.physicalReadWait = t_readWait; % Arthur manual definition
    case 'mental'
        taskTimes.t_min_scalingFactor = 150/100; % multiply calibrated minimal time by this value
end
t_fbk = 1; % feedback display
t_fail_and_repeat_fbk = 3; % feedback after a failure => repeat the effort after that
taskTimes.preChoiceCross.mainTask = t_choiceCross;
taskTimes.preEffortCross.mainTask = t_effortCross;
taskTimes.choice = t_choice;
taskTimes.dispChoice = t_dispChoice;
taskTimes.feedback = t_fbk;
taskTimes.finalCross = t_finalCross;
taskTimes.fail_and_repeat_fbk = t_fail_and_repeat_fbk;

%% training timings
trainingTimes.instructions = 5;
trainingTimes.trainingEnd   = 5;
% jitters for fixation cross during training
n_trainingCond = length(trainingConditions);
for iTraining = 1:n_trainingCond
    trainingCond = trainingConditions{iTraining};
    [~, n_trainingTrials] = training_options(trainingCond, n_R_levels, n_E_levels);
    % initial cross before choice
    jittersTrainingChoice = linspace(jitterMin_choiceCross, jitterMax_choiceCross, n_trainingTrials);
    jitterTrainingChoiceRdmPerm = randperm(n_trainingTrials);
    t_trainingChoiceCross = jittersTrainingChoice(jitterTrainingChoiceRdmPerm);
    % cross before effort
    jittersTrainingEffort = linspace(jitterMin_effortCross, jitterMax_effortCross, n_trainingTrials);
    jitterTrainingEffortRdmPerm = randperm(n_trainingTrials);
    t_trainingEffortCross = jittersTrainingEffort(jitterTrainingEffortRdmPerm);
    % store timings
    trainingTimes.preChoiceCross.(trainingCond) = t_trainingChoiceCross;
    trainingTimes.preEffortCross.(trainingCond) = t_trainingEffortCross;
end
% other times
switch effort_type
    case 'physical' % in case you use different numbers for each effort type
        % max effort time
        trainingTimes.max_effort = t_max_effort;
        
        % store frame rate for physical effort task
        % query the frame duration (inter-frame interval)
        %         trainingTimes.ifi = Screen('GetFlipInterval', scr.window);
        trainingTimes.ifi = t_ifi; % manual definition to match with read frame rate
        
        % define pause duration after read to make it work without losing
        % too much in the display
        trainingTimes.physicalReadWait = t_readWait; % Arthur manual definition
    case 'mental'
        trainingTimes.t_min_scalingFactor = 150/100; % multiply calibrated minimal time by this value
end
trainingTimes.choice                = t_choice;
trainingTimes.dispChoice            = t_dispChoice;
trainingTimes.feedback              = t_fbk;
trainingTimes.fail_and_repeat_fbk   = t_fail_and_repeat_fbk;

%% time feedback end of a block
mainTimes.endfMRI = 20; % to get end of fMRI response and avoid artifacts
mainTimes.endSession = 3; % display of the total amount of money obtained

end % function