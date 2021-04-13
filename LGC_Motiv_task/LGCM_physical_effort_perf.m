function[physicalE_perf, trial_success, onsets] = LGCM_physical_effort_perf(scr, stim,...
    E_chosen,...
    E_time_levels,...
    F_threshold, F_tolerance,...
    time_limit, t_max_effort)
% [physicalE_perf, trial_success, onsets] = LGCM_physical_effort_perf(scr, stim,...
%     E_chosen,...
%     E_time_levels,...
%     F_threshold, F_tolerance,...
%     time_limit, t_max_effort)
% LGCM_physical_effort_perf will check if physical effort was performed as
% requested in the available time.
%
% INPUTS
% scr: structure with screen data
%
% stim: structure with data around the stimuli to display
%
% E_chosen: effort level of the chosen option
%
% E_time_levels: structure which stores the timing to reach for each effort
% level (in seconds)
%
% F_threshold: threshold of force (% of MVC) to maintain
%
% F_tolerance: percentage of MVC tolerance we accept around the threshold
%
% time_limit:
% (false): no time limit = keep doing until performance reached
% (true): time limit = fixed duration for this period
%
% t_max_effort: maximum time until when the required effort can be performed
%
% OUTPUTS
% physicalE_perf: structure with summary of physical effort performance
%   .nTrials: number of trials it took to reach a correct
% performance
%   .rt: reaction time (rt) for each question
%   .taskType: task type (0: odd/even; 1: lower/higher than 5)
%   .totalTime_success: time passed between trial onset and last correct
%   answer (if the trial was a success), NaN if trial was a failure
%
% trial_success:
% false: failure
% true: effort was correctly performed in the requested time
%
% onsets: structure with onset of each number on the screen
%   .nb_i: onset of (i) question during the current trial
%
%
% See also LGCM_choice_task_main.m

%% initialize the variables of interest
trial_success = 0;
% angle values
startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen)]);
currentAngle = startAngle;
endAngle = stim.difficulty.arcEndAngle;
totalAngleDistance = endAngle - startAngle;

% effort time to keep
onsets.effort_starts = [];
t_effort_to_keep = E_time_levels.(['level_',num2str(E_chosen)]);

%% display all relevant variables on the screen for the effort period

% display effort scale on the left for the live display of the effort level
% add force threshold to maintain
disp_realtime_force(scr, F_threshold, F_tolerance, 0);

% display difficulty level on top of the reward as an arc
Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
    stim.difficulty.middle_center,...
    startAngle,...
    endAngle - startAngle);

[~,onsetEffortPhase] = Screen('Flip',window);
onsets.effort_phase = onsetEffortPhase;
timeNow = onsetEffortPhase;

%% start timer: if they don't achieve the effort in time they will have to repeat without the reward
force_levels = [];
if time_limit == true
    while timeNow <= onsetEffortPhase + t_max_effort
        [timeNow, force_levels, onsets, currentAngle, was_F_above_threshold] = LGCM_physical_effort_check_livePerf(scr, totalAngleDistance,...
            force_levels,...
            F_threshold, F_tolerance, t_effort_to_keep,...
            was_F_above_threshold, has_F_threshold_ever_been_reached, onsets);
    end % time loop
    
    %% check whether performance was or not achieved in the time course available
    if currentAngle >= endAngle
        trial_success = 1;
    end
    
else % no trial time limit (for training for example)
    while trial_success == 0
        
        % keep performing the effort while the trial is not a success
        [timeNow, force_levels, onsets, currentAngle, was_F_above_threshold] = LGCM_physical_effort_check_livePerf(scr, totalAngleDistance,...
            force_levels,...
            F_threshold, F_tolerance, t_effort_to_keep,...
            was_F_above_threshold, has_F_threshold_ever_been_reached, onsets);
        
        % check whether performance was or not achieved
        if currentAngle >= endAngle
            trial_success = 1;
            onsets.effort_success = timeNow;
        end
    end
end % time limit

%% record vars of interest
physicalE_perf.trial_success = trial_success;
physicalE_perf.onsets = onsets;
% record all the force levels during the performance
physicalE_perf.force_levels = force_levels;

end % function