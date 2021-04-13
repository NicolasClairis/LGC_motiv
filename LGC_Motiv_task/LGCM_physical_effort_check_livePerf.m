function[timeNow, force_levels, onsets, currentAngle, has_F_threshold_ever_been_reached, was_F_above_threshold] = LGCM_physical_effort_check_livePerf(scr, totalAngleDistance,...
    force_levels,...
    F_threshold, F_tolerance, t_effort_to_keep,...
    onsets)
% [timeNow, force_levels, onsets, currentAngle, has_F_threshold_ever_been_reached, was_F_above_threshold] = LGCM_physical_effort_check_livePerf(scr, totalAngleDistance,...
%     force_levels,...
%     F_threshold, F_tolerance, t_effort_to_keep,...
%     onsets)
% LGCM_physical_effort_check_livePerf checks live force being done and
% displays it accordingly
%
% INPUTS
% scr: structure with screen parameters
%
% totalAngleDistance: total angle distance
%
% force_levels: vector with all force levels of the current trial
%
% F_threshold: force threshold to reach
%
% F_tolerance: tolerance around the threshold
%
% t_effort_to_keep: time you need to keep the perf above the threshold to
% consider that the trial is a success
%
% onsets: structure with all the onsets of the current trial
%
% OUTPUTS
% timeNow: time now in seconds (using GetSecs)
%
% force_levels: vector with all force levels during the current trial
% (would be nice to also record corresponding timepoints)
%
% onsets: updated onsets variable
% 
% currentAngle: current angle
%
% has_F_threshold_ever_been_reached
% (0) 
% (1) 
%
% was_F_above_threshold
% (0)
% (1)

%% initialize indicators whether force threshold has been reached or not
has_F_threshold_ever_been_reached = 0;
was_F_above_threshold = 0;

%% check timing
timeNow = GetSecs;

%% read current level of force
F_now;
% store force levels in the output
force_levels = [force_levels; F_now];

% display effort scale on the left for the live display of the effort level
% add force threshold to maintain and current level of force being made
disp_realtime_force(scr, F_threshold, F_tolerance, F_now);

%% update the center display according to if force above or below the threshold
if F_now >= F_threshold - F_tolerance % while force > threshold, update the timer
    
    
    %% if force was not above the threshold
    % => store when the force started above
    if was_F_above_threshold == 0
        %% if force threshold had not been reached yet
        % => update this
        % => store the timing
        if has_F_threshold_ever_been_reached == 0
            has_F_threshold_ever_been_reached = 1;
            onsets.effort.start1 = timeNow;
        else
            n_starts = fieldnames(onset.Effort);
            onsets.effort.(['start',num2str(n_starts + 1)]) = timeNow;
        end
        onsets.effort_starts = [timeNow; onsets.effort_starts];
    end % identify if transition from below to above threshold or not
    
    %% update current angle based on time passed and time remaining until achieving performance
    timeEffortStart_tmp = max(onsets.effort_starts);
    percentageTimeForceAlreadyMaintained = (timeNow - timeEffortStart_tmp)/t_effort_to_keep;
    currentAngle = totalAngleDistance*percentageTimeForceAlreadyMaintained;
    
    %% display principal informations
    
    % display real-time force level
    disp_realtime_force(scr, F_threshold, F_tolerance, F_now);
    
    % display performance achieved level on top of the reward as an arc
    Screen('FillArc', window,...
        stim.difficulty.currLevelColor,...
        stim.difficulty.middle_center,...
        currentAngle,...
        endAngle - currentAngle);
    
    [~,timeNow] = Screen('Flip',window);
    
    %% update indicator if F is or not above threshold
    was_F_above_threshold = 1;
    
else
    %% update angle restart at zero
    currentAngle = startAngle;
    
    %% update indicator if F is or not above threshold
    was_F_above_threshold = 0;
end
end % function