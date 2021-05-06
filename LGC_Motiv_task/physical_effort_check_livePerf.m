function[timeNow, force_levels, onsets, currentAngle,...
    has_F_threshold_ever_been_reached] = physical_effort_check_livePerf(scr, dq,...
    MVC,...
    totalAngleDistance,...
    currentAngle, force_levels, has_F_threshold_ever_been_reached,...
    F_threshold, F_tolerance, t_effort_to_keep,...
    onsets)
% [timeNow, force_levels, onsets, currentAngle,...
%     has_F_threshold_ever_been_reached, was_F_above_threshold] = physical_effort_check_livePerf(scr, dq,...
%     MVC,...
%     totalAngleDistance,...
%     currentAngle, force_levels,...
%     F_threshold, F_tolerance, t_effort_to_keep,...
%     onsets)
% physical_effort_check_livePerf checks live force being done and
% displays it accordingly
%
% INPUTS
% scr: structure with screen parameters
%
% dq: information about the handgrip
%
% MVC: maximum voluntary contraction force done during the first
% calibration (expressed in Voltage)
%
% totalAngleDistance: total angle distance
%
% currentAngle: current angle of the arc (to keep track of how long the
% performance has or not been achieved)
%
% force_levels: vector with all force levels of the current trial
%
% has_F_threshold_ever_been_reached:
% (0) the threshold has not yet been reached ever in the course of the
% trial
% (1) the threshold was reached at least once in the course of the trial
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
% onsets: updated onsets variable; structure with all the relevant
% information of the trial
% 
% currentAngle: current angle for the arc displaying the performance
%
% has_F_threshold_ever_been_reached
% (0) the threshold has not been reached yet
% (1) the threshold was reached
%

%% read current level of force
% check timing now when we extract the level of force
timeNow = GetSecs;
% check force being applied now
F_now_Voltage_table = read(dq);
F_now_Voltage = F_now_Voltage_table.Variables;
% convert force level from Voltage to a percentage of MVC
F_now = F_now_Voltage/MVC;
% store force levels in the output
force_levels = [force_levels;...
    [F_now, timeNow, F_now_Voltage]]; % store F in % of MVC, time and F in Volts

%% update the center display according to if force above or below the threshold
if F_now >= (F_threshold - F_tolerance) % while force > threshold, update the timer
    % note: all force levels here are expressed in percentage of MVC
    
    %% if force was not above the threshold
    % => store when the force started above
    if currentAngle > startAngle
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
    
else
    %% reset angle
    currentAngle = startAngle;
end % is current level of force above the requested threshold?

%% display on screen accordingly

% display real-time force level
disp_realtime_force(scr, F_threshold, F_tolerance, F_now, 'task');

% display performance achieved level on top of the reward as an arc
Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
    stim.difficulty.middle_center,...
    currentAngle,...
    endAngle - currentAngle);

[~,timeNow] = Screen('Flip',window);

end % function