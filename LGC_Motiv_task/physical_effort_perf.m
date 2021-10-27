function[physicalE_perf, trial_success, onsets] = physical_effort_perf(scr, stim, dq,...
    MVC,...
    R_or_P, R_chosen, E_chosen,...
    E_time_levels,...
    F_threshold, F_tolerance,...
    time_limit, timings)
% [physicalE_perf, trial_success, onsets] = physical_effort_perf(scr, stim, dq,...
%     MVC,...
%     E_chosen,...
%     E_time_levels,...
%     F_threshold, F_tolerance,...
%     time_limit, timings)
% physical_effort_perf will check if physical effort was performed as
% requested in the available time.
%
% INPUTS
% scr: structure with screen data
%
% stim: structure with data around the stimuli to display
%
% dq: information about the handgrip
%
% MVC: maximal voluntary contraction force done during first calibration
% (expressed in Voltage)
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
% timings: structure with relevant timings for the task (pause timing, time
% to wait for the effort performance, etc.)
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
% See also choice_task_main.m

%% initialize the variables of interest
% screen and PTB parameters
window = scr.window;
% extract timings
ifi = timings.ifi;
t_max_effort = timings.max_effort;
t_readWait = timings.physicalReadWait;
% by default when the trial starts, the trial cannot be a success already
% => initialize to zero
trial_success = 0;
% initialize indicators whether force threshold has been reached or not
has_F_threshold_ever_been_reached = 0;
% angle values
startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen)]);
currentAngle = startAngle;
endAngle = stim.difficulty.arcEndAngle;
totalAngleDistance = endAngle - startAngle;

% effort time to keep
[onsets.effort.start, onsets.effort.stop] = deal([]);
t_effort_to_keep = E_time_levels.(['level_',num2str(E_chosen)]);

%% start acquiring the data in the background (if you don't use this
% function, everytime you call the read function, it will take a
% long time to process)
start(dq,"continuous");
pause(0.125);
% will need data = read(dq) function only to read the signal

%% display all relevant variables on the screen for the effort period

% display effort scale on the left for the live display of the effort level
% add force threshold to maintain
force_initialDisplay = 0;
disp_realtime_force(scr, F_threshold, F_tolerance, force_initialDisplay, 'task');

% display difficulty level on top of the reward as an arc
Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
    stim.difficulty.middle_center,...
    startAngle,...
    endAngle - startAngle);

[lastFrameTime, onsetEffortPhase] = Screen('Flip',window);
onsets.effort_phase = onsetEffortPhase;
timeNow = onsetEffortPhase;

%% effort
[dispInfos.time,...
    dispInfos.currentAngle,...
    dispInfos.forceLevel] = deal([]);
stateSqueezeON = false;
force_levels = [];

[percSqueeze, percStoppedSqueezing] = deal(0);

% initialize read
[F_now_Voltage, timeCheck, sampleOk_tmp] = F_read(dq, t_readWait);
% if read failed replace by zero for start
if isnan(F_now_Voltage)
    F_now_Voltage = 0;
end
% convert force level from Voltage to a percentage of MVC
F_now = (F_now_Voltage/MVC)*100;
% store force levels in the output
force_levels = [force_levels;...
    [F_now, timeCheck, F_now_Voltage, sampleOk_tmp]]; % store F in % of MVC, time and F in Volts

while (trial_success == 0) &&...
        ( (time_limit == false) ||...
        ( (time_limit == true) && (timeNow <= (onsetEffortPhase + t_max_effort)) ))
    % you either stop if the trial was successful (both learning and actual
    % task) OR if there is a time_limit and this time_limit was reached (actual task only)
    
    %% read current level of force
    [F_now_Voltage, timeNow, sampleOk_tmp] = F_read(dq, t_readWait);
    % update F_now only when the force sample was read properly (otherwise
    % keep previous value)
    if sampleOk_tmp == 1
        % convert force level from Voltage to a percentage of MVC
        F_now = (F_now_Voltage/MVC)*100;
    end
    % store force levels in the output
    force_levels = [force_levels;...
        [F_now, timeNow, F_now_Voltage, sampleOk_tmp]]; % store F in % of MVC, time and F in Volts
    
    %% update the center display according to if force above or below the threshold
    if F_now >= (F_threshold - F_tolerance) % while force > threshold, update the timer
        
        % update the effort index
        percSqueeze = percSqueeze + (force_levels(end,2) - force_levels(end-1,2))/t_effort_to_keep;
        
        if stateSqueezeON == false % the participant was not squeezing above threshold => new start
            n_starts = length(onsets.effort.start);
            onsets.effort.start(n_starts + 1) = timeNow;
            % update state of squeeze
            stateSqueezeON = true;
            % for the first occurrence update the indicator that the force
            % has been reached at least once
            if n_starts == 0
                has_F_threshold_ever_been_reached = 1;
            end
        end
        
    else % force below threshold
        % update the stop index only when angle is higher than zero
%         if currentAngle > startAngle
%             percStoppedSqueezing = percStoppedSqueezing + (timeNow - force_levels(end-1,2))/t_effort_to_keep;
%         end
        
        % switch from squeeze above threshold to squeeze below threshold
        if stateSqueezeON == true % the participant was squeezing above threshold but now he squeezes below threshold
            n_stops = length(onsets.effort.stop);
            onsets.effort.stop(n_stops + 1) = timeNow;
            % update state of squeeze
            stateSqueezeON = false;
        end
    end % Force above or below threshold?
    
    % update the angle for display
    if has_F_threshold_ever_been_reached == 1 % no need to update the angle as long as the threshold has not been reached at least once
        percentageTimeForceAlreadyMaintained = percSqueeze - percStoppedSqueezing;
        if percentageTimeForceAlreadyMaintained >= 0 && percentageTimeForceAlreadyMaintained <= 1
            currentAngle = startAngle + totalAngleDistance*percentageTimeForceAlreadyMaintained;
        elseif percentageTimeForceAlreadyMaintained < 0 % bound the angle to the start of the trial
            currentAngle = startAngle;
        elseif percentageTimeForceAlreadyMaintained > 1 % bound the top for visual display (otherwise ugly)
            currentAngle = endAngle;
        end
    end

    %% display on screen accordingly
    
    % display real-time force level
    disp_realtime_force(scr, F_threshold, F_tolerance, F_now, 'task');
    
    % display performance achieved level on top of the reward as an arc
    Screen('FillArc', window,...
        stim.difficulty.currLevelColor,...
        stim.difficulty.middle_center,...
        currentAngle,...
        endAngle - currentAngle);
    
%     [~,timeDispNow] = Screen('Flip',window);
    [lastFrameTime, timeDispNow]  = Screen('Flip', window, lastFrameTime + (0.5*ifi));
    
    % record effort display informations (timing and angle + force level
    dispInfos.time          = [dispInfos.time, timeDispNow];
    dispInfos.currentAngle  = [dispInfos.currentAngle, currentAngle];
    dispInfos.forceLevel    = [dispInfos.forceLevel, F_now];
    
    %% check whether performance was or not achieved
    if currentAngle >= endAngle
        trial_success = 1;
        onsets.effort_success = timeNow;
    end
end % time loop

%% stop acquisition of biopac handgrip
% stop acquiring data in the grip buffer
stop(dq);
% empty the grip buffer
flush(dq);

%% compute gains
physicalE_perf.performance = 100*((currentAngle - startAngle)/(360 - startAngle));

%% record vars of interest
physicalE_perf.trial_success = trial_success;
% performance is the amount in percentage of completion of the trial
physicalE_perf.performance = ((currentAngle-startAngle)/(360-startAngle))*100;
physicalE_perf.onsets = onsets;
% record all the force levels during the performance
physicalE_perf.t_max_effort = t_max_effort;
physicalE_perf.t_effort_to_keep = t_effort_to_keep;
physicalE_perf.force_levels = force_levels;
physicalE_perf.startAngle = startAngle;
physicalE_perf.finalAngle = currentAngle;
physicalE_perf.endAngle = endAngle;
physicalE_perf.has_F_threshold_ever_been_reached = has_F_threshold_ever_been_reached;
physicalE_perf.displayInformations = dispInfos;

end % function