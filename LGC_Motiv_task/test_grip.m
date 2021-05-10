[scr, xScreenCenter, yScreenCenter, window, baselineTextSize] = ScreenConfiguration(0, 1);
[key, dq] = relevant_key_definition('physical', 0);
[stim] = stim_initialize(scr, 3, 3, 'D:\LGC_motiv\LGC_Motiv_task\Coin_PNG\');

start(dq,"continuous");
pause(0.125)
% will need data = read(dq) function only to read the signal

MVC=3.5;
F_threshold = 50;
t_effort_to_keep=3;
F_tolerance=5;
nFrames=6;
ifi = Screen('GetFlipInterval', window);
startAngle = stim.difficulty.startAngle.(['level_',num2str(2)]);
currentAngle = startAngle;
endAngle = stim.difficulty.arcEndAngle;
totalAngleDistance = endAngle - startAngle;
%% display all relevant variables on the screen for the effort period

% display effort scale on the left for the live display of the effort level
% add force threshold to maintain
force_initialDisplay = 0;
disp_realtime_force(scr, F_threshold, F_tolerance, force_initialDisplay, 'task');

% display difficulty level on top of the reward as an arc
lastFrameTime=Screen('Flip',window);
trial_success=false;
%% effort
[force_levels,...
    dispInfos.time,...
    dispInfos.currentAngle,...
    dispInfos.forceLevel] = deal([]);
force_levels = [0,NaN,0];
[timeEffortStart_tmp, timeEffortStop_tmp] = deal(0);
stateSqueezeON = false;
[onsets.effort.start, onsets.effort.stop] = deal([]);

% initialize read
timeCheck = GetSecs;
timeStart = timeCheck;
F_now_Voltage_tmp = read(dq,'all','OutputFormat','Matrix');
F_now_Voltage = F_now_Voltage_tmp(end);
% convert force level from Voltage to a percentage of MVC
F_now = (F_now_Voltage/MVC)*100;
F_now_interpolate_tmp = linspace(force_levels(end, 1), F_now, nFrames);
% store force levels in the output
force_levels = [force_levels;...
    [F_now, timeCheck, F_now_Voltage]]; % store F in % of MVC, time and F in Volts
% short pause to avoid empty matrix reading by read.m function
pause(0.1);
iFrame = 1;
jTest=0;
while (trial_success == false) || (timeNow < (timeStart + 5))
% for j=1:1000
    % you either stop if the trial was successful (both learning and actual
    % task) OR if there is a time_limit defined and this time_limit was
    % reached (actual task only)
    
    %% read current level of force
    % check timing now when we extract the level of force
    timeNow = GetSecs;
    a = timeNow -(timeStart + 5)
    % check force being applied now
    if timeNow >= (timeCheck + 0.3)
        F_now_Voltage_tmp = read(dq,'all','OutputFormat','Matrix');
        if ~isempty(F_now_Voltage_tmp)
            F_now_Voltage = F_now_Voltage_tmp(end);
            timeCheck = timeNow;
            iFrame = 1;
            jTest = jTest + 1;
            disp(jTest);
        end
        % convert force level from Voltage to a percentage of MVC
        F_now = (F_now_Voltage/MVC)*100;
        F_now_interpolate_tmp = linspace(force_levels(end, 1), F_now, nFrames);
        % store force levels in the output
        force_levels = [force_levels;...
            [F_now, timeNow, F_now_Voltage]]; % store F in % of MVC, time and F in Volts
        % reset frame for visual display
    end
    
    %% update the center display according to if force above or below the threshold
    if F_now >= (F_threshold - F_tolerance) % while force > threshold, update the timer
        
        if stateSqueezeON == false % the participant was not squeezing above threshold => new start
            n_starts = length(onsets.effort.start);
            onsets.effort.start(n_starts + 1) = timeNow;
            % update state of squeeze
            stateSqueezeON = true;
        end
        
        % update the angle for the display
        if ~isempty(onsets.effort.start)
            timeEffortStart_tmp = onsets.effort.start(end);
        end
        if ~isempty(onsets.effort.stop)
            timeEffortStop_tmp = onsets.effort.stop(end);
        end
        
    else % force below threshold
        % switch from squeeze above threshold to squeeze below threshold
        if stateSqueezeON == true % the participant was squeezing above threshold but now he squeezes below threshold
            n_stops = length(onsets.effort.stop);
            onsets.effort.stop(n_stops + 1) = timeNow;
            % update state of squeeze
            stateSqueezeON = false;
        end
        
        % update the angle for the display
        if ~isempty(onsets.effort.start)
            timeEffortStart_tmp = onsets.effort.start(end);
        end
        if ~isempty(onsets.effort.stop)
            timeEffortStop_tmp = onsets.effort.stop(end);
        end
    end % Force above or below threshold?
    
    % update the angle for display
    if timeEffortStart_tmp > 0
        percentageTimeForceAlreadyMaintained = (timeNow - timeEffortStart_tmp + timeEffortStop_tmp)/t_effort_to_keep;
        if percentageTimeForceAlreadyMaintained >= 0
            currentAngle = totalAngleDistance*percentageTimeForceAlreadyMaintained;
        else
            currentAngle = startAngle;
        end
    end
    
    %% display on screen accordingly
    
    % display real-time force level
    disp_realtime_force(scr, F_threshold, F_tolerance, F_now_interpolate_tmp(iFrame), 'task');
    if iFrame < nFrames
        iFrame = iFrame + 1;
    end
    
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
        trial_success = true;
        onsets.effort_success = timeNow;
    end
end % time loop


%% stop acquisition of biopac handgrip
stop(dq);
sca;