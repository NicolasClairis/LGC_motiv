[scr, xScreenCenter, yScreenCenter, window, baselineTextSize] = ScreenConfiguration(0, 1);

[key] = relevant_key_definition('physical', 0, 2);
R_money = R_amounts(3,'yes');
[stim] = stim_initialize(scr, 3, 'fr', R_money);


%% initialize BHAPI
% initialize the grip
samplingRate = 100; % = how many samples per second
[libname] = loading_grip(samplingRate);
% launch acquisition Daemon
retval = calllib(libname, 'startMPAcqDaemon');

if ~strcmp(retval,'MPSUCCESS')
    fprintf(1,'Failed to Start Acquisition Daemon.\n');
    calllib(libname, 'disconnectMPDev');
    sca;
    return
end

% launch actual acquisition
retval = 'MPHasToStart';
z = 0;
nRestartTry = 25;
while ~strcmp(retval,'MPSUCCESS') && z < nRestartTry
    retval = calllib(libname, 'startAcquisition');
    z = z + 1;
end
if ~strcmp(retval,'MPSUCCESS')
    fprintf(1, ['Failed to start acquisition for grip measurement after ',num2str(nRestartTry),' tries.\n']);
    calllib(libname, 'disconnectMPDev');
    sca;
    return
end
pause(0.125);

%% general parameters
MVC = 0.0001;
F_threshold = 50;
t_effort_to_keep = 3;
F_tolerance = 3;
nFrames = 6;
% ifi = Screen('GetFlipInterval', window);
ifi = 1/15;
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
lastFrameTime = Screen('Flip',window);
trial_success = false;
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
F_now_Voltage_tmp = read_bhapi_grip(libname);
F_now_Voltage = F_now_Voltage_tmp(end);
% convert force level from Voltage to a percentage of MVC
F_now = (F_now_Voltage/MVC)*100;
% F_now_interpolate_tmp = linspace(force_levels(end, 1), F_now, nFrames);
% store force levels in the output
force_levels = [force_levels;...
    [F_now, timeCheck, F_now_Voltage]]; % store F in % of MVC, time and F in Volts
% short pause to avoid empty matrix reading by read.m function
pause(0.1);
iFrame = 1;
jTest=0;
z=0;
[percSqueeze, percStoppedSqueeze] = deal(0);
timeNow = GetSecs;
while (trial_success == false) && (timeNow < (timeStart + 5))
    z = z+ 1;
    % for j=1:1000
    % you either stop if the trial was successful (both learning and actual
    % task) OR if there is a time_limit defined and this time_limit was
    % reached (actual task only)
    
    %% read current level of force
    % check timing now when we extract the level of force
    timeNow = GetSecs;
    a = timeNow -(timeStart + 5);
    % check force being applied now
    %     if timeNow >= (timeCheck + 0.3)
    F_now_Voltage_tmp = read_bhapi_grip(libname);
    pause(0.075)
    if ~isempty(F_now_Voltage_tmp)
        F_now_Voltage = F_now_Voltage_tmp(end);
        timeCheck = timeNow;
        iFrame = 1;
        jTest = jTest + 1;
        disp(jTest);
    end
    % convert force level from Voltage to a percentage of MVC
    F_now = (F_now_Voltage/MVC)*100;
    %         F_now_interpolate_tmp = linspace(force_levels(end, 1), F_now, nFrames);
    % store force levels in the output
    force_levels = [force_levels;...
        [F_now, timeNow, F_now_Voltage]]; % store F in % of MVC, time and F in Volts
    % reset frame for visual display
    %     end
    
    %% update the center display according to if force above or below the threshold
    if F_now >= (F_threshold - F_tolerance) % while force > threshold, update the timer
        
        percSqueeze = percSqueeze + (force_levels(end,2) - force_levels(end-1,2))/t_effort_to_keep;
        
        if stateSqueezeON == false % the participant was not squeezing above threshold => new start
            n_starts = length(onsets.effort.start);
            onsets.effort.start(n_starts + 1) = timeNow;
            % update state of squeeze
            stateSqueezeON = true;
        end
        
    else % force below threshold
        
        % update the stop only when angle is higher than zero
        if currentAngle > startAngle
            percStoppedSqueeze = percStoppedSqueeze + (timeNow - force_levels(end-1,2))/t_effort_to_keep;
        end
        
        % switch from squeeze above threshold to squeeze below threshold
        if stateSqueezeON == true % the participant was squeezing above threshold but now he squeezes below threshold
            n_stops = length(onsets.effort.stop);
            onsets.effort.stop(n_stops + 1) = timeNow;
            % update state of squeeze
            stateSqueezeON = false;
        end
    end % Force above or below threshold?
    
    % update the angle for display
    if ~isempty(onsets.effort.start)
        percentageTimeForceAlreadyMaintained = percSqueeze - percStoppedSqueeze;
        if percentageTimeForceAlreadyMaintained >= 0 && percentageTimeForceAlreadyMaintained <= 1
            currentAngle = startAngle + totalAngleDistance*percentageTimeForceAlreadyMaintained;
        elseif percentageTimeForceAlreadyMaintained < 0
            currentAngle = startAngle;
        elseif percentageTimeForceAlreadyMaintained > 1
            currentAngle = endAngle;
        end
    end
    
    %% display on screen accordingly
    
    % display real-time force level
    disp_realtime_force(scr, F_threshold, F_tolerance, F_now, 'task');
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
% stop acquiring data
retval = calllib(libname, 'stopAcquisition');
if ~strcmp(retval,'MPSUCCESS')
    fprintf(1,'Failed to Stop\n');
    calllib(libname, 'disconnectMPDev');
    sca;
    return
end

% disconnect MPDEV
fprintf(1,'Disconnecting...\n')
retval = calllib(libname, 'disconnectMPDev');
% close PTB
sca;