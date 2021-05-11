function[] = disp_realtime_force(scr, F_threshold, F_tolerance, F_now, calib_or_Task)
%[] = disp_realtime_force(scr, F_threshold, F_tolerance, F_now, calib_or_Task)
% disp_realtime_force will display the force being exerted in real-time on
% the left of the screen. The top of the scale will correspond to the MVC
% of the participant. A red bar will be placed on the F_threshold that the
% participant needs to reach. F_tolerance allows some tolerance around this
% threshold.
%
% INPUTS
% scr: structure with screen informations
%
% F_threshold: threshold that the participants need to stay above for
% considering that the performance is ok
%
% F_tolerance: tolerated range around the threshold that we will still
% consider to be ok (maybe display with a similar colour?)
%
% F_now: actual level of force (expressed in percentage of MVC, not Voltage!)
%
% calib_or_Task: will display the scale on the left (to allow space for the
% oval) or on the center of the screen depending on whether it's the actual
% task or the calibration
% 'calib': calibration phase (display on the middle)
% 'task': task phase (display on the left)
%
% OUTPUTS
%

%% extract screen relevant parameters
window = scr.window;
xScreenCenter = scr.xCenter;
yScreenCenter = scr.yCenter;

%% color parameters
red = [255 0 0];
weakRed = [100 0 0];
white = [255 255 255];
orange = [255 153 0];

%% screen coordinates for effort scale
bottomScaleLimit    = yScreenCenter*(3/2); % bottom limit of the scale
topScaleLimit       = yScreenCenter*(1/2); % upper limit of the scale
graphYSize = bottomScaleLimit - topScaleLimit;
% size and coordinates of half of the effort scale
yMetrics = yScreenCenter/2;
% distance between graduations
bigGrad = yMetrics/5;
smallGrad = bigGrad/4;
switch calib_or_Task
    case 'calib' % center the scale
        leftScaleLimit      = xScreenCenter*(3.5/4); % left limit of the scale
        rightScaleLimit     = xScreenCenter*(4.5/4); % right limit of the scale
    case 'task' % put the scale on the left
        leftScaleLimit      = xScreenCenter*(1/4); % left limit of the scale
        rightScaleLimit     = xScreenCenter*(2/4); % right limit of the scale
end
graphXSize = rightScaleLimit - leftScaleLimit;
% screen coordinates for the bar representing the realtime force level
leftBarLimit = leftScaleLimit + graphXSize*(1/4);
rightBarLimit = leftScaleLimit + graphXSize*(3/4);

%% draw a line on the left of the scale (vertical bar)
% Screen('DrawLine', window, white, leftScaleLimit, topScaleLimit, leftScaleLimit, bottomScaleLimit, 3);

%% draw the scale (horizontal bars)
for yaxis = -yMetrics:smallGrad:yMetrics
    Screen('DrawLine', window, weakRed,...
        leftScaleLimit,...
        (yScreenCenter+yaxis),...
        rightScaleLimit,...
        (yScreenCenter+yaxis), 1);
end

%% ??
for yaxis = (-4*yMetrics/5):(yMetrics/5):yMetrics
    Screen('DrawLine', window, white,...
        leftScaleLimit,...
        (yScreenCenter+yaxis),...
        rightScaleLimit,...
        (yScreenCenter+yaxis), 3);
end

%% draw the threshold to reach

% % draw a red line to indicate the threshold to reach
% yThreshold = bottomScaleLimit - graphYSize*(F_threshold/100);
% Screen('DrawLine', window, red, leftScaleLimit, yThreshold, rightScaleLimit, yThreshold,3);

% draw a red rectangle representing the window (including the tolerance
% threshold) to reach to have a good performance
yLowThreshold = bottomScaleLimit - graphYSize*((F_threshold - F_tolerance)/100);
yTopThreshold = bottomScaleLimit - graphYSize*((F_threshold + F_tolerance)/100);
Screen('FillRect', window, red,...
    [leftScaleLimit,...
    yTopThreshold,...
    rightScaleLimit,...
    yLowThreshold]);

% %% you can also add the tolerance threshold but maybe better avoid it so
% that participant do not see this as the actual threshold to reach
% yThresholdTolerance = bottomScaleLimit - graphYSize*((F_threshold - F_tolerance)/100);
% Screen('DrawLine', window, weakRed, leftScaleLimit, yThresholdTolerance, rightScaleLimit, yThresholdTolerance,3);

%% draw an orange bar with the actual level of force
yActualLevelBottom = bottomScaleLimit + 10;
if F_now > 0 && F_now < 100
    yActualLevelTop = bottomScaleLimit - (F_now/100)*graphYSize;
elseif F_now <= 0 % bound to bottom of the scale
    yActualLevelTop = bottomScaleLimit;
elseif F_now >= 100 % bount to top of the scale
    yActualLevelTop = bottomScaleLimit - graphYSize;
end
Screen('FillRect', window, orange,...
    [leftBarLimit,...
    yActualLevelTop,...
    rightBarLimit,...
    yActualLevelBottom]);

end % function