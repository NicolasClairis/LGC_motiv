function[] = disp_realtime_force(scr, F_threshold, F_tolerance, F_now)
%[] = disp_realtime_force(scr, F_threshold, F_tolerance, F_now)
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
% F_now: actual level of force
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
leftScaleLimit      = xScreenCenter*(1/4); % left limit of the scale
rightScaleLimit     = xScreenCenter*(2/4); % right limit of the scale
bottomScaleLimit    = yScreenCenter*(3/2); % bottom limit of the scale
topScaleLimit       = yScreenCenter*(1/2); % upper limit of the scale
graphYSize = bottomScaleLimit - topScaleLimit;
graphXSize = rightScaleLimit - leftScaleLimit;
% screen coordinates for bar at the center
leftBarLimit = leftScaleLimit + graphXSize*(1/4);
rightBarLimit = leftScaleLimit + graphXSize*(3/4);
% size and coordinates of half of the effort scale
yMetrics = yScreenCenter/2;
% distance between graduations
bigGrad = yMetrics/5;
smallGrad = bigGrad/4;

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
yFmaxUntilNow = bottomScaleLimit - graphYSize*(F_threshold/100);
Screen('DrawLine', window, red, leftScaleLimit, yFmaxUntilNow, rightScaleLimit, yFmaxUntilNow,3);

%% draw an orange bar with the actual level of force
yActualLevelBottom = bottomScaleLimit + 10;
yActualLevelTop = bottomScaleLimit - (F_now/100)*graphYSize;
Screen('FillRect', window, orange,...
    [leftBarLimit,...
    yActualLevelTop,...
    rightBarLimit,...
    yActualLevelBottom]);

%% display on screen
Screen(window,'Flip');
end % function