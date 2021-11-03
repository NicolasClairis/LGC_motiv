function[scr, xScreenCenter, yScreenCenter, window, baselineTextSize] = ScreenConfiguration(IRM, testing_script)
% [scr, xScreenCenter, yScreenCenter, window, baselineTextSize] = ScreenConfiguration(IRM, testing_script)
% function with common parameters for starting psychtoolbox for any of the
% three tasks used in fMRI (taskGripRP, taskMentalRP and taskLearning75).
% It could be reused by any other task.
%
% INPUTS
% IRM: is it for training (outside fMRI: only one screen) (0) or inside the
% scanner (should have 2 screens: one for us and one where you can follow
% what the subject sees in the scanner. PTB has to be opened in the latter)
% testing_script: precise whether you don't care at all about timings (1)
% or if you are actually testing a real subject and thus you care about it
% (0)
%
% OUTPUTS
% scr: structure with main screen informations (screen window, screen size,
% x and y center coordinates, etc.)
%
% xScreenCenter,yScreenCenter: x and y coordinates of the center of the screen
%
% window: window where PTB stims are displayed
%
% baselineTextSize: baseline size of the text displayed on the screen
%
% Developed by Nicolas Clairis - february 2017

%% select on which screen PTB will be displayed
% (particularly when you have several screens, 
% like at the MRI scanner where you need to be sure that things are 
% displayed on the scanner screen)
screens = Screen('Screens');
if IRM == 0
    whichScreen = max(screens);
elseif IRM == 1 % CIBM computer
    if testing_script == 0 % for fMRI experiment, display on the projector
        whichScreen = max(screens); % 1 if 2 screens, 0 if one screen
    elseif testing_script == 2 % for fMRI training display on screen 1
        whichScreen = 1;
    elseif testing_script == 1 % for debugging
        whichScreen = 1;
    end
end

%% set screen colour
black = [0 0 0];
white = [255 255 255];
orange = [255 153 0];
red = [255 0 0];
weakRed = [100 0 0];
grey = [128 128 128];
screenColour = grey;

%% open PTB window + set debug parameters

% remove initial Psychtoolbox welcome and warning window
Screen('Preference','VisualDebugLevel', 0);

% require the computer to have best timing performance only when not
% testing
switch testing_script
    case {0,2} % CIBM
        Screen('Preference', 'SkipSyncTests', 0); % needs all other processes shut off
    case 1 % my own computer
        Screen('Preference', 'SkipSyncTests', 1); % can work even if other softwares are on but displays an ugly red triangle at start
end

% open PTB window
window = Screen('OpenWindow',whichScreen,screenColour);

%% hide mouse cursor
HideCursor();

%% extract x and y coordinates of the center of the screen
[L, H] = Screen('WindowSize',whichScreen);

%% for fMRI CIBM shitty screen, fix the location of your stimuli
if IRM == 0
    leftBorder = 0;
    upperBorder = 0;
    rightBorder = L;
    lowerBorder = H;
elseif IRM == 1
    leftBorder = 0;
    upperBorder = 160;
    rightBorder = L;
    lowerBorder = 934;
end
xScreenCenter = leftBorder + ( rightBorder - leftBorder )/2;
yScreenCenter = upperBorder + ( lowerBorder - upperBorder )/2;
visibleYsize = lowerBorder - upperBorder;
visibleXsize = rightBorder - leftBorder;
visibleWindowRectCoord = [leftBorder upperBorder rightBorder lowerBorder];

%% text display properties
% rescaling factor for shitty CIBM scanner display screen
if IRM == 0
    rescaleCIBM = 1;
elseif IRM == 1
    rescaleCIBM_x = (rightBorder - leftBorder)/L;
    rescaleCIBM_y = (lowerBorder - upperBorder)/H;
    rescaleCIBM = max(rescaleCIBM_x, rescaleCIBM_y);
end
baselineTextSize = round(50*rescaleCIBM);
Screen('TextSize', window, baselineTextSize);
Screen('TextFont', window, 'arial');
textSize.baseline = baselineTextSize;
textSize.mentalNumber = round(130*rescaleCIBM);
textSize.big = round(160*rescaleCIBM);
textSize.middle = round(110*rescaleCIBM);
textSize.reward = round(80*rescaleCIBM);
textSize.taskPeriodsTitles = round(80*rescaleCIBM);

%% store main informations inside scr structure
scr.realXsize = L;
scr.realYsize = H;
scr.rescaleCIBM = rescaleCIBM;
scr.screenNumber = whichScreen;
scr.textSize = textSize;
scr.window = window;
scr.xCenter = xScreenCenter;
scr.yCenter = yScreenCenter;
scr.leftBorder = leftBorder;
scr.lowerBorder = lowerBorder;
scr.upperBorder = upperBorder;
scr.rightBorder = rightBorder;
scr.visibleYsize = visibleYsize;
scr.visibleXsize = visibleXsize;
scr.visibleWindow = visibleWindowRectCoord;
scr.background_colour = screenColour;
scr.colours.grey = grey;
scr.colours.white = white;
scr.colours.black = black;
scr.colours.orange = orange;
scr.colours.red     = red;
scr.colours.weakRed = weakRed;
scr.wrapat = 50; % limit characters for drawformattedtext
end % function