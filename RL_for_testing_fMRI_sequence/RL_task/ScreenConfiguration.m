function[x, y, window, baselineTextSize] = ScreenConfiguration(screen_nb, testing_script)
% [x, y, window] = ScreenConfiguration(IRM, testing_script)
% function with common parameters for starting psychtoolbox for any of the
% three tasks used in fMRI (taskGripRP, taskMentalRP and taskLearning75).
% It could be reused by any other task.
%
% INPUTS:
% IRM: is it for training (outside fMRI: only one screen) (0) or inside the
% scanner (should have 2 screens: one for us and one where you can follow
% what the subject sees in the scanner. PTB has to be opened in the latter)
% testing_script: precise whether you don't care at all about timings (1)
% or if you are actually testing a real subject and thus you care about it
% (0)
%
% OUTPUTS
% x,y: x and y coordinates of the center of the screen
% window: window where PTB stims are displayed
%
% Developed by Nicolas Clairis - february 2017
% updated november 2020 for LGC experiment

%% screen colour
sc_colour = [0.5 0.5 0.5];

%% define screen where task will be displayed
screens = Screen('Screens');
if screen_nb == 1
    whichScreen = max(screens);
elseif screen_nb == 2
    if testing_script == 0
        whichScreen = 1; % 1 if 2 screens, 0 if one screen
    elseif testing_script == 1 % my own computer
        whichScreen = max(screens);
    end
end

%% avoid initial Psychtoolbox window
Screen('Preference','VisualDebugLevel', 1);

%% timings
if testing_script == 0
    % PTB priority option ON
    Screen('Preference', 'SkipSyncTests', 0); % needs all other processes shut off
    window = Screen('OpenWindow',whichScreen,sc_colour);
    
    % Windows Maximum priority level so that all the computer power is given to this
    topPriorityLevel = MaxPriority(window);
    Priority(topPriorityLevel);
    
elseif testing_script == 1 % my own computer
    Screen('Preference', 'SkipSyncTests', 1); % can work even if other softwares are on but displays an ugly red triangle at start
    window = Screen('OpenWindow',whichScreen,sc_colour);
end
HideCursor();

baselineTextSize = 40;
Screen('TextSize', window, baselineTextSize);
Screen('TextFont', window, 'arial');
[L, H] = Screen('WindowSize',whichScreen);

x = L/2;
y = H/2;


end