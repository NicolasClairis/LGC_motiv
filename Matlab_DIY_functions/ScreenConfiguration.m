function[x, y, window, baselineTextSize] = ScreenConfiguration(IRM, testing_script)
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

screens = Screen('Screens');
if IRM == 0
    whichScreen = max(screens);
elseif IRM == 1
    if testing_script == 0 % MRI
        whichScreen = 1; % 1 if 2 screens, 0 if one screen
    elseif testing_script == 1 % my own computer
        whichScreen = max(screens);
    end
end

Screen('Preference','VisualDebugLevel', 1); % avoid initial Psychtoolbox window
if testing_script == 0 % MRI
    Screen('Preference', 'SkipSyncTests', 0); % needs all other processes shut off
    window = Screen('OpenWindow',whichScreen,[0 0 0]);
elseif testing_script == 1 % my own computer
    Screen('Preference', 'SkipSyncTests', 1); % can work even if other softwares are on but displays an ugly red triangle at start
    % when testing, set resolution equal to MRI computer for testing
%     L_MRI_Screen = 1024;
%     H_MRI_Screen = 768;
%     Screen('Resolution', whichScreen, L_MRI_Screen, H_MRI_Screen); % sets a new whichScreen resolution for PTB
%     window = Screen('OpenWindow',whichScreen,[0 0 0], [1, 1, L_MRI_Screen, H_MRI_Screen]); % default window is full whichScreen screen

% window = Screen('OpenWindow', whichScreen, [], [1, 1, 800, 600]); % for debugging
    window = Screen('OpenWindow',whichScreen,[0 0 0]);
end
HideCursor()

baselineTextSize = 40;
Screen('TextSize', window, baselineTextSize);
Screen('TextFont', window, 'arial');
% if testing_script == 0 % MRI
[L, H] = Screen('WindowSize',whichScreen);
% elseif testing_script == 1
%        L = L_MRI_Screen;
%        H = H_MRI_Screen;
% end
x = L/2;
y = H/2;
% nber of characters authorized before breaking the line:
% wrapat_nb_char = 25;

end