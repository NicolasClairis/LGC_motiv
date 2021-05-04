function [speed] = LGCM_stim_speed(window)
%[speed] = LGCM_stim_speed(window)
% LGCM_stim_speed
%
% INPUTS
% window: identification of the PTB window
%
% OUTPUTS
% speed: structure with main informations about stimulus speed display
%
% See also LGCM_main_experiment.m

%% Query the frame duration
speed.ifi = Screen('GetFlipInterval', window);

%% Prepare Rest interval in seconds
speed.isShortBreak = true;
speed.shortBreak = 10;
speed.longBreak = 60;

%% Sync us and get a time stamp
speed.vbl = Screen('Flip', window);
speed.waitframes = 1;
% Set the amount we want our square to move on each button press
speed.pixelsPerFrame = 6;

end % function