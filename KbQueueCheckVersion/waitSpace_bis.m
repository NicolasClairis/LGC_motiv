function [timeDispWait, timePress] = waitSpace_bis(scr, stim, window, keys)
% [timeDispWait, timePress] = waitSpace_bis(scr, stim, window, keys)
% function to wait for space press before moving on
%
% INPUTS
% scr: screen colours, etc
%
% window: PTB window index
%
% stim: structure with relevant key information
%
% keys: code for key space stored in this structure
%
% OUTPUTS
% timeDispWait: time when message to wait for space press is displayed on
% screen
%
% timePress: time of space press

%% display wait for space on screen
DrawFormattedText(window,...
            stim.pressSpace.text,...
            stim.pressSpace.x, stim.pressSpace.y,...
            scr.colours.white, scr.wrapat);
[~, timeDispWait] = Screen(window,'Flip');
disp('Please press space');

%% wait for space press
% [~, timePress, keyCode] = KbCheck();
[~, ~, ~, lastPress, ~] = KbQueueCheck;
while(lastPress(keys.space) < timeDispWait)
    % wait until the key has been pressed
    [~, ~, ~, lastPress, ~] = KbQueueCheck;
end

timePress = lastPress(keys.space);

end % function