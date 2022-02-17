function [timeDispWait, timePress] = waitSpace(scr, stim, window, keys)
% [timeDispWait, timePress] = waitSpace(scr, stim, window, keys)
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

% kbcheck version
% [~, timePress, keyCode] = KbCheck();
% while(keyCode(keys.space) ~= 1)
%     % wait until the key has been pressed
%     [~, timePress, keyCode] = KbCheck();
% end

% kbqueuecheck version
[pressed, ~, ~, lastPress, ~] = KbQueueCheck();
while (pressed == 0) ||...
        (lastPress(keys.space) <= timeDispWait)
    % wait until the key has been pressed
    [pressed, ~, ~, lastPress, ~] = KbQueueCheck();
end
% wait until key is released
KbReleaseWait;

end % function