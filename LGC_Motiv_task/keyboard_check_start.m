function[] = keyboard_check_start(key, IRM)
%[] = keyboard_check_start(key, IRM)
% keyboard_check_start starts recording key presses (and TTL inputs)
%
% INPUTS
% key: structure with subfield with the corresponding key code for left and
% right key presses
%
% IRM:
% (0) no MRI: no need to check for TTL
% (1) MRI: wait for first TTL trigger from fMRI + record all the TTL sent by
% the scanner
%
% See also keyboard_check_end.m

%% empty the buffer in case that KbQueue was not closed properly
KbQueueRelease;

%% initialize keyboard
keysOfInterest = zeros(1,256);

%% record all subsequent TTL in the whole task
if IRM == 1
    keysOfInterest(key.trigger_id) = 1; % check TTL
end
% check also all relevant keyboard presses
keysOfInterest(key.left) = 1;
keysOfInterest(key.right) = 1;
if key.n_buttonsChoice == 4
    keysOfInterest(key.leftSure) = 1;
    keysOfInterest(key.leftUnsure) = 1;
    keysOfInterest(key.rightUnsure) = 1;
    keysOfInterest(key.rightSure) = 1;
end
% record space for training (but useless for fMRI)
if IRM == 0
    keysOfInterest(key.space) = 1;
end
% record escape in case experiment needs to be crashed
keysOfInterest(key.escape) = 1;

%% create buffer with keys to check
KbQueueCreate(0,keysOfInterest); % checks TTL and keys of pad
%% start filling the buffer
KbQueueStart; % starts checking

end % function