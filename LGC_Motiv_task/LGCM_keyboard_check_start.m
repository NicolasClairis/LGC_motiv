function[T0, TTL] = LGCM_keyboard_check_start(dummy_scan, trigger_id, key)
%[T0, TTL] = LGCM_keyboard_check_start(dummy_scan, trigger_id, key)
% LGCM_keyboard_check_start starts recording key presses (and TTL inputs)
%
% INPUTS
% dummy_scan: number of TTL to wait before starting the task
%
% trigger_id: number corresponding to the TTL trigger
%
% key: structure with subfield with the corresponding key code for left and
% right key presses
%
% OUTPUTS
% T0: time of the first TTL on which all onsets will be referenced
%
% TTL: vector which will serve to keep all TTL
%
% See also LGCM_keyboard_check_end.m

next = 0;
TTL = []; % TTL TIMES
% wait dummy_scan number of volumes before starting the task
while next < dummy_scan
    [keyisdown, T0IRM, keycode] = KbCheck;
    
    if keyisdown == 1 && keycode(trigger_id) == 1
        if next == 0
            T0 = T0IRM;
        end
        next = next+1;
        disp(next);
        TTL = [TTL; T0IRM];
        while keycode(trigger_id) == 1
            [keyisdown, T, keycode] = KbCheck;
        end
    end
end

%% record all subsequent TTL in the whole task
keysOfInterest = zeros(1,256);
keysOfInterest(trigger_id) = 1; % check TTL
keysOfInterest(key.left) = 1; % check all left key press
keysOfInterest(key.right) = 1; % check all right key press
KbQueueCreate(0,keysOfInterest); % checks TTL and keys of pad
KbQueueStart; % starts checking

end % function