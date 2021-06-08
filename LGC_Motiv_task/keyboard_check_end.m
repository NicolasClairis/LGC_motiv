function[TTL, keyLeft, keyRight] = keyboard_check_end(TTL, key)
%[TTL, keyLeft, keyRight] = keyboard_check_end(TTL, trigger_id)
%
% INPUTS
% TTL: vector with initial fMRI TTL timings
%
% key: structure with relevant keys
%   .trigger_id : number corresponding to the key associated to the TTL
%   .left
%   .right
%
% OUTPUTS
% TTL: TTL vector updated with all fMRI TTL timings
%
% keyLeft: timing for all presses of left key
%
% keyRight: timing for all presses of right key
%
% See also keyboard_check_start.m

%% stop recording key presses
KbQueueStop;

%% initialize keys of interest to store
keyLeft.Start = []; % time when starts pressing left key
keyLeft.Release = []; % time when releases right key
keyRight.Start = []; % time when starts pressing right key
keyRight.Release = []; % time when releases right key
%% release all keys and associated timings
while KbEventAvail
    [event, n] = KbEventGet;
    if event.Keycode == key.trigger_id
        TTL = [TTL; event.Time];
    elseif event.Keycode == key.left % if left key pressed
        if event.Pressed == 1 % record start of press
            keyLeft.Start = [keyLeft.Start; event.Time];
        elseif event.Pressed == 0 % record time when release
            keyLeft.Release = [keyLeft.Release; event.Time];
        end
    elseif event.Keycode == key.right % if right key pressed
        if event.Pressed == 1 % record start of press
            keyRight.Start = [keyRight.Start; event.Time];
        elseif event.Pressed == 0 % record time when release
            keyRight.Release = [keyRight.Release; event.Time];
        end
    end
end

%% stop KbQueueCreate and clear cache
KbQueueRelease;

end % function