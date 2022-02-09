function[TTL, keyLeft, keyRight,...
    keyLeftUnsure, keyLeftSure, keyRightUnsure, keyRightSure] = keyboard_check_end(key)
%[TTL, keyLeft, keyRight,...
%     keyLeftUnsure, keyLeftSure, keyRightUnsure, keyRightSure] = keyboard_check_end(TTL, trigger_id)
%
% INPUTS
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
% keyLeftSure, keyLeftUnsure: timing for all presses of left sure and left
% unsure keys
%
% keyRightSure, keyRightUnsure: timing for all presses of right sure and right
% unsure keys
%
% See also keyboard_check_start.m

%% stop recording key presses
KbQueueStop;

%% initialize keys of interest to store
keyLeft.Start = []; % time when starts pressing left key
keyLeft.Release = []; % time when releases right key
keyRight.Start = []; % time when starts pressing right key
keyRight.Release = []; % time when releases right key

% for quadratic case:
keyLeftUnsure.Start = [];
keyLeftUnsure.Release = [];
keyLeftSure.Start = [];
keyLeftSure.Release = [];
keyRightUnsure.Start = [];
keyRightUnsure.Release = [];
keyRightSure.Start = [];
keyRightSure.Release = [];
TTL = [];
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
        
    elseif key.n_buttonsChoice == 4 && event.Keycode == key.leftSure
        if event.Pressed == 1 % record start of press
            keyLeftSure.Start = [keyLeftSure.Start; event.Time];
        elseif event.Pressed == 0 % record time when release
            keyLeftSure.Release = [keyLeftSure.Release; event.Time];
        end
        
    elseif key.n_buttonsChoice == 4 && event.Keycode == key.leftUnsure
        if event.Pressed == 1 % record start of press
            keyLeftUnsure.Start = [keyLeftUnsure.Start; event.Time];
        elseif event.Pressed == 0 % record time when release
            keyLeftUnsure.Release = [keyLeftUnsure.Release; event.Time];
        end
        
    elseif key.n_buttonsChoice == 4 && event.Keycode == key.rightUnsure
        if event.Pressed == 1 % record start of press
            keyRightUnsure.Start = [keyRightUnsure.Start; event.Time];
        elseif event.Pressed == 0 % record time when release
            keyRightUnsure.Release = [keyRightUnsure.Release; event.Time];
        end
        
    elseif key.n_buttonsChoice == 4 && event.Keycode == key.rightSure
        if event.Pressed == 1 % record start of press
            keyRightSure.Start = [keyRightSure.Start; event.Time];
        elseif event.Pressed == 0 % record time when release
            keyRightSure.Release = [keyRightSure.Release; event.Time];
        end
    end
end

%% stop KbQueueCreate and clear cache
KbQueueRelease;

end % function