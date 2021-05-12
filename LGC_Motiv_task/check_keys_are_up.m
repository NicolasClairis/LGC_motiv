function[was_a_key_pressed_bf_trial,...
    onsets_keyReleaseMessage] = LGCM_check_keys_are_up(scr, key)
% [was_a_key_pressed_bf_trial,...
%     onsets_keyReleaseMessage] = LGCM_check_keys_are_up(scr, key)
% LGCM_check_keys_are_up checks whether all relevant keys are up before
% starting the trial. If one of the relevant keys was being pressed before
% starting, displays an error message and waits for the participant to
% release all relevant keys before continuing the task.
%
% INPUTS
% scr: structure with display parameters
%
% key: structure with left and right key codes
%
% OUTPUTS
% was_a_key_pressed_bf_trial
% (0) no key was pressed before the start of the trial
% (1) one of the relevant keys was being pressed before the start of the
% trial
% 
% onsets_keyReleaseMessage: onset of the error message (if displayed),
% otherwise NaN value

%% by default no key was pressed before the trial
was_a_key_pressed_bf_trial = 0;
onsets_keyReleaseMessage = NaN;

%% check key presses
[keyIsDown, secs, keyCode] = KbCheck();

%% if one of the relevant keys is being pressed => display release message
while (keyIsDown == 1) &&...
        ( (keyCode(key.left) == 1) || (keyCode(key.right) == 1 ))
    was_a_key_pressed_bf_trial = 1;
    DrawFormattedText(scr.window, 'Relâchez les boutons svp','center','center')
    [~, onsets_keyReleaseMessage] = Screen(scr.window,'Flip');
    % keep checking the buttons to know if the keyboard has been released
    % or not
    [keyIsDown, secs, keyCode] = KbCheck();
end % some relevant key is being pressed

end % function