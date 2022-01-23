function[was_a_key_pressed_bf_trial,...
    onsets_keyReleaseMessage,...
    dur_keyReleaseMessage] = check_keys_are_up_bis(scr, stim, key)
% [was_a_key_pressed_bf_trial,...
%     onsets_keyReleaseMessage,...
%     dur_keyReleaseMessage] = check_keys_are_up_bis(scr, stim, key)
% check_keys_are_up checks whether all relevant keys are up before
% starting the trial. If one of the relevant keys was being pressed before
% starting, displays an error message and waits for the participant to
% release all relevant keys before continuing the task.
%
% This script is like check_keys_are_up.m but it uses KbQueueCheck instead
% of KbCheck
%
% INPUTS
% scr: structure with display parameters
%
% stim: structure with stimuli informations
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
%
% dur_keyReleaseMessage: duration of the display of the error message (if
% displayed)

%% by default no key was pressed before the trial
was_a_key_pressed_bf_trial = 0;
[onsets_keyReleaseMessage, dur_keyReleaseMessage] = deal(NaN);

window = scr.window;
%% check key presses
% [keyIsDown, ~, keyCode] = KbCheck();
[keyIsDown, ~, ~, lastPress, lastRelease] = KbQueueCheck;

%% if one of the relevant keys is being pressed => display release message
% while (keyIsDown == 1) &&...
%         
while (keyIsDown == 1) &&...
        (key.n_buttonsChoice == 2 &&...
        ( (lastRelease(key.left) < lastPress(key.left)) ||...
        (lastRelease(key.right) < lastPress(key.right)) )) ||...
        (key.n_buttonsChoice == 4 &&...
        ( (lastRelease(key.leftSure) < lastPress(key.leftSure)) ||...
        (lastRelease(key.leftUnsure) < lastPress(key.leftUnsure)) ||...
        (lastRelease(key.rightUnsure) < lastPress(key.rightUnsure)) ||...
        (lastRelease(key.rightSure) < lastPress(key.rightSure)) ) )
    was_a_key_pressed_bf_trial = 1;
    DrawFormattedText(window, stim.releaseButtonsMsg.text,...
        stim.releaseButtonsMsg.x, stim.releaseButtonsMsg.y,...
        stim.releaseButtonsMsg.colour);
    [~, onsets_keyReleaseMessage] = Screen(window,'Flip');
    % keep checking the buttons to know if the keyboard has been released
    % or not
    [keyIsDown, ~, ~, lastPress, lastRelease] = KbQueueCheck;
    dur_keyReleaseMessage = GetSecs - onsets_keyReleaseMessage;
end % some relevant key is being pressed

end % function