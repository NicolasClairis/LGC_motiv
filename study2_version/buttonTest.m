function[]=buttonTest(key, window)
% function waiting for 2 button presses to move on
%
% INPUTS
% key: structure with code name for each key of the response button
%
% window: index for the window to use
%
% Nicolas Clairis

% wait all keys are released before starting
KbReleaseWait;

% start checking the keys
iPressed = 0;
DrawFormattedText(window,'Use one of the buttons to move on.','center','center'); % no bug experienced at this stage
Screen(window,'Flip');

% check key presses 
onset_question = GetSecs;
while iPressed < 2
    [keyisdown, ~, ~, lastPress, ~] = KbQueueCheck;
    
    [extremeLeftButtonPressed, middleLeftButtonPressed,...
        middleRightButtonPressed, extremeRightButtonPressed] = deal(false);
    if keyisdown == 1 % was a button pressed?
        
        % wait for participant to release the keys before updating the variables of interest
        % without this, it can create weird bugs as for the next question
        % it could reuse the previous answer
        KbReleaseWait();
        
        % one single button should be pressed at a time
        if lastPress(key.leftSure) > onset_question &&...
                lastPress(key.leftUnsure) < onset_question &&...
                lastPress(key.rightSure) < onset_question &&...
                lastPress(key.rightUnsure) < onset_question % extreme left answer (blue button)
            extremeLeftButtonPressed = true;
        elseif lastPress(key.leftSure) < onset_question &&...
                lastPress(key.leftUnsure) > onset_question &&...
                lastPress(key.rightSure) < onset_question &&...
                lastPress(key.rightUnsure) < onset_question % middle left answer (yellow button)
            middleLeftButtonPressed = true;
        elseif lastPress(key.leftSure) < onset_question &&...
                lastPress(key.leftUnsure) < onset_question &&...
                lastPress(key.rightSure) > onset_question &&...
                lastPress(key.rightUnsure) < onset_question % middle right answer (green button)
            middleRightButtonPressed = true;
        elseif lastPress(key.leftSure) < onset_question &&...
                lastPress(key.leftUnsure) < onset_question &&...
                lastPress(key.rightSure) < onset_question &&...
                lastPress(key.rightUnsure) > onset_question % extreme right answer (red button)
            extremeRightButtonPressed = true;
        end % which key was pressed
        
        if extremeLeftButtonPressed == true ||...
                middleLeftButtonPressed == true ||...
                middleRightButtonPressed == true ||...
                extremeRightButtonPressed == true
            iPressed = iPressed + 1;
        end
    end % was a key pressed?
end % wait for 2 buttons press before moving on

end % function