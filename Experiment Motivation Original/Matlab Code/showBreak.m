function totalMoney = showBreak(scr,stim,speed,totalMoney)

% Prepare general text size
Screen('TextSize', scr.window, 130);

% If we want a short break or long break
switch speed.isShortBreak
    case 1
        
        % Show money participant owns and the resting time
        DrawFormattedText(scr.window, ['Rest for ' num2str(speed.shortBreak) ' seconds'], 'center',...
            scr.screenYpixels * 0.3, 1);
        DrawFormattedText(scr.window, strcat('Your total amount of money is\n ',num2str(totalMoney),' Fr'), 'center',...
            scr.screenYpixels * 0.7, 1);
        
        % Flip to the screen
        speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
        pause(speed.shortBreak)
        
        % Announce next type of block coming
        switch stim.isPositiveStimuli
            case 1
                DrawFormattedText(scr.window, 'Next up : MISS or WIN Money ', 'center',...
                    scr.screenYpixels * 0.7, 1);
            case 0
                DrawFormattedText(scr.window, 'Next up : SAVE or LOSE Money ', 'center',...
                    scr.screenYpixels * 0.7, 1);
            case 2
                DrawFormattedText(scr.window, 'Next up : Neutral Value ', 'center',...
                    scr.screenYpixels * 0.7, 1);
        end
        
    case 0
        % Show money participant owns and the resting time
        DrawFormattedText(scr.window, ['Rest for ' num2str(speed.longBreak/60) ' minutes'], 'center',...
            scr.screenYpixels * 0.3, 1);
        DrawFormattedText(scr.window, strcat('Your total amount of money is\n ',num2str(totalMoney),' Fr'), 'center',...
            scr.screenYpixels * 0.7, 1);
        
        % Flip to the screen
        speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
        pause(speed.longBreak)
        Screen('TextSize', scr.window, 100);
        DrawFormattedText(scr.window, 'Next up : Maximum Power ', 'center',...
            scr.screenYpixels * 0.7, 1);
        
end

%Prepare participant for incoming block
Screen('TextSize', scr.window, 150)
DrawFormattedText(scr.window, 'Get Ready', 'center', scr.screenYpixels * 0.3, 1);

% Flip to the screen
speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);

pause(5)





end

