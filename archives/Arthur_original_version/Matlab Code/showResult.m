function totalMoney = showResult(scr,stim,speed,totalMoney,doWin,amountMoney)
switch stim.isPositiveStimuli
    case 1
        if doWin == 1
            Screen('TextSize', scr.window, 120);
            DrawFormattedText(scr.window, ['You won: ' num2str(amountMoney) ' Fr'], 'center',...
                scr.screenYpixels * 0.3, [0 0.8 0 ]);
            totalMoney = totalMoney + amountMoney;
        elseif doWin == 0
            Screen('TextSize', scr.window, 120);
            DrawFormattedText(scr.window, ['You missed: ' num2str(amountMoney) ' Fr'], 'center',...
                scr.screenYpixels * 0.3, [1 0.5 0 ]);
        end
        
        % Draw text in the bottom of the screen in Times in blue
        Screen('TextSize', scr.window, 100);
        DrawFormattedText(scr.window, strcat('Your total amount of money is\n ',num2str(totalMoney),' Fr'), 'center',...
            scr.screenYpixels * 0.7, 1);
        
    case 0
        
        % Draw text in the upper portion of the screen with the default font in red
        if doWin == 1
            Screen('TextSize', scr.window, 120);
            DrawFormattedText(scr.window, [ 'You saved: ' num2str(amountMoney) ' Fr'], 'center',...
                scr.screenYpixels * 0.3, [1 0.5 0 ]);
            totalMoney = totalMoney + amountMoney;
        elseif doWin == 0
            Screen('TextSize', scr.window, 120);
            DrawFormattedText(scr.window, ['You lost: ' num2str(amountMoney) ' Fr'], 'center',...
                scr.screenYpixels * 0.3, [0.8 0 0 ]);
        end
        
        % Draw text in the bottom of the screen in Times in blue
        Screen('TextSize', scr.window, 100);
        DrawFormattedText(scr.window, strcat('Your total amount of money is\n ',num2str(totalMoney),' Fr'), 'center',...
            scr.screenYpixels * 0.7, 1);
        
end


% Flip to the screen
speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);

% Now we have drawn to the screen we wait for a keyboard button press (any
% key) to terminate the demo
pause(1)





end

