function [measuredMVC,signal] = measureMVC(scr,stim,speed)

% Reset some variable locally
stim.VCsignal = 0;
signal = [];
stim.centeredCore = CenterRectOnPointd(stim.threshold_1,scr.xCenter, scr.yCenter);

% Set screen text size
Screen('TextSize', scr.window, 140);

%Measure MVC 3 times
for i = 1:3
    
    % Draw text and circle
    DrawFormattedText(scr.window, 'Ready', 'center', scr.screenYpixels * 0.8, [0 0.8 0 ]);
    Screen('FillOval', scr.window, [0.8 0 0 ], stim.centeredCore);
    speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
    pause(3)
    DrawFormattedText(scr.window, 'GO !', 'center', scr.screenYpixels * 0.8, 1);
    speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
    
    % During 4 second, show signal power and record it
    for frame_i = 1:60*4
        
        % Read the signal
        stim.VCsignal = fread(stim.u_out,1,'double');
        signal(i,frame_i) = stim.VCsignal;
        
        % Make a base Rect of 100 pixel width
        stim.baseRect = [scr.xCenter-50 scr.yCenter-85-stim.VCsignal*4 scr.xCenter+50 scr.yCenter-85];
        
        % We set bounds to make sure our square doesn't go completely off of the screen
        if stim.baseRect < 0
            stim.baseRect = 0;
        end
        if scr.squareY < 0
            scr.squareY = 0;
        elseif scr.squareY > scr.screenYpixels
            scr.squareY = scr.screenYpixels;
        end
        
        % Draw the rectangle and the oval with the text
        Screen('FillOval', scr.window, [0.8 0 0], stim.centeredCore);
        Screen('FillRect', scr.window, [0.8 0 0], stim.baseRect);
        DrawFormattedText(scr.window, 'GO !', 'center', scr.screenYpixels * 0.8,1);
        
        % Flip to the screen
        speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
        
    end
    
    % Show a rest text and give some rest
    DrawFormattedText(scr.window, 'Rest', 'center', scr.screenYpixels * 0.8, [0 0.8 0 ]);
    speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
    measuredMVC(i) = max(signal(i,:))
    pause(7)
    
end

end

