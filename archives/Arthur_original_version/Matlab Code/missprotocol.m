function [stim,speed,frame_i,signal] = missprotocol(scr,stim,speed,signal,frame_i)

% prepare sequence of when we have to draw the circles or not to create a flickering effect
showCircles = repmat([0 0 0 0 0 1 1 1 1 1],9);
showCircles = showCircles(1,:);

% The money will flicker too so we have to center it
stim.centeredMoney  = CenterRectOnPointd(stim.money,scr.xCenter, scr.yCenter);

% During 3 seconds
for i = 1:90
    
    % Continue keeping track of the signal, even if we don't show it
    stim.VCsignal = fread(stim.u_out,1,'double');
    
    %keep the signal to save it
    signal(frame_i) = stim.VCsignal;
    frame_i = frame_i +i;
    
    % draw the first threshold and the rest half of the time
    Screen('FrameOval', scr.window,1, stim.centeredthreshold_1,5);
    
    % if they are in the frame we draw them, then do it.
    if showCircles(i) == 1
        Screen('FrameOval', scr.window,1, stim.centeredthreshold_2,5);
        Screen('FrameOval', scr.window,[1 0 0], stim.missThreshold,5);
    end
    Screen('FrameOval', scr.window,[1 0.5 0 ], stim.missTarget,5);
    
    %The coin is flickering too
    if showCircles(i) == 1
        Screen('DrawTexture', scr.window, stim.imageTextures(stim.imageTextureIdx), [],stim.centeredMoney);
    end
    
    % Flip to the screen
    speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
    
end

% Stop showing threshold_2 and the red limit, they failed
stim.showCircles = false;
end

