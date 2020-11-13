function [doWin,signal,firstT2] = stimulusPresentation(scr,stim,speed,sound)

% reset some variables at each stimulus presentation. Initialize others
frame_i =1;
frame_center_i =1;
stim.VCsignal = 0;
stim.atCenter = 0;
stim.showCircles = true;
stim.triggerCatch = false;
enlargedMoneyIdx = round(linspace(stim.threshold_1(3),stim.threshold_2(3)+50,3*60));
signal = [];
t2 = true;
firstT2 = 0;

% In case of a problem, win is not a boolean but 2 so that we can discard the trial
doWin = 2;

% Eternal loop, get out only for specific conditions
while 1
    
    % measure time
    
    tic
    
    % Make Read the data and prepare a circle corresponding to the signal
    stim.VCsignal = fread(stim.u_out,1,'double');
    stim.baseSignal = [0 0 round((500/70)*stim.VCsignal)+1 round((500/70)*stim.VCsignal)+1];
    stim.centeredMoney  = CenterRectOnPointd(stim.money,0+scr.squareX, scr.yCenter);
    
    %keep the signal to save it
    signal(frame_i) = stim.VCsignal;
    
    % Check different case condition, in which thing might happen
    % If the stimulus is in the center of the screen, but is underneath the minimum amount required then
    if stim.atCenter == 1 && stim.VCsignal < 60 && doWin ~= false
        
        % Lost the coin and show animation of failure to aquire the coin
        doWin = 0;
        [stim,speed,frame_i,signal] = missprotocol(scr,stim,speed,signal,frame_i);
        
        % If the stimulus is at the center (+-10 pixels) and that the participant either triggered  catching the coin for the first time, or is trying to catch it.
    elseif ((stim.VCsignal) > 70 && (scr.xCenter-10 < scr.squareX) && (scr.squareX < scr.xCenter+10) && (doWin ~= false)) || ((stim.VCsignal) > 60 && (stim.triggerCatch == true) && (scr.xCenter-10 < scr.squareX) && (scr.squareX < scr.xCenter+10) && (doWin ~= false))
        
        % the coin arrived at the center
        stim.atCenter =1;
        
        % Record first time they triggered the second threshold
        if t2 == true
            firstT2 = frame_i/60;
            t2 = false;
        end
        
        % the coin is trying to be caught
        stim.triggerCatch = true;
        
        % fix position of the coin to it's center (overide other position)
        scr.squareX = scr.xCenter;
        
        %Increase the coin size over time
        enlargedMoney = [0 0 enlargedMoneyIdx(frame_center_i) enlargedMoneyIdx(frame_center_i)];
        stim.centeredMoney  = CenterRectOnPointd(enlargedMoney,scr.xCenter, scr.yCenter);
        frame_center_i = frame_center_i +1;
        
        % If the coin had it's max size, then trigger victory
        if frame_center_i == size(enlargedMoneyIdx,2) && stim.isPositiveStimuli ~= 2
            % Fill the audio playback buffer with the audio data, doubled for stereo presentation
            PsychPortAudio('FillBuffer', sound.pahandle, [sound.audio_win'; sound.audio_win']);
            % Start audio playback no repetition (1), start instantly (0), wait for device to really start (1)
            PsychPortAudio('Start', sound.pahandle, 1, 0, 1);
            
            % Win and get out
            doWin = 1;
            break
        end
        
        % If we have the neutral stimuli, win but without sound
        if frame_center_i == size(enlargedMoneyIdx,2) && stim.isPositiveStimuli == 2
            doWin = 1;
            break
        end
        % If coin reaches the other end of the screen inside the miss threshold
    elseif (scr.windowRect(3)-stim.moneySize/2 -10 < scr.squareX) && (scr.squareX < scr.windowRect(3)-stim.moneySize/2+10)
        
        % loss the coin
        doWin = 0;
        
        % If we were in the punishment block, play lose sound
        if stim.isPositiveStimuli == false
            
            PsychPortAudio('FillBuffer', sound.pahandle, [sound.audio_lose(:,1)'; sound.audio_lose(:,2)']);
            
            % Start audio playback no repetition (1), start instantly (0), wait for device to really start (1)
            PsychPortAudio('Start', sound.pahandle, 1, 0, 1);
        end
        
        % end the trial as they lost
        break
        
        % If we are above the first threshold, move the coin
    elseif stim.VCsignal> 30
        scr.squareX = scr.squareX + speed.pixelsPerFrame;
    end
    
    
    % We set bounds to make sure our square doesn't go completely off of the screen
    if scr.squareX < 0
        scr.squareX = 0;
    elseif scr.squareX > scr.screenXpixels
        scr.squareX = scr.screenXpixels;
    elseif scr.squareY < 0
        scr.squareY = 0;
    elseif scr.squareY > scr.screenYpixels
        scr.squareY = scr.screenYpixels;
    end
    
    % Center the signal
    stim.centeredSignal = CenterRectOnPointd(stim.baseSignal,scr.xCenter, scr.yCenter);
    
    % Draw the circles to the screen
    Screen('FillOval', scr.window, stim.signalColor, stim.centeredSignal);
    
    % If it is the first frame of stimulus presentation, wait 1.5 second so that participant can see the value of the coin
    if stim.first_appearence == false
        
        % Only draw on the screen the coin they are playing for
        stim.first_appearence = true;
        Screen('DrawTexture', scr.window, stim.imageTextures(stim.imageTextureIdx), [],stim.centeredMoney);
        speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
        pause(1.5)
    end
    
    % If we did not win or lose, draw all our threshold and stimulus
    Screen('FrameOval', scr.window,1, stim.centeredthreshold_1,5);
    
    % If they lost trying to catch the coin, don't show those threshold anymore
    if stim.showCircles == true
        Screen('FrameOval', scr.window,1, stim.centeredthreshold_2,5);
        Screen('FrameOval', scr.window,[1 0 0], stim.missThreshold,5);
    end
    
    % Show stimulus and threshold. The order in which we draw to the screen is important. last element is the one in the front
    Screen('FrameOval', scr.window,[1 0.5 0 ], stim.missTarget,5);
    Screen('DrawTexture', scr.window, stim.imageTextures(stim.imageTextureIdx), [],stim.centeredMoney);
    
    % Flip to the screen
    speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
    
    % A frame passed, increment it
    frame_i = frame_i+1;
    
    % Keep track of time if wanted.
    aux(frame_i)=toc;
    
end

% Show a grey screen in between
speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);

end

