function [doWin,forceSignal,onsets] = LGCM_physical_effort_stimulusPresentation(scr, stim, speed, audio_fbk_yn, sound)
%[doWin,signal,onsets] = LGCM_physical_effort_stimulusPresentation(scr, stim, speed, audio_fbk_yn, sound)
% LGCM_physical_effort_stimulusPresentation will show the physical effort
% task
%
% INPUTS
% scr: structure with main screen parameters
%
% stim
%
% speed:
%
% audio_fbk_yn:
% 'yes': include audio feedback
% 'no': no audio feedback
%
% sound: structure with sound to use for win/lose
%
% OUTPUTS
% doWin: did the subject win the trial (1) or did he miss it (0)?
%
% signal: force signal
%
% onsets: structure with timings

%% extract main variables
% screen
window = scr.window;
xCenter = scr.xCenter;
yCenter = scr.yCenter;
% stimulus
signalColor = stim.signalColor;
threshold_1 = stim.threshold_1;
threshold_2 = stim.threshold_2;
% % timings
% t_effortPhase = t_wait.effortTotalPhase;

% reset some variables at each stimulus presentation. Initialize others
%% need to add description of the role of each of these variables
jFrame = 1; % keep track of frame number ~equivalent to time sample
frame_center_i =1; % serves to know if the size of the money is ok
stim.atCenter = 0; % is the money at the center yet?
stim.showCircles = true; % if missed will switch to false
stim.triggerCatch = false; % did subject started to try to reach threshold2?
enlargedMoneyIdx = round(linspace(threshold_1(3), threshold_2(3)+50, 3*60)); % index of money enlargement speed
forceSignal = []; % force signal
t2 = true;
firstT2 = 0;% time when threshold 2 reached

%% In case of a problem, win is not a boolean but 2 so that we can discard the trial
doWin = 2;

%% start effort period
t_effortStart = GetSecs;

% Eternal loop, get out only for specific conditions
while 1
    
    
    % Make Read the data and prepare a circle corresponding to the signal
    forceSignal(jFrame) = fread(stim.u_out,1,'double');
    stim.baseSignal = [0 0 round((500/70)*forceSignal)+1 round((500/70)*forceSignal)+1];
    stim.centeredMoney  = CenterRectOnPointd(stim.money,0+scr.squareX, scr.yCenter);
    
    % Check different case condition, in which thing might happen
    % If the stimulus is in the center of the screen, but is underneath the minimum amount required then
    if stim.atCenter == 1 && forceSignal < 60 && doWin ~= false
        %% center + force lower than 60% => participant lost
        
        % Lost the coin and show animation of failure to aquire the coin
        doWin = 0;
        [stim,speed,jFrame,forceSignal] = missprotocol(scr,stim,speed,forceSignal,jFrame);
        
        % If the stimulus is at the center (+-10 pixels) and that the participant either triggered  catching the coin for the first time, or is trying to catch it.
    elseif ((forceSignal) > 70 && (scr.xCenter-10 < scr.squareX) && (scr.squareX < scr.xCenter+10) && (doWin ~= false)) || ((forceSignal) > 60 && (stim.triggerCatch == true) && (scr.xCenter-10 < scr.squareX) && (scr.squareX < scr.xCenter+10) && (doWin ~= false))
        %% if t2 reached (=70%) AND CENTER
        % => while force still higher than 60%
        % => piece money starts enlarging
        %
        % if one of these conditions is not reached => missed
        
        % the coin arrived at the center
        stim.atCenter =1;
        
        % Record first time they triggered the second threshold
        if t2 == true
            firstT2 = jFrame/60;
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
        %% end of the trial in case of miss:
        % money moves towards the right of the screen
        
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
    elseif forceSignal> 30
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
    jFrame = jFrame+1;
    
    % Keep track of time if wanted.
    aux(jFrame)=toc;
    
end


%% update total amounf of money depending on block type (reward/punishment)
% and on if the trial was successfull (doWin = 1) or not (doWin = 0)
switch reward_or_punishment
    case 'R'
        if doWin == 1
            totalMoney = totalMoney + amountMoney;
        end
    case 'P'
        if doWin == 0
            totalMoney = totalMoney - amountMoney;
        end
end

%% sound feedback
if strcmp(audio_fbk_yn,'yes')
    switch doWin
        case 1 % won the coin => positive feedback
            % Fill the audio playback buffer with the audio data, doubled for stereo presentation
            PsychPortAudio('FillBuffer', sound.pahandle, [sound.audio_win'; sound.audio_win']);
            % Start audio playback no repetition (1), start instantly (0), wait for device to really start (1)
            onsets.audio_fbk = PsychPortAudio('Start', sound.pahandle, 1, 0, 1);
            
        case 0 % lost the coin => negative feedback
            PsychPortAudio('FillBuffer', sound.pahandle, [sound.audio_lose(:,1)'; sound.audio_lose(:,2)']);
            
            % Start audio playback no repetition (1), start instantly (0), wait for device to really start (1)
            onsets.audio_fbk = PsychPortAudio('Start', sound.pahandle, 1, 0, 1);
    end
end % audio YES/NO

end % function