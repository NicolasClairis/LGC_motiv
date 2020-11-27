function [doWin, forceSignal, onsets] = LGCM_physical_effort_stimulusPresentation(scr, stim, speed, audio_fbk_yn, sound)
%[doWin, forceSignal, onsets] = LGCM_physical_effort_stimulusPresentation(scr, stim, speed, audio_fbk_yn, sound)
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
% forceSignal: force signal during this whole period
%
% onsets: structure with timings

%% extract main variables
% screen
window = scr.window;
xCenter = scr.xCenter;
yCenter = scr.yCenter;
% stimulus
forceSignalColor = stim.signalColor;
threshold_1 = stim.threshold_1;
threshold_2 = stim.threshold_2;
% threshold to consider around the center
screenXcenter_inc_threshold = stim.inc_threshold;
% % timings
% t_effortPhase = t_wait.effortTotalPhase;

% reset some variables at each stimulus presentation. Initialize others
%% need to add description of the role of each of these variables
jFrame = 0; % keep track of frame number ~equivalent to time sample
frame_center_i = 1; % serves to know if the size of the money is ok
stim.atCenter = 0; % is the money at the center yet?
stim.showCircles = true; % if missed will switch to false
stim.triggerCatch = false; % did subject started to try to reach threshold2?
enlargedMoneyIdx = round(linspace(threshold_1(3), threshold_2(3)+50, 3*60)); % index of money enlargement speed
forceSignal = []; % initialize force signal
t2 = true;
firstT2 = 0;% time when threshold 2 reached

is_threshold_1_reached = 0;
is_threshold_2_reached = 0;

% load incentive image
inc_img = stim.inc;

%% In case of a problem, win is not a boolean but 2 so that we can discard the trial
doWin = 2;

%% start effort period
t_effortStart = GetSecs;

% Eternal loop, get out only for specific conditions
while 1
    
    
    %% read BioPac force output and prepare a circle corresponding to the signal
    jFrame = jFrame + 1; % update time index
    current_ForceLevel_tmp = fread(stim.u_out,1,'double');
    forceSignal(jFrame) = current_ForceLevel_tmp;
    
    %% display corresponding level of force being made in live with a green circle
    baseSignal = [0,...
        0,...
        round((500/70)*current_ForceLevel_tmp) + 1,...
        round((500/70)*current_ForceLevel_tmp) + 1]; % what are those 500/70 values?
    centeredSignal = CenterRectOnPointd(baseSignal, xCenter, yCenter);
    Screen('FillOval', window, forceSignalColor, centeredSignal);
    
    %% draw the 2 thresholds that have to be reached
    Screen('FrameOval', window, 1, stim.centeredthreshold_1, 5);
    Screen('FrameOval', window, 1, stim.centeredthreshold_2, 5);
    
    %% display the incentive level on the left until first threshold has been reached
    % then start making it move towards the center
    % then make it grow if successfull trial or leave to the right if
    % failure trial
    switch is_threshold_1_reached
        %% threshold 1 has not been reached yet
        % (=> means that infinite time could pass here if the subject does
        % not reach the threshold, we should add some maximal timing)
        case 0
            
            % draw the incentive money on the left of the screen
            Screen('DrawTexture', window, inc_img, [],leftMoneyRect);
            
            %     speed.vbl  = Screen('Flip', window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
            [~,timenow] = Screen('Flip',window);
            
            % when threshold 1 reached => extract the timing and update the
            % index
            if current_ForceLevel_tmp >= threshold_1_forceLevel
                % extract the timing
                onsets.threshold_1_reached = timenow;
                
                % update variable since threshold 1 has been reached
                is_threshold_1_reached = 1;
            end
        
        %% threshold 1 has been reached
        % => money starts moving from left to right
        case 1
            
            %% threshold 2 has not been reached yet
            switch is_threshold_2_reached
                case 0
                    % money moves from the left to the center of the
                    % screen
                    squareX = squareX + speed.pixelsPerFrame;
                    % display the new position of the money incentive
                    movingMoneyRect  = CenterRectOnPointd(stim.money,squareX, yCenter);
                    Screen('DrawTexture', window, inc_img, [],movingMoneyRect);
                    [~,timenow] = Screen('Flip',window);
                    
                    % the money is still on the left of the screen
                    if squareX <= xCenter - screenXcenter_inc_threshold
                        % extract the timing and the position of the money
                        onsets.(['money_left_X_coord_',num2str(squareX)]) = timenow;
                        
                        %% the money is at the center
                    elseif squareX > xCenter - screenXcenter_inc_threshold &&...
                            squareX <= xCenter + screenXcenter_inc_threshold
                        % extract the timing and the position of the money
                        onsets.(['money_center_X_coord_',num2str(squareX)]) = timenow;
                        
                        % if current force is higher than threshold 2
                        % => then the trial is a win and you go onto the
                        % last phase of the trial (money gets bigger while
                        % the force is higher than threshold 2 and you win
                        % the trial if you keep it like that) otherwise =
                        % miss
                        if current_ForceLevel_tmp > threshold_2_forceLevel
                            
                        end
                        
                        %% threshold 2 has not been reached yet and the money left the center
                        % => missed trial
                    elseif squareX > xCenter + screenXcenter_inc_threshold
                        % extract the timing and the position of the money
                        onsets.(['money_missed_right_X_coord_',num2str(squareX)]) = timenow;
                        
                        % store that the trial was missed
                        doWin = 0;
                    end
                    
                %% threshold 2 has been reached while money was at the center
                case 1
                %% the subject keeps performing a force higher than
                % threshold 2 => the money gets larger and larger until the
                % trial ends
                if current_ForceLevel_tmp > threshold_2_forceLevel - threshold_2_forceLevel_tolerance_wdw && doWin == 2
                    if moneySize < moneyEndSize
                        moneySize = moneySize + 1;
                        movingMoneyRect  = CenterRectOnPointd(stim.money,squareX, yCenter);
                        Screen('DrawTexture', window, inc_img, [],movingMoneyRect);
                        [~,timenow] = Screen('Flip',window);
                        
                    else % end of the trial when money reached maximum size
                        movingMoneyRect  = CenterRectOnPointd(stim.money,squareX, yCenter);
                        Screen('DrawTexture', window, inc_img, [],movingMoneyRect);
                        [~,timenow] = Screen('Flip',window);
                        onsets.trial_end = timenow;
                        % trial is a success
                        doWin = 1;
                    end
                    
                    %% the subject's performance dropped below the tolerance threshold
                    % => the trial is a loss
                elseif current_ForceLevel_tmp < threshold_2_forceLevel - threshold_2_forceLevel_tolerance_wdw
                    
                    % money moves from the center to the right of the
                    % screen
                    squareX = squareX + speed.pixelsPerFrame;
                    % display the new position of the money incentive
                    movingMoneyRect  = CenterRectOnPointd(stim.money,squareX, yCenter);
                    Screen('DrawTexture', window, inc_img, [],movingMoneyRect);
                    [~,timenow] = Screen('Flip',window);
                    
                    % trial is a failure
                    doWin = 0;
                end
            end % has threshold 1 has been reached during the trial?
            
    end % has threshold 1 been reached?
    
end % eternal loop

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