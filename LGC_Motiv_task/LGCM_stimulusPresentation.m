function [doWin,signal,firstT2] = LGCM_physical_effort_stimulusPresentation(scr, stim, speed, audio_fbk_yn, sound)
%[doWin,signal,firstT2] = LGCM_physical_effort_stimulusPresentation(scr, stim, speed, audio_fbk_yn, sound)
% LGCM_physical_effort_stimulusPresentation will show the physical effort
% task
%
% INPUTS
% audio_fbk_yn
% 'yes': include audio feedback
% 'no': no audio feedback
%
% sound: structure with sound to use for win/lose
%
% OUTPUTS
%

%% extract main variables
% screen
window = scr.window;
xCenter = scr.xCenter;
yCenter = scr.yCenter;
% stimulus
signalColor = stim.signalColor;
threshold_1 = stim.threshold_1;
threshold_2 = stim.threshold_2;
% timings
t_effortPhase = t_wait.effortTotalPhase;

% reset some variables at each stimulus presentation. Initialize others
%% need to add description of the role of each of these variables
frame_i = 1;
frame_center_i = 1;
VCsignal = 0;
atCenter = 0;
showCircles = true;
triggerCatch = false;
centeredMoney = ;
enlargedMoneyIdx = round(linspace(threshold_1(3), threshold_2(3) + 50, 3*60));
signal = [];
t2 = true;
firstT2 = 0;

%% In case of a problem, win is not a boolean but 2 so that we can discard the trial
doWin = 2;

 
%% display incentive before effort starts
centeredSignal = CenterRectOnPointd(baseSignal, xCenter, yCenter);
Screen('FillOval', window, signalColor, centeredSignal);
Screen('DrawTexture', window, image, [], centeredMoney);
% speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
[~, timenow] = Screen('Flip',window);
onsets.incentive = timenow;
WaitSecs(t_wait_stimDisp);


%% start effort period
[~, t_effortStart] = Screen('Flip',window);
while timenow < t_effortStart + t_effortPhase
    % 
    
    % show 2 circle thresholds and miss area
    Screen('FrameOval', window,1, centeredthreshold_1, 5);
    Screen('FrameOval', window,1, centeredthreshold_2, 5);
    Screen('FrameOval', window,[1 0 0], missThreshold, 5);
    
    
    % extract time to see when to get out of the loop
    timenow = GetSecs;
end % time loop

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
            PsychPortAudio('Start', sound.pahandle, 1, 0, 1);
            
        case 0 % lost the coin => negative feedback
            PsychPortAudio('FillBuffer', sound.pahandle, [sound.audio_lose(:,1)'; sound.audio_lose(:,2)']);
            
            % Start audio playback no repetition (1), start instantly (0), wait for device to really start (1)
            PsychPortAudio('Start', sound.pahandle, 1, 0, 1);
    end
end % audio YES/NO

end % function