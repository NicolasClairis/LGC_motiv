function[onsets] = LGCM_effort_phase(scr, speed, stim,...
    wait, audio_fbk_yn, n_trials, reward_or_punishment)
% will perform the main LGC motivational task

%% extract main variables for display
screenYpixels = scr.Ycenter*2;
window = scr.window;
t_block_instructions = wait.block_instructions;
t_cross = wait.fixation_cross;
t_wait_stimDisp = wait.stimulus_display;

%% Announce the next block
Screen('TextSize',window,100)
switch reward_or_punishment
    case 'R'
        DrawFormattedText(window, 'Next up : MISS or WIN Money ', 'center', screenYpixels*0.7, 1);
    case 'P'
        DrawFormattedText(window, 'Next up : SAVE or LOSE Money ', 'center', screenYpixels*0.7, 1);
end
Screen('TextSize', window, 150)
DrawFormattedText(window, 'Get Ready', 'center', scr.screenYpixels * 0.3, 1);

speed.vbl  = Screen('Flip', window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
warning('extract onsets here');
WaitSecs(t_block_instructions)


%% launch the task
for iTrial = 1:n_trials
    %% fixation cross display
    Screen('FillRect',window,white,stim.cross_verticalLine); % vertical line
    Screen('FillRect',window,white,stim.cross_horizontalLine); % horizontal line
    [~,timenow] = Screen('Flip',window); % display the cross on screen
    onsets.cross(iTrial) = timenow;
    WaitSecs(t_cross);
    
    %% display the incentive before effort starts
    centeredSignal = CenterRectOnPointd(baseSignal, xCenter, yCenter);
    Screen('FillOval', window, signalColor, centeredSignal);
    Screen('DrawTexture', window, image, [], centeredMoney);
    % speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
    [~, timenow] = Screen('Flip',window);
    onsets.incentive(iTrial) = timenow;
    WaitSecs(t_wait_stimDisp);
    
    %% effort phase
    [doWin, signal, onsets_tmp] = stimulusPresentation(scr,stim,iTrial,speed,sound);
    onsets.success(iTrial) = onsets_tmp.success;
    onsets.miss(iTrial) = onsets_tmp.miss;
    
    %% display trial feedback
    onsets.fbk(iTrial) = LGCM_trial_feedback(scr,stim,speed,totalMoney,doWin,stim.incentive(stim.incentiveIdx(iTrial)));
    
    %% display how many trials have been done (for the experimenter)
    disp(['Trial ',num2str(iTrial),'/',num2str(n_trials),' done']);
end % trial loop



end % function