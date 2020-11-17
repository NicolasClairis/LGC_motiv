function[onsets] = LGCM_effort_phase(scr, speed, stim, wait, audio, n_trials, reward_or_punishment)
% will perform the main LGC motivational task

%% extract main variables for display
screenYpixels = scr.Ycenter*2;
window = scr.window;
t_block_instructions = wait.block_instructions;
t_cross = wait.fixation_cross;

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
    warning('extract onsets here');
    WaitSecs(t_cross);
    
    %% display the incentive
    
    %% effort phase
    stim.imageTextureIdx = stim.incentiveIdx(iTrial);
    [doWin, signal, firstT2] = stimulusPresentation(scr,stim,speed,sound);
    
    %% display trial feedback
    totalMoney = showResult(scr,stim,speed,totalMoney,doWin,stim.incentive(stim.incentiveIdx(iTrial)));
    
    %% display how many trials have been done (for the experimenter)
    disp(['Trial ',num2str(iTrial),'/',num2str(n_trials),' done']);
end % trial loop



end % function