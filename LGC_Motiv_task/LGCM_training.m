function [all, stim] = LGCM_training(scr, stim, speed,...
    audio_fbk_yn, sound,...
    session_effort_type)
%[all, stim] = LGCM_training(scr, stim, speed,...
%     audio_fbk_yn, sound,...
%     session_effort_type)
% LGCM_training will allow to perform the training trials before starting
% the main task


%% Show the training text
Screen('TextSize',scr.window,100)
DrawFormattedText(scr.window, 'TRAINING', 'center',scr.screenYpixels * 0.7, 1);
Screen('TextSize', scr.window, 150)
DrawFormattedText(scr.window, 'Get Ready', 'center', scr.screenYpixels * 0.3, 1);

%% Flip to the screen
speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
warning('add onsets here');
pause(t_); % 3 seconds

%% Train participant on the task, so they know their limit before hand.
n_trials = 3;
for iTrainingTrial = 1:n_training_trials
    stim.imageTextureIdx = stim.incentiveIdx(iTrainingTrial);
    [doWin,signal,firstT2] = stimulusPresentation(scr,stim,speed,sound);
    all.signals{1,iTrainingTrial} = signal;
    all.win{1,iTrainingTrial}  = doWin;
    all.VCmax{1,iTrainingTrial} = max(signal);
    all.trialLength{1,iTrainingTrial} = length(signal);
    all.incentive{1,iTrainingTrial} = stim.incentive(stim.incentiveIdx(iTrainingTrial));
    all.firstT2{1,iTrainingTrial} = firstT2;
end


end % function