function [onsets] = LGCM_training(scr, stim, speed,...
    audio_fbk_yn, sound,...
    session_effort_type, n_training_trials)
%[all, stim] = LGCM_training(scr, stim, speed,...
%     audio_fbk_yn, sound,...
%     session_effort_type, n_training_trials)
% LGCM_training will allow to perform the training trials before starting
% the main task
%
% INPUTS
%
% OUTPUTS
%
% See also LGCM_main_experiment.m

%% extract variables
% screen parameters
window = scr.window;
screenYpixels = scr.Ycenter*2;
bigTextSize = scr.textSize.big;
middleTextSize = scr.textSize.middle;
% timings
t_training_instruction = t_wait.training_instructions;

%% Show the training text
Screen('TextSize',window,middleTextSize)
DrawFormattedText(window, 'TRAINING', 'center',screenYpixels * 0.7, 1);
Screen('TextSize', window, bigTextSize)
DrawFormattedText(window, 'Get Ready', 'center', screenYpixels * 0.3, 1);

% speed.vbl  = Screen('Flip', window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
[~,timenow] = Screen('Flip',window);
onsets.training_instruction = timenow;
pause(t_training_instruction); % 3 seconds

%% Train participant on the task, so they know their limit before hand.
for iTrainingTrial = 1:n_training_trials
    [doWin,signal,firstT2] = stimulusPresentation(scr,stim,speed,sound);
end % training trial loop


end % function