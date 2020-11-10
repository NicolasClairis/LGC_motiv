function [all,stim] = training(scr,stim,speed,sound)

% Quick text to introduce first task
Screen('TextSize', scr.window, 100);
DrawFormattedText(scr.window, 'First : Maximum Power', 'center',scr.screenYpixels * 0.7, 1);
Screen('TextSize', scr.window, 150)
DrawFormattedText(scr.window, 'Get Ready', 'center', scr.screenYpixels * 0.3, 1);

% Flip to the screen
speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
pause(5)

% Measure MVC as normal
[measuredMVC,signalMVC] = measureMVC(scr,stim,speed);
all.measuredMVC{1,1} = measuredMVC;
all.MVCsignal{1,1} = signalMVC(1,:);
all.MVCsignal{1,2} = signalMVC(2,:);
all.MVCsignal{1,3} = signalMVC(3,:);
stim.measuredMVC = max(measuredMVC);
all.usedMVC{1,1} = stim.measuredMVC;

% Show the training text
Screen('TextSize',scr.window,100)
DrawFormattedText(scr.window, 'TRAINING', 'center',scr.screenYpixels * 0.7, 1);
Screen('TextSize', scr.window, 150)
DrawFormattedText(scr.window, 'Get Ready', 'center', scr.screenYpixels * 0.3, 1);

% Flip to the screen
speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
pause(3)

% Train participant on the task, so they know their limit before hand.
n_trials = 3;
for training_i = 1:n_trials
    stim.imageTextureIdx = stim.incentiveIdx(training_i);
    [doWin,signal,firstT2] = stimulusPresentation(scr,stim,speed,sound);
    all.signals{1,training_i} = signal;
    all.win{1,training_i}  = doWin;
    all.VCmax{1,training_i} = max(signal);
    all.trialLength{1,training_i} = length(signal);
    all.incentive{1,training_i} = stim.incentive(stim.incentiveIdx(training_i));
    all.firstT2{1,training_i} = firstT2;
end

% Announce the next block
Screen('TextSize',scr.window,100)
switch stim.isPositiveStimuli
    case 1
        DrawFormattedText(scr.window, 'Next up : MISS or WIN Money ', 'center', scr.screenYpixels * 0.7, 1);
    case 0
        DrawFormattedText(scr.window, 'Next up : SAVE or LOSE Money ', 'center', scr.screenYpixels * 0.7, 1);
end
Screen('TextSize', scr.window, 150)
DrawFormattedText(scr.window, 'Get Ready', 'center', scr.screenYpixels * 0.3, 1);

% Flip to the screen
speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
pause(3)


end

