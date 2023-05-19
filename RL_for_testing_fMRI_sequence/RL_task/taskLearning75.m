% function [] = taskLearning75(subid,nrun)
%  MBB motivational battery (version INSEP).
%  taskLearning
%
%  Instrumental learning with monetary gain and loss. 
%  3 pairs (positive, negative, each 24 trials, and neutral, 12 trials).
%  Behavioural test
%  Mathias Pessiglione December 2009
%  Adapted for the motivational battery - Sept 2013. Raphael Le Bouc
%  Re-adapted for fMRI - February 2017 - Nicolas Clairis


clear all; close all; clc;

testing_script = 0; % parameters of screen to check on my own computer (1) or (0) when at CENIR

%% Identification
subid = input('subject identification number? ','s');
subject = str2double(subid);
nbTotalRuns = 4; % total number of runs (all tasks considered)
nrun = input('Run number ?');
runname = num2str(nrun);
IRM = 1;
nbTotalSessions = 4; % 8 runs 4 learning75
taskName = 'learning';

% identification_batmotiv: define name of the main results file = resultname
[resultname] = identification_batmotiv(taskName, subid, runname);

% setDir: create directory with results for all subjects (resultdir) if did
% not exist + folder for this subject (subdir) + folder where behavioral
% results stored for this subject (behaviordir)
% + define "root" folder
[root, resultdir, subdir, behaviordir, fMRIScansDir] = setDir(subid, IRM);

%% check if run already made or not (for the task only, not training)
[MadeOrNot] = checkIfRunMade(behaviordir, subid, nrun, nbTotalSessions);
if MadeOrNot == 1
    return; % stop task if run number already made
else
    clear('MadeOrNot')
end

%% Screen configuration
[x, y, window, baselineTextSize] = ScreenConfiguration(IRM, testing_script);

%% Loading images
% cross for ITI
picCross = Screen('MakeTexture',window,imread('Cross.bmp'));
rectCross = CenterRectOnPoint(Screen('Rect',picCross),x,y);
% feedback
cd('learningim') % go to folder with pictures of this task
gainFdbk = Screen('MakeTexture',window,imread('gain.bmp'));
lookFdbk = Screen('MakeTexture',window,imread('neutral.bmp'));
lossFdbk = Screen('MakeTexture',window,imread('loss.bmp'));
cd .. % back to main folder
rectGain = CenterRectOnPoint(Screen('Rect',gainFdbk),x,y);
rectLook = CenterRectOnPoint(Screen('Rect',lookFdbk),x,y);
rectLoss = CenterRectOnPoint(Screen('Rect',lossFdbk),x,y);

% assigning stimuli images to conditions (Reward/Neutral/Loss)
% alternate the meaning of symbols (rewarding/nothing, loss/nothing) in a
% pair across subjects: pair/impair subjects

% define order of pairs of stimuli across runs
pair_order = 1:nbTotalRuns;
pairForThisRun = pair_order(nrun);

% define attribution of pairs across subjects (not always the same pair for
% R/r, N/n and P/p)
pairPermut = perms(1:3);
pairPermut_idx = randi(size(pairPermut,1),1);
nstim = pairPermut(pairPermut_idx,:);

% define attribution of each item of a pair across subjects
if (subject/2) == floor(subject/2) % for pair subject A=A image, B=B image
    randomizeA = 'A';
    randomizeB = 'B';
else % for impair subjects A = images with B name, B= images with A name
    randomizeA = 'B';
    randomizeB = 'A';
end

% load corresponding images for each condition (Reward/Neutral/Loss)
cd([root, '\learningim\']) % go to folder with pictures of this task
[stimBest, stimWorse] = deal(cell(1,3));
for iCondition = 1:3
    stimBest{iCondition} = Screen('MakeTexture',window,imread(['Stim', num2str(pairForThisRun), num2str(nstim(iCondition)), randomizeA,'.bmp'])); % "good" stimulus (Reward/ no-loss)
    stimWorse{iCondition} = Screen('MakeTexture',window,imread(['Stim', num2str(pairForThisRun), num2str(nstim(iCondition)), randomizeB,'.bmp'])); % "bad" stimulus (no reward/loss)
end
cd .. % back to main folder
% define position of the stimuli on the screen
stimRectSize{1} = CenterRectOnPoint(Screen('Rect',stimBest{1}),   x/2, y); % Left position
stimRectSize{2} = CenterRectOnPoint(Screen('Rect',stimBest{1}), 3*x/2, y); % Right position

%% main parameters for pairs of stimuli
% number of trials per pair condition
ntrialsGain = 24;
ntrialsLoss = 24;
ntrialsNeutral = 12;
% total number of trials
totalTrial = ntrialsGain + ntrialsLoss + ntrialsNeutral;
% use independent probabilities for the two items (1) or dependent ones (0)
indpProb = 0;
% create trial vectors. This script defines the following variables:
% lottery (probability of winning/losing for each trial)
% npair (gain, neutral or loss pair trial)
% side (left or right side for the best option (reward/no-loss)
[side, npair, lottery] = designLearning75_bis(ntrialsGain, ntrialsLoss, ntrialsNeutral, indpProb);
% designLearning75;
%% data to save
trial = 1:totalTrial;
rt_fp = NaN(1,totalTrial); % RT first press
choice = NaN(1,totalTrial);
[response, feedback, gain] = deal(zeros(1,totalTrial));

%% Generator reset
rand('state',sum(100*clock));

%% Timings
instruction_T = 3;
cross_ITI_T = 0.5;
% ITT = 5; % inter-task time (pause between each task)

% final cross
ending_T = 10;
release_Wait = 0.5; % duration of cross after "please release"
PressTime_Wait = 3; % max time for waiting a response
dispChoice_Wait = 0.5; % 500ms display chosen option
% feedback jitter display between 1 to 3.95 seconds
jitters = 1:0.05:3.95;
varFdbkPerms = randperm(60);
varFdbk_Wait = jitters(varFdbkPerms);
% display of total gain at the end of the task
totalGain_Wait = 2;

%% events onsets and durations definition
list = {'cross_ITI',... % fixation cross for ITI
    'displayOptions',... % display options on screen
    'choice',... % moment when choice of an option is made
    'displayChoice',... % display chosen option
    'gainFeedback',... % feedback win/lose/neutral
    'lossFeedback',...
    'neutralFeedback',...
    'pleaseRelease',... % for cases where subject presses some button during the ITI
    'crossRelease',... % fixation cross after a "Please Release" message
    'tooSlowTrial_displayOptions',... % too slow trials
    'tooSlow',...
    'tooSlowTrial_feedback',...
    'pleaseKeepTrial_displayOptions',... % trials when subject doesn't maintain one key until the end of the trial
    'pleaseKeep',...
    'pleaseKeepTrial_feedback',...
    'onlyOneButtonTrial_displayOptions',...% trials when subject pressed two buttons at the end
    'onlyOneButton',...
    'onlyOneButtonTrial_feedback'};

for l = 1:length(list)
    onset.(list{l}) = NaN(totalTrial,1);
    duration.(list{l}) = NaN(totalTrial,1);
end
onset.TTL = []; % record TTL in that variable also to have all timings here
onset.instruction = [];
duration.instruction = [];
onset.cross = []; % fixation cross when no ITI (start of task, end of task)
duration.cross = [];
onset.totalGain = []; % display of total gain at the end of the run
duration.totalGain = [];

%% give/remind instructions before initial cross
[timenow1] = disp_instructions(window, x,y, IRM, taskName);
onset.instruction = timenow1;
if IRM == 0
    KbWait;
elseif IRM == 1
    WaitSecs(instruction_T);
end
timenow2 = GetSecs;
duration.instruction = timenow2 - timenow1;

%% keyboard keys configuration + waiting and recording first TTL for fMRI
if IRM == 0
    %% key configuration
    KbName('UnifyKeyNames');
    key.left = KbName('LeftArrow');
    key.right = KbName('RightArrow');
    key.escape= KbName('escape');
elseif IRM == 1
    %% fMRI key configuration
    KbName('UnifyKeyNames');
    key.left = 49;        % 49 % DROITE  bleu, left press
    key.right = 50;       %50   %% GAUCHE JAUNE
    key.escape = KbName('escape');
    
    % TTL trigger
    trigger = 54;
    
    %% add a cross at the beginning
    Screen('DrawTexture',window,picCross,[],rectCross);
    % onset
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.cross = [onset.cross; timenow1];
    
    %% TRIGGER RAPHAEL
    if IRM == 1 % useless for training
        disp('Now waiting for first TTL to start');
        dummy_scans = 1; % number of TTL to wait before starting the task (dummy scans are already integrated in CIBM scanner)
        [T0] = TTL_wait(dummy_scans, trigger);
        
        % record all subsequent TTL in the whole task
        keysOfInterest = zeros(1,256);
        keysOfInterest(trigger) = 1; % check TTL
        keysOfInterest(key.left) = 1; % check all left key press
        keysOfInterest(key.right) = 1; % check all right key press
        KbQueueCreate(0,keysOfInterest); % checks TTL and keys of pad
        KbQueueStart; % starts checking
    end
    
    % duration starting cross
    timenow2 = GetSecs;
    dur = timenow2 - timenow1;
    duration.cross = [duration.cross; dur];
end

%% behavioural task

% if testing, allow to escape from the task, but avoid this while real fMRI
% scan, because you don't want someone to press it by mistake and stop the
% task....
if testing_script == 1
    stoptask = 0;
end
for iTrial = 1:totalTrial

    %% fixation cross
    Screen('DrawTexture',window,picCross,[],rectCross);
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    onset.cross_ITI(iTrial) = timenow1;
    WaitSecs(cross_ITI_T);
    % duration
    timenow2 = GetSecs;
    duration.cross_ITI(iTrial) = timenow2 - timenow1;
    
    %% check if all keys are up before starting the trial
    keys_are_not_up = 0;
    [keyisdown, secs, keycode] = KbCheck;
    if keyisdown == 1 && (keycode(key.left) == 1 || keycode(key.right) == 1) % if subjects presses one of the two keys
        keys_are_not_up = 1;
        while keys_are_not_up == 1
            [keyisdown, secs, keycode] = KbCheck;
            DrawFormattedText(window,'Relachez les boutons svp','center','center', [255 255 255]);
            % onset
            [~,timenow1,~,~,~] = Screen(window,'Flip');
            onset.pleaseRelease(iTrial) = timenow1;
            if keycode(key.left) == 0 && keycode(key.right) == 0 % check if these keys are free
                keys_are_not_up = 0;
            end
        end
        % duration
        timenow2 = GetSecs;
        duration.pleaseRelease(iTrial) = timenow2 - timenow1;
        
        % display a cross again before starting the trial
        Screen('DrawTexture',window,picCross,[],rectCross);
        % onset
        [~,timenow1,~,~,~] = Screen(window,'Flip');
        onset.crossRelease(iTrial) = timenow1;
        WaitSecs(release_Wait);
        % duration
        timenow2 = GetSecs;
        duration.crossRelease(iTrial) = timenow2 - timenow1;
    end
    
    %% display stimuli on screen
    Screen('DrawTexture', window, stimBest{npair(iTrial)}, [], stimRectSize{( side(iTrial) + 3)/2});
    Screen('DrawTexture', window, stimWorse{npair(iTrial)}, [], stimRectSize{(-side(iTrial) + 3)/2});
    [~,timenow1,~,~,~] = Screen(window,'Flip');
    
    %% Check Response
    % should only take into account the last key pressed when the timing
    % ends.
    % If they didn't press any key => tooslow trial, corresponding feedback
    % If they pressed but didn't keep it until the end of the trial, then
    % display a corresponding feedback
    
    noKeyPressAtAll = 1; % check if at least one key was pressed during the trial
    keyMaintenance = 0; % check if key pressed was maintained pressed until the end of the trial
    bothKeysPressed = 0; % if two keys pressed at the end of the trial
    
    % check key press during 3 secs, focus on last key press only (for
    % the rest KbQueue will do the job!)
    trialTime = 0;
    while trialTime <= 3
        timenow3 = GetSecs;
        trialTime = timenow3 - timenow1;
        % Record Response
        [keyisdown, pressTime, keycode] = KbCheck;
        
        % check first time a key is pressed in this trial
        if noKeyPressAtAll == 1 && keyisdown == 1 && (keycode(key.left) == 1 || keycode(key.right) == 1)
            noKeyPressAtAll = 0; % at least one key was pressed during the trial, even if it was released before the end
            onset.choice(iTrial) = pressTime; % record time of first press in this trial
            rt_fp(iTrial) = pressTime - timenow1; % record RT between first press and start of trial
        end
    end
    
    % check last key pressed
    if noKeyPressAtAll == 1 % no key press at all = too slow/sleepy
        choice(iTrial) = 0; % no choice was made because no key pressed at all
        onset.tooSlowTrial_displayOptions(iTrial) = timenow1;
        duration.tooSlowTrial_displayOptions(iTrial) = pressTime - timenow1;
        
    elseif noKeyPressAtAll == 0 % at least one key was pressed during this trial
        if keyisdown == 0 || (keycode(key.left) == 0 && keycode(key.right) == 0) % pressed left or right key before the end of the trial, but didn't keep it until the end
            choice(iTrial) = 0; % no choice was made because didn't keep key pressed
            keyMaintenance = 0; % didn't maintain the pressure on the key
            onset.pleaseKeepTrial_displayOptions(iTrial) = timenow1;
            duration.pleaseKeepTrial_displayOptions(iTrial) = pressTime - timenow1;
            
        elseif keyisdown == 1 % one or two keys pressed until the end of the trial
            keyMaintenance = 1; % at least one key was maintained until the end of the trial
            if keycode(key.left) == 0 && keycode(key.right) == 1 % pressed right key
                onset.displayOptions(iTrial) = timenow1;
                duration.displayOptions(iTrial) = pressTime - timenow1;
                choice(iTrial) = 1; % 1 = right chosen
                
            elseif keycode(key.left) == 1 && keycode(key.right) == 0 % pressed left key
                onset.displayOptions(iTrial) = timenow1;
                duration.displayOptions(iTrial) = pressTime - timenow1;
                choice(iTrial) = -1; % -1 = left chosen
                
            elseif keycode(key.left) == 1 && keycode(key.right) == 1 % pressed both keys
                onset.onlyOneButtonTrial_displayOptions(iTrial) = timenow1;
                duration.onlyOneButtonTrial_displayOptions(iTrial) = pressTime - timenow1;
                choice(iTrial) = 0; % no choice was made because 2 buttons pressed
                bothKeysPressed = 1;
                
            elseif testing_script == 1 && keyisdown == 1 && keycode(key.escape) == 1 % stop the task (pressed escape)
                % if testing, allow to escape from the task, but avoid this while real fMRI
                % scan, because you don't want someone to press it by mistake and stop the
                % task....
                stoptask = 1;
                if stoptask % stop the task
                    break
                end
            end
        end
    end
    
    %% Show Response
    if noKeyPressAtAll == 0 % subject pressed at least one key during the trial
        if keyMaintenance == 1 % subject pressed one or two keys until the end of the trial
            if bothKeysPressed == 0 % subject pressed only one key until the end of the trial
                % case where subject answered normally to one of the two cues
                % (answered to only one cue and maintained pression until the
                % end of the trial)
                Screen('TextSize', window, 30);
                Screen('DrawTexture', window, stimBest{npair(iTrial)}, [], stimRectSize{(side(iTrial)+3)/2});
                Screen('DrawTexture', window, stimWorse{npair(iTrial)}, [], stimRectSize{(-side(iTrial)+3)/2});
                
                switch choice(iTrial)
                    case 1 % right option chosen => red square around it
                        L = 3*x/2 - 100;
                        T = y - 100;
                        R = 3*x/2 + 100;
                        B = y + 100;
                    case -1 % left option chosen => red square around it
                        L = x/2 - 100;
                        T = y - 100;
                        R = x/2 + 100;
                        B = y + 100;
                end
                
                Screen('DrawLine', window, [255 0 0], L, T, L, B, 3); % left
                Screen('DrawLine', window, [255 0 0], R, T, R, B, 3); % right
                Screen('DrawLine', window, [255 0 0], L, T, R, T, 3); % top
                Screen('DrawLine', window, [255 0 0], L, B, R, B, 3); % bottom
                [~,timenow1,~,~,~] = Screen(window, 'Flip');
                onset.displayChoice(iTrial) = timenow1;
                WaitSecs(dispChoice_Wait);
                % duration
                timenow2 = GetSecs;
                duration.displayChoice(iTrial) = timenow2 - timenow1;
                
                % if two buttons were pressed at the end of the trial
            elseif bothKeysPressed == 1
                DrawFormattedText(window,'Un seul bouton!','center','center', [255 255 255]);
                [~,timenow1,~,~,~] = Screen(window, 'Flip');
                onset.onlyOneButton(iTrial) = timenow1;
                WaitSecs(dispChoice_Wait);
                % duration
                timenow2 = GetSecs;
                duration.onlyOneButton(iTrial) = timenow2 - timenow1;
            end
            
            % if a key was pressed during the trial, but no key pressed at
            % the end
        elseif keyMaintenance == 0
            DrawFormattedText(window,'Maintenez votre appui!','center','center', [255 255 255]);
            [~,timenow1,~,~,~] = Screen(window, 'Flip');
            onset.pleaseKeep(iTrial) = timenow1;
            WaitSecs(dispChoice_Wait);
            % duration
            timenow2 = GetSecs;
            duration.pleaseKeep(iTrial) = timenow2 - timenow1;
        end
        
    % if no key at all was pressed during this trial
    elseif noKeyPressAtAll == 1
        DrawFormattedText(window,'Reveillez-vous!','center','center', [255 255 255]);
        [~,timenow1,~,~,~] = Screen(window, 'Flip');
        onset.tooSlow(iTrial) = timenow1;
        WaitSecs(dispChoice_Wait);
        % duration
        timenow2 = GetSecs;
        duration.tooSlow(iTrial) = timenow2 - timenow1;
    end
    
    %% Show feedback
    response(iTrial) = side(iTrial)*choice(iTrial); % -1 = incorrect, 1 = correct
    % 1 = side chosen and side of generally best option are the same
    % -1 = opposed sides for generally best option and choice made
    feedback(iTrial) = response(iTrial)*lottery(iTrial); % -1 = bad feedback, 1 = good feedback
    % 1 = generally worst option chosen but unlikely proba or generally
    % best option and likely proba (success)
    % -1= generally worst option and likely proba or generally best option and unlikely proba (fail)

    switch npair(iTrial)
        % win
        case 1
            if feedback(iTrial) == 1
                gain(iTrial) = 1;
                Screen('DrawTexture', window, gainFdbk, [], rectGain);
            else
                Screen('DrawTexture', window, lookFdbk, [], rectLook);
            end
            [~,timenow1,~,~,~] = Screen(window,'Flip');
            onset.gainFeedback(iTrial) = timenow1;
            
        % neutral
        case 2
            Screen('DrawTexture', window, lookFdbk, [], rectLook); % neutral feedback whatever option of the pair chosen (100%)
            [~,timenow1,~,~,~] = Screen(window,'Flip');
            onset.neutralFeedback(iTrial) = timenow1;
            
        % loss
        case 3
            if feedback(iTrial) == -1 || feedback(iTrial) == 0 % if answered loss element or if didn't answer => loss
                gain(iTrial) = -1;
                Screen('DrawTexture', window, lossFdbk, [], rectLoss);
            else
                Screen('DrawTexture', window, lookFdbk, [], rectLook);
            end
            [~,timenow1,~,~,~] = Screen(window,'Flip');
            onset.lossFeedback(iTrial) = timenow1;
    end
    WaitSecs(varFdbk_Wait(iTrial));
    % duration
    timenow2 = GetSecs;
    dur = timenow2 - timenow1;
    if npair(iTrial) == 1 && feedback(iTrial) == 1 % win
        duration.gainFeedback(iTrial) = dur;
    elseif npair(iTrial) == 3 && feedback(iTrial) == -1 % loss
        duration.lossFeedback(iTrial) = dur;
    else % neutral
        duration.neutralFeedback(iTrial) = dur;
    end
    
    % display number of the trial so that the experimenter can keep track
    % on how much there is left at the end of each trial
    disp(['Trial ',num2str(iTrial),'/',num2str(totalTrial)]);
    % save all stuff at the end of each trial in the subject's results directory
    save([behaviordir,filesep,'global_sub_',subid,'_session_',runname,'_',taskName,'.mat'])
end

%% fixation cross (or black screen) at the end of the script
Screen('DrawTexture',window,picCross,[],rectCross);
% onset
[~,timenow1,~,~,~] = Screen(window,'Flip');
onset.cross = [onset.cross; timenow1];
WaitSecs(ending_T);
% duration
timenow2 = GetSecs;
dur = timenow2 - timenow1;
duration.cross = [duration.cross; dur];

%% announce total gain after this run here
%% here or at the beginning you should load previous totalGain if this is not the first run
if nrun == 1 % for the first run, define the totalgain variable
    totalGain = nan(1,nbTotalRuns);
    totalGain(nrun) = sum(gain)*10; % multiply by 10 to match with actual values on screen
elseif nrun > 1 % otherwise load totalGain and add gain of this run to it
    load([behaviordir,'\totalGain_sub_',subid,'.mat'],'totalGain')
    if nrun == 4
        totalGain(nrun) = totalGain(1) + sum(gain)*10;
    elseif nrun == 7
        totalGain(nrun) = totalGain(4) + sum(gain)*10;
    end
end

save([behaviordir,'\totalGain_sub_',subid,'.mat'],'totalGain')
timenow1 = totalGainDisp(window, totalGain(nrun));
onset.totalGain = timenow1;
WaitSecs(totalGain_Wait);
% duration
timenow2 = GetSecs;
duration.totalGain = timenow2 - timenow1;

%% save data
cd(behaviordir);
session(1:totalTrial) = nrun;
data = [session; trial; npair; side; lottery; rt_fp; choice; response; feedback; gain].';
save(resultname,'subid','nrun','data','npair','side','lottery','rt_fp','choice','response','feedback','gain',...
    'pairForThisRun','nstim','randomizeA','randomizeB');

%% save all TTL and onsets in a results file
if IRM == 1
    KbQueueStop; % stop recording key presses
    keyLeft.Start = []; % time when starts pressing left key
    keyLeft.Release = []; % time when releases right key
    keyRight.Start = []; % time when starts pressing right key
    keyRight.Release = []; % time when releases right key
    while KbEventAvail
        [event, n] = KbEventGet;
        if event.Keycode == trigger
            TTL = [TTL; event.Time];
        elseif event.Keycode == key.left % if left key pressed
            if event.Pressed == 1 % record start of press
                keyLeft.Start = [keyLeft.Start; event.Time];
            elseif event.Pressed == 0 % record time when release
                keyLeft.Release = [keyLeft.Release; event.Time];
            end
        elseif event.Keycode == key.right % if right key pressed
            if event.Pressed == 1 % record start of press
                keyRight.Start = [keyRight.Start; event.Time];
            elseif event.Pressed == 0 % record time when release
                keyRight.Release = [keyRight.Release; event.Time];
            end
        end
    end
    KbQueueRelease;
    save([behaviordir,'\TTL_sub',subid,'_sess',runname,'.mat'],'T0');
    save([behaviordir,'\MBB_battery_',taskName,'_onsets_sub',subid,'_sess',runname,'.mat'],'onset','duration');
end

%% save all stuff
save([behaviordir,'\global_sub_',subid,'_session_',runname,'_',taskName,'.mat']) % save all variables here
cd(root); % go back to scripts folder
Screen('CloseAll'); % close PTB

%% check whether performance is above random for R/r and L/l trials, ignore neutral ones
% check R/r trials
bestChosenRr = (sum(response(npair == 1))/sum(npair==1))*100;
% check L/l trials
bestChosenLl = (sum(response(npair == 3))/sum(npair==3))*100;

% display all in one figure
figure
% R/r trials
if bestChosenRr > 50
    bar(1,bestChosenRr,'g'); % correctly understood
else
    bar(1,bestChosenRr,'r'); % slow learner or silly or too slow
end
% L/l trials
hold on
if bestChosenLl > 50
    bar(2,bestChosenLl,'g'); % correctly understood
else
    bar(2,bestChosenLl,'r'); % slow learner or silly or too slow
end
xlim([0 2.5]);
ylim([0 110]);
xticks(1:2);
xticklabels({'Gains','Losses'});
ylabel('% correct choices');

% display temporal performance also
figure
% R/r trials
RrPlot = subplot(2,1,1);
plot(response(npair == 1),'b');
ylim([-1.2 1.2]);
xlabel('Gain pair');
ylabel('% correct choices');
% L/l trials
LlPlot = subplot(2,1,2);
ylim([-1.2 1.2]);
plot(response(npair == 3),'b');
ylim([-1.2 1.2]);
xlabel('Loss pair');
ylabel('% correct choices');