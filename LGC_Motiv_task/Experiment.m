% Clear the workspace and the screen, instrreset resets the udp channels
sca;
close all;
clearvars;
instrreset

%% working directories
savePath = ['D:' filesep 'Matlab_codes' filesep 'LGC_Motiv_task' filesep];

%% Open a connexion to the biopack datastream
  stim.u_out = udp('127.0.0.1', 2012, 'LocalPort', 15010);
%  stim.u_out = udp('127.0.0.1', 2012);
%      stim.u_out = udp('128.178.188.240', 16212, 'LocalPort', 15010);
     
% stim.u_out = udp('128.178.188.240', 16212);
% stim.u_out = udp('','LocalHost', '','LocalPort', 16213);
 fopen(stim.u_out);
% for i = 1:1000
%         fread(stim.u_out,1,'double')
%         i
% end
% stim.u_out = tcpip('128.178.188.240', 16212, 'NetworkRole','server');
% fopen(stim.u_out);
% for i = 1:1000
%         fread(stim.u_out,1,'double')
%         i
% end

% t = tcpip('localhost', 15010, 'NetworkRole', 'server');
% % set(t, 'InputBufferSize', 16); 
% fopen(t);
% data1 = fread(t, 16,'float');

%% subject identification
% Insert the initials, the number of the participants
info = inputdlg({'Initials', 'Subject ID'});
[init, subjectID] = info{[1,2]};

% Create subjectCodeName which is used as a file saving name
subjectCodeName = strcat(init,'_s',subjectID);

% Escape the experiment in case of no answer to either question
if isempty(init) == true || isempty(subjectID) == true
    errordlg('You didn''t answer everything ! We need all information to continue.')
end

%% Initialize variables

% PTB initialization + anti bug due to sync
% Here we call some default settings for setting up Psychtoolbox. So it doesn't crash when we call it
PsychDefaultSetup(2);
Screen('Preference','VisualDebugLevel', 0);
Screen('Preference', 'SkipSyncTests', 1);




% PTB sound
% Initialize Sounddriver and prepare 2 audio input we will use. Resample the second one
InitializePsychSound(1);
[sound.audio_win,sound.Fs_win] = audioread(['D:' filesep 'Matlab codes' filesep 'Experiment Motivation Original' filesep 'Sounds' filesep 'Win.mp3']);
[sound.audio_lose,sound.Fs_lose] = audioread(['D:' filesep 'Matlab codes' filesep 'Experiment Motivation Original' filesep 'Sounds' filesep filesep 'loose.mp3']);
[p,q] = rat(sound.Fs_win / sound.Fs_lose);
sound.audio_lose = resample(sound.audio_lose,p,q);
sound.numberChannels = 2;

% Open Psych-Audio port, with the follow arguements
% (1) [] = default sound device
% (2) 1 = sound playback only
% (3) 1 = default level of latency
% (4) Requested frequency in samples per second
% (5) 2 = stereo putput
sound.pahandle = PsychPortAudio('Open', [], 1, 1, sound.Fs_win, sound.numberChannels);
% sound.pahandle_lose = PsychPortAudio('Open', [], 1, 1, sound.Fs_lose, sound.nrchannels);

% Set the volume to half
PsychPortAudio('Volume', sound.pahandle, 0.5);




% Initialize variable only important at this layer
totalMoney = 20;
nbTrialPerCoinType = 1;
nbTrialPerPhase = [3 nbTrialPerCoinType*3 nbTrialPerCoinType*3 nbTrialPerCoinType];





% Screen variables
% GET classical screen information : screen numbers,
scr.screens = Screen('Screens');

% Draw to the external screen if avaliable
scr.screenNumber = max(scr.screens);
stim.VSsignal = 0;

% Define black and white
scr.white = WhiteIndex(scr.screenNumber);
scr.black = BlackIndex(scr.screenNumber);

scr.grey = 0.5;
% Open an on screen window
[scr.window, scr.windowRect] = PsychImaging('OpenWindow', scr.screenNumber, scr.grey);

% Get the size of the on screen window
[scr.screenXpixels, scr.screenYpixels] = Screen('WindowSize', scr.window);

% Get the centre coordinate of the window
[scr.xCenter, scr.yCenter] = RectCenter(scr.windowRect);

% Enable alpha blending for anti-aliasing of all our textures
Screen('BlendFunction', scr.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');





% Speed related variables
% Query the frame duration
speed.ifi = Screen('GetFlipInterval', scr.window);

% Prepare Rest interval in seconds
speed.isShortBreak = true;
speed.shortBreak = 10;
speed.longBreak = 60;

% Sync us and get a time stamp
speed.vbl = Screen('Flip', scr.window);
speed.waitframes = 1;
% Set the amount we want our square to move on each button press
speed.pixelsPerFrame = 6;




% Stimulus related variables
% The color used to represent the signal
stim.signalColor = [0 0.75 0.3];

% Money variables
stim.moneySize = 214;

stim.money = [0 0 stim.moneySize stim.moneySize];
% 30% of MVC defined to this size so that 70 can be at 500 pixels
stim.threshold_1 = [0 0 stim.moneySize stim.moneySize];
stim.missTarget = [0 0 stim.moneySize stim.moneySize];

% 60% of MVC
stim.missThreshold = [0 0 429 429];

% 70% of MVC
stim.threshold_2 = [0 0 500 500];

% First appearance of the stimulus to wait and let the participant see the value of the coin
stim.first_appearence = false;


% Is it going to be reward or punishment first
isPositiveStimuli = [randperm(numel([1 0]))-1 2];
stim.isPositiveStimuli = isPositiveStimuli(1);

% Set the intial position of the coin for screen. (even if not in the good place)
scr.squareX = stim.moneySize/2;
scr.squareY = scr.yCenter;

% Position of each ring/coin/signal on the screen
stim.centeredthreshold_1 = CenterRectOnPointd(stim.threshold_1,scr.xCenter, scr.yCenter);
stim.centeredthreshold_2 = CenterRectOnPointd(stim.threshold_2,scr.xCenter, scr.yCenter);
stim.missThreshold = CenterRectOnPointd(stim.missThreshold,scr.xCenter, scr.yCenter);
stim.missTarget = CenterRectOnPointd(stim.missTarget,scr.windowRect(3) - scr.squareX, scr.yCenter);

%Import the images
imageLocation_20Cent = ['D:' filesep 'Matlab codes' filesep 'Experiment Motivation Original' filesep 'Coin PNG' filesep 'SMT_Coin_20RP.png'];
[image_20Cent,~,alpha_20Cent] = imread(imageLocation_20Cent);
imageLocation_50Cent = ['D:' filesep 'Matlab codes'  filesep 'Experiment Motivation Original' filesep 'Coin PNG'  filesep 'SMT_Coin_50RP.png'];
[image_50Cent,~,alpha_50Cent] = imread(imageLocation_50Cent);
imageLocation_1FR = ['D:' filesep 'Matlab codes' filesep 'Experiment Motivation Original' filesep 'Coin PNG'  filesep 'SMT_Coin_1FR.png'];
[image_1FR,~,alpha_1FR] = imread(imageLocation_1FR);
imageLocation_0FR = ['D:' filesep 'Matlab codes' filesep 'Experiment Motivation Original' filesep 'Coin PNG'  filesep 'SMT_Coin_0FR.png'];
[image_0FR,~,alpha_0FR] = imread(imageLocation_0FR);

% Make the image into a texture
image_20Cent(:,:,4) = alpha_20Cent;
stim.imageTexture_20Cent = Screen('MakeTexture', scr.window, image_20Cent);
image_50Cent(:,:,4) = alpha_50Cent;
stim.imageTexture_50Cent = Screen('MakeTexture', scr.window, image_50Cent);
image_1FR(:,:,4) = alpha_1FR;
stim.imageTexture_1FR = Screen('MakeTexture', scr.window, image_1FR);
image_0FR(:,:,4) = alpha_0FR;
stim.imageTexture_0FR = Screen('MakeTexture', scr.window, image_0FR);

% Prepare an array of incentive, each type of coin being shown the same amount of time in a randomized order
incentiveIdx = cat(1,repelem([1 2 3],nbTrialPerCoinType*3));
stim.incentiveIdx = incentiveIdx(randperm(nbTrialPerCoinType*3));

% Prepare incentive texture
stim.incentive = [1,0.5,0.2,0];
stim.imageTextures = [stim.imageTexture_1FR, stim.imageTexture_50Cent, stim.imageTexture_20Cent, stim.imageTexture_0FR];





% Maximum priority level so that all the computer power is given to this
topPriorityLevel = MaxPriority(scr.window);
Priority(topPriorityLevel);

%% Stimulus Presentation Training - Reward/Punishment - Neutral

% Launch training
[all,stim] = training(scr,stim,speed,sound);

% Launch Reward/Punishment
for blockType = 1:2
    if stim.isPositiveStimuli == true
        
        % Reward Protocol followed by showing the results of the trial
        for i = 1:nbTrialPerCoinType*3
            stim.imageTextureIdx = stim.incentiveIdx(i);
            [doWin,signal,firstT2] = stimulusPresentation(scr,stim,speed,sound);
            all.signals{2,i} = signal;
            all.win{2,i}  = doWin;
            all.VCmax{2,i} = max(signal);
            all.trialLength{2,i} = length(signal);
            all.incentive{2,i} = stim.incentive(stim.incentiveIdx(i));
            all.firstT2{2,i} = firstT2;
            totalMoney = showResult(scr,stim,speed,totalMoney,doWin,stim.incentive(stim.incentiveIdx(i)));
        end
        
        % Select wanted variable and Launch a big break after a block
        stim.isPositiveStimuli = false;
        speed.isShortBreak = false;
        showBreak(scr,stim,speed,totalMoney)
        
        % Measure maximum power again ! save the data in the all struct.
        [measuredMVC,signal] = measureMVC(scr,stim,speed);
        all.measuredMVC{2,1} = measuredMVC;
        all.MVCsignal{2,1} = signal(1,:);
        all.MVCsignal{2,2} = signal(2,:);
        all.MVCsignal{2,2} = signal(3,:);
        stim.measuredMVC = max(measuredMVC);
        all.usedMVC{2,1} = stim.measuredMVC;
        
        % If it is the first block, do 10 sec break before launching next block
        if blockType == 1
            speed.isShortBreak = true;
            showBreak(scr,stim,speed,totalMoney)
        end
        
    elseif stim.isPositiveStimuli == false
        
        % Launch Punishment protocol and show results of the trial
        for i = 1:nbTrialPerCoinType*3
            stim.imageTextureIdx = stim.incentiveIdx(i);
            [doWin,signal,firstT2] = stimulusPresentation(scr,stim,speed,sound);
            all.signals{3,i} = signal;
            all.win{3,i}  = doWin;
            all.VCmax{3,i} = max(signal);
            all.trialLength{3,i} = length(signal);
            all.incentive{3,i} = stim.incentive(stim.incentiveIdx(i));
            all.firstT2{3,i} = firstT2;
            totalMoney = showResult(scr,stim,speed,totalMoney,doWin,stim.incentive(stim.incentiveIdx(i)));
        end
        
        % Select wanted variable and Launch a big break after a block
        stim.isPositiveStimuli = true;
        speed.isShortBreak = false;
        showBreak(scr,stim,speed,totalMoney)
        
        % Measure maximum power again ! save the data in the all struct.
        [measuredMVC,signal] = measureMVC(scr,stim,speed);
        all.measuredMVC{3,1} = measuredMVC;
        all.MVCsignal{3,1} = signal(1,:);
        all.MVCsignal{3,2} = signal(2,:);
        all.MVCsignal{3,2} = signal(3,:);
        stim.measuredMVC = max(measuredMVC);
        all.usedMVC{3,1} = stim.measuredMVC;
        if blockType == 1
            speed.isShortBreak = true;
            showBreak(scr,stim,speed,totalMoney)
        end
    end
    
end

% Select wanted variable and Launch a big break after a block
speed.isShortBreak = false;
showBreak(scr,stim,speed,totalMoney)

% Measure maximum power again ! save the data in the all struct one last time
[measuredMVC,signal] = measureMVC(scr,stim,speed);
all.measuredMVC{4,1} = measuredMVC;
all.MVCsignal{4,1} = signal(1,:);
all.MVCsignal{4,2} = signal(2,:);
all.MVCsignal{4,3} = signal(3,:);
stim.measuredMVC = max(measuredMVC);
all.usedMVC{4,1} = stim.measuredMVC;

%% neutral block
if strcmp(neutral_block_yn,'yes')
    [all] = LGCM_neutral_intrinsic_motiv_block(all, scr, stim, speed, sound, totalMoney, nbTrialPerCoinType);
end

%% Clear the screen
sca;
%% Save Data

% Double save to finish. one .mat (easier to take and do ML) and one excel just in case move them in the saving folde
save([savePath, subjectCodeName '_data.mat'],'-struct','all')
saveDataExcel(all,nbTrialPerPhase,subjectCodeName);