% Clear the workspace and the screen, instrreset resets the udp channels
sca; % close all PTB screens
close all; % close all windows
clearvars; % clear variables from memory
instrreset; % Disconnect and delete all instrument objects

%% working directories
main_folder                 = ['D:' filesep 'Matlab_codes' filesep];
main_task_folder            = [main_folder, 'LGC_Motiv_task' filesep];
savePath                    = [main_folder, 'LGC_Motiv_results' filesep];
BioPac_folder               = [main_folder, 'BioPac_functions' filesep];
pics_folder                 = [main_task_folder, 'Coin_PNG', filesep];
Matlab_DIY_functions_folder = [main_folder, 'Matlab_DIY_functions', filesep];

% add personal functions (needed for PTB opening at least)
addpath(Matlab_DIY_functions_folder);

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

%% session number (crucial for fMRI mapping btw fMRI session and behavioral data)
session_nm = input('Session number?');

%% task type (physical/mental effort)


%% Initialize main variables
totalMoney = 20; % initial monetary amount for participating in the study(?)
warning('Is the comment for the variable totalMoney correct Arthur?');
nbTrialPerCoinType = 1; % what is this?
nbTrialPerPhase = [3 nbTrialPerCoinType*3 nbTrialPerCoinType*3 nbTrialPerCoinType]; % what is this?

% fMRI/behavioral version?
task_type = questdlg('behavioral or fMRI version?',...
    'task version',...
    'behavioral','fMRI','behavioral');

% depending on the version (fMRI/behavioral), adapt the default options
% accordingly
switch task_type
    case 'behavioral' % behavioral version (Arthur protocol)
        IRM = 0;
        
        neutral_block_yn = 'yes'; % include neutral block at the end for intrinsic motivation
        audio_fbk_yn = 'yes'; % include audio feedback (?)
    case 'fMRI' % fMRI version (Nicolas protocol)
        IRM = 1;
        
        neutral_block_yn = 'no'; % include neutral block at the end for intrinsic motivation
        audio_fbk_yn = 'no'; % include audio feedback (?)
end

% physical or mental effort task

% physical effort task => requires the BioPac
% (mental effort task might also require the BioPac)
BioPac_yn = 'yes';

BioPac_yn = 'no';

file_nm = [subjectCodeName,'_session',];

%% check subject, run and task are new (= no overwritting of previous files)
if exist([savePath, file_nm,'.mat'],'file')
    error(['The file name ',file_nm,'.mat already exists. Please relaunch with a new file name or delete the data.']);
end
if exist([savePath, file_nm,'.xlsx'],'file')
    error(['The file name ',file_nm,'.xlsx already exists. Please relaunch with a new file name or delete the data.']);
end

%% start PTB
% Screen variables
% PTB initialization + anti bug due to sync
% Here we call some default settings for setting up Psychtoolbox. So it doesn't crash when we call it
[xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(IRM, testing_script);

%% check this is equivalent to WindowSize in pixels (Arthur was extracting both as if they were different measures)
[scr.xCenter, scr.yCenter] = RectCenter(scr.windowRect);
%%
% Draw to the external screen if avalaible
scr.screenNumber = max(scr.screens);
stim.VSsignal = 0;

% Open an on screen window
[scr.window, scr.windowRect] = PsychImaging('OpenWindow', scr.screenNumber, scr.grey);

% Enable alpha blending for anti-aliasing of all our textures
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% store main info in scr structure
scr.window = window;
% Get the centre coordinate of the window
scr.xCenter = xScreenCenter;
scr.yCenter = yScreenCenter;

%% Speed related variables
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
scr.squareY = yScreenCenter;

% Position of each ring/coin/signal on the screen
stim.centeredthreshold_1 = CenterRectOnPointd(stim.threshold_1,xScreenCenter, yScreenCenter);
stim.centeredthreshold_2 = CenterRectOnPointd(stim.threshold_2,xScreenCenter, yScreenCenter);
stim.missThreshold = CenterRectOnPointd(stim.missThreshold,xScreenCenter, yScreenCenter);
stim.missTarget = CenterRectOnPointd(stim.missTarget,scr.windowRect(3) - scr.squareX, yScreenCenter);

%Import the images
[image_20Cent,~,alpha_20Cent] = imread([pics_folder 'SMT_Coin_20RP.png']);
[image_50Cent,~,alpha_50Cent] = imread([pics_folder 'SMT_Coin_50RP.png']);
[image_1FR,~,alpha_1FR] = imread([pics_folder 'SMT_Coin_1FR.png']);
[image_0FR,~,alpha_0FR] = imread([pics_folder 'SMT_Coin_0FR.png']);

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

%% PTB sound
if strcmp(audio_fbk_yn,'yes')
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
end

%% Open a connection to the biopack datastream
if strcmp(BioPac_yn,'yes')
    % add biopac functions if missing so that they can be used for grip
    % acquisition
    addpath(BioPac_folder);
    % start BioPac recording
    [stim] = BioPac_start();
end

%% perform the task physical/mental effort

%% Launch training
[all,stim] = training(scr,stim,speed,sound);

%% Launch Reward/Punishment
for iBlockType = 1:2
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
        if iBlockType == 1
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
        if iBlockType == 1
            speed.isShortBreak = true;
            showBreak(scr,stim,speed,totalMoney)
        end
    end
    
end

% Select wanted variable and Launch a big break after a block
speed.isShortBreak = false;
showBreak(scr,stim,speed,totalMoney)

%% Measure maximum power again ! save the data in the all struct one last time
[measuredMVC,signal] = measureMVC(scr,stim,speed);
all.measuredMVC{4,1} = measuredMVC;
all.MVCsignal{4,1} = signal(1,:);
all.MVCsignal{4,2} = signal(2,:);
all.MVCsignal{4,3} = signal(3,:);
stim.measuredMVC = max(measuredMVC);
all.usedMVC{4,1} = stim.measuredMVC;

%% neutral block
if strcmp(neutral_block_yn,'yes')
    [all] = LGCM_neutral_intrinsic_motiv_block(all,...
    scr, stim, speed, sound, totalMoney,...
    pause_dur,...
    nbTrialPerCoinType);
end

%% shouldn't there be a function to close BioPac stuff?

%% Clear the PTB screen
sca;

%% Save Data
% Double save to finish: .mat and .csv format
save([savePath, file_nm,'.mat'],'-struct','all');
saveDataExcel(all, nbTrialPerPhase, file_nm);