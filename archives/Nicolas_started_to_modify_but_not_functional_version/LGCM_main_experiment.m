% Clear the workspace and the screen, instrreset resets the udp channels
sca; % close all PTB screens
close all; % close all windows
clearvars; % clear variables from memory
instrreset; % Disconnect and delete all instrument objects

%% working directories
% main_folder                 = ['D:' filesep 'Matlab_codes' filesep]; % Arthur's laptop path
main_folder                 = ['C:',filesep,'Users',filesep,'Loco',filesep,...
    'Documents',filesep,'GitHub',filesep,...
    'LGC_motiv', filesep]; % Nico's laptop Github path
main_task_folder            = [main_folder, 'LGC_Motiv_task' filesep];
savePath                    = [main_folder, 'LGC_Motiv_results' filesep];
BioPac_folder               = [main_folder, 'BioPac_functions' filesep];
pics_folder                 = [main_task_folder, 'Coin_PNG', filesep];
Matlab_DIY_functions_folder = [main_folder, 'Matlab_DIY_functions', filesep];

% add personal functions (needed for PTB opening at least)
addpath(genpath(main_task_folder));
addpath(BioPac_folder);
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

%% physical or mental effort task
session_effort_type = questdlg('Which task: physical or mental effort?',... % question
    'effort type',... % window title
    'physical','mental',... % options
    'mental'); % default option

%% if physical effort task => requires the BioPac
% (mental effort task might also require the BioPac)
switch session_effort_type
    case 'physical'
        BioPac_yn = 'yes';
    case 'mental'
        BioPac_yn = 'no';
end

%% fMRI/behavioral version of the task?
task_type = questdlg('behavioral or fMRI version?',...
    'task version',...
    'behavioral','fMRI','behavioral');

% depending on the version of the task (fMRI/behavioral),
% adapt the default options accordingly
switch task_type
    case 'behavioral' % behavioral version (Arthur protocol)
        IRM = 0; % will not include fMRI
        punishment_block_yn = 'yes'; % include punishment block after the reward block?
        neutral_block_yn = 'yes'; % include neutral block at the end for intrinsic motivation
        audio_fbk_yn = 'yes'; % include audio feedback (?)
    case 'fMRI' % fMRI version (Nicolas protocol)
        IRM = 1; % will include fMRI => record TTL pulses + adapt key recordings accordingly
        punishment_block_yn = 'no'; % include punishment block after the reward block?
        neutral_block_yn = 'no'; % include neutral block at the end for intrinsic motivation
        audio_fbk_yn = 'no'; % include audio feedback (?)
end

%% Initialize main variables
totalMoney = 20; % initial monetary amount for participating in the study(?)
warning('Is the comment for the variable totalMoney correct Arthur?');
nbTrialPerCoinType = 1; % what is this?
nbTrialPerPhase = [3,...
    nbTrialPerCoinType*3,...
    nbTrialPerCoinType*3,...
    nbTrialPerCoinType]; % what is this?

% define order of each incentive display
%
% METHOD 1: it will randomize the order across all trials
% ** advantage = reduces predictability for the next trial
% ** inconvenient = might create patterns inside the task with periods with
% a lot of high amount incentive and others with a lot of low amount
% incentive
%
% % Prepare an array of incentive, each type of coin being shown the same amount of time in a randomized order
% incentiveIdx = cat(1,...
%     repelem([1 2 3],...
%     nbTrialPerCoinType*3));
% incentive_rdm_order_idx = randperm(nbTrialPerCoinType*3);
% incentiveIdx = incentiveIdx(incentive_rdm_order_idx);
%
% METHOD 2: it will randomize the order within each mini-block. One
% mini-block contains all possible incentives one single time. 
% ** advantage: the frequency of each incentive is stable in the whole
% task => no high vs low reward periods
% ** inconvenient: if the number of incentives is low (2 or 3 typically)
% then it becomes more predictable (if you have seen 2/3 of the possible 
% incentives, then you might expect the 3rd one) => more efficient if you
% have more than 2/3 incentive possibilities


% training variables
n_training_trials = 3;

% MVC
n_MVC_measures = 2; % 1 initial measure + 1 at the end
n_MVC_repeat = 3; % number of repetitions of the MVC measurement
all.MVC = NaN(1,n_MVC_measures);
all.MVC_allValues_perTest = NaN(n_MVC_repeat, n_MVC_measures);

%% effort stimulus characteristics

%% timings
t_wait.training_instructions = 3;
t_wait.fixation_cross = 0.5;
% design jitter for period varying


%% file name
file_nm = [subjectCodeName,'_session',num2str(session_nm),'_',session_effort_type,'_task'];

%% check subject, run and task are new (= no overwritting of previous files)
% check mat file
if exist([savePath, file_nm,'.mat'],'file')
    error(['The file name ',file_nm,'.mat already exists.',...
        ' Please relaunch with a new file name or delete the data.']);
end
% check excel file
if exist([savePath, file_nm,'.xlsx'],'file')
    error(['The file name ',file_nm,'.xlsx already exists.',...
        ' Please relaunch with a new file name or delete the data.']);
end

%% PTB initialization
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = LGCM_ScreenConfiguration(IRM, testing_script);

%% Stimulus related variables
[stim] = LGCM_stim_initialize(scr);

% stimulus motion speed
[speed] = LGCM_stim_speed(window);

%% PTB sound
if strcmp(audio_fbk_yn,'yes')
    [sound] = LGCM_fbk_sound(main_task_folder);
else
    sound = [];
end

%% Open a connection to the biopack datastream
% (for measuring physical force + maybe also other variables like Skin
% Conductance or heart rate)
if strcmp(BioPac_yn,'yes')
    % add biopac functions if missing so that they can be used for grip
    % acquisition
    addpath(BioPac_folder);
    
    % start BioPac recording
    [u_out] = BioPac_start();
end

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
end

%% start recording fMRI TTL
%% TRIGGER RAPHAEL
if IRM == 1
    dummy_scan = 4; % number of TTL to wait before starting the task
    trigger = 53; % trigger value corresponding to TTL code
    [T0, TTL] = LGCM_keyboard_check_start(dummy_scans, trigger);
end % fMRI check

%% perform the task
%% initial MVC measurement
[MVC_initial, onsets_MVC_initial] = LGCM_MVC_measurement(scr, session_effort_type, speed, stim, n_MVC_repeat);

%% Launch training trials
[onsets_training] = LGCM_training(scr, stim, speed,...
    audio_fbk_yn, sound,...
    session_effort_type, n_training_trials);

%% Launch main task

%% Announce the next block
%�%
% reward block

% Select wanted variable and Launch a big break after a block
showBreak(scr,stim,speed,totalMoney)

% punishment block (always performed after reward block if included so that
% it does not influence reward block parameters)
if strcmp(punishment_block_yn,'yes')
    
end

%% Measure maximum power again ! save the data in the all struct one last time
[MVC_last, onsets_MVC_last] = LGCM_MVC_measurement(scr, session_effort_type, speed, stim, n_MVC_repeat);

%% neutral block
if strcmp(neutral_block_yn,'yes')
    [output_vars] = LGCM_neutral_intrinsic_motiv_block(all,...
    scr, stim, speed, sound, totalMoney,...
    pause_dur,...
    nbTrialPerCoinType);
end

%% get all TTL from the task
if IRM == 1
    switch session_effort_type
        case 'physical'
            TTL = LGCM_keyboard_check_end(TTL);
        case 'mental'
            [TTL, keyLeft, keyRight] = LGCM_keyboard_check_end(TTL);
            % key storage of when left/right key have been pressed
            all.keys.keyLeft    = keyLeft;
            all.keys.keyRight   = keyRight;
    end
    % store T0 and TTL timings in onsets structure
    onsets.T0 = T0;
    onsets.TTL = TTL;
end

%% shouldn't there be a function to close BioPac stuff?

%% Clear the PTB screen
sca;

%% Save Data
% store all relevant data in final output variable
all.MVC(1) = MVC_initial.MVC; % store initial MVC
all.MVC(2) = MVC_last.MVC; % store final MVC
all.MVC_allValues_perTest(:,1) = MVC_initial.MVC_perCalibSession;
all.MVC_allValues_perTest(:,2) = MVC_last.MVC_perCalibSession;
all.MVC_initial = MVC_initial;
all.MVC_last = MVC_last;
onsets.initial_MVC = onsets_MVC_initial;
onsets.final_MVC = onsets_MVC_final;

% store timings in all structure
all.onsets = onsets;

% Double save to finish: .mat and .csv format
save([savePath, file_nm,'.mat'],'-struct','all');
saveDataExcel(all, nbTrialPerPhase, file_nm);

save([savePath, file_nm,'_messyAllStuff.mat']);