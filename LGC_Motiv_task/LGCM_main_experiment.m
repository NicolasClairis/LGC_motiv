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
% Prepare an array of incentive, each type of coin being shown the same amount of time in a randomized order
incentiveIdx = cat(1,...
    repelem([1 2 3],...
    nbTrialPerCoinType*3));
incentiveIdx = incentiveIdx(randperm(nbTrialPerCoinType*3));

%% timings


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
    window, baselineTextSize] = ScreenConfiguration(IRM, testing_script);

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
    [stim] = BioPac_start();
end

%% perform the task
%% initial MVC measurement
[MVC, onsets_MVC] = LGCM_MVC_measurement(scr, session_effort_type, speed, stim);
% store all in final output variable
all.MVC(1) = MVC.MVC; % store initial MVC
all.MVC_allValues_perTest(:,1) = MVC.MVC_perCalibSession;
onsets.initial_MVC = onsets_MVC;

%% Launch training trials
[all, stim] = LGCM_training(scr, stim, speed,...
    audio_fbk_yn, sound,...
    session_effort_type);

%% Launch main task
% reward block


% punishment block (always performed after reward block if included so that
% it does not influence reward block parameters)
if strcmp(punishment_block_yn,'yes')
    
end

%% Measure maximum power again ! save the data in the all struct one last time
warning('re-use same function but replace variable names here');

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