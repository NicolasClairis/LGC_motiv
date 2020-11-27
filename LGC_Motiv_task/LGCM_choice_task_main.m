% main task of the experiment:
% 1) fixation cross
% 2) choice between 2 options
% 3) display chosen option
% 4) perform the effort corresponding to the chosen option
%
% developed by Arthur Barakat & Nicolas Clairis - 2020/2021
%
% See also LGCM_ScreenConfiguration.m, LGCM_MVC_measurement.m, etc.

%% Clear the workspace and the screen, instrreset resets the udp channels
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

% create results folder if no subject has been acquired yet
if ~exist(savePath,'dir')
    mkdir(savePath);
end

%% task type (physical/mental effort) + subject + session identification
[effort_type_letter,...
    sub_initials,...
    sub_nber,...
    session_nm] = deal('');

% effort type?
while ~ismember(effort_type_letter,{'p','m'})
    effort_type_letter = input('effort type? For physical, press ''p'', for mental press ''m''.','s');
end
switch effort_type_letter
    case 'p'
        effort_type = 'physical';
    case 'm'
        effort_type = 'mental';
end

% subject
while isempty(sub_initials)
    sub_initials = input('Subject initials?','s');
end
while isempty(sub_nber)
    sub_nber = input('Subject id number?','s');
end
% Create subjectCodeName which is used as a file saving name
subjectCodeName = strcat(sub_initials,'_s',sub_nber);

% ask session number (especially for fMRI mapping with behavioral data)
n_sess_max = 10; % maximal number of 
while isempty(session_nm) ||...
        str2double(session_nm) < 1 ||...
        str2double(session_nm) > n_sess_max
    session_nm = input('Session number?','s');
end

% file name
file_nm = [subjectCodeName,'_session',session_nm,'_',effort_type,'_task'];
% verify the files do not already exist
if exist([savePath, file_nm,'.mat'],'file')
    error(['The file name ',file_nm,'.mat already exists.',...
        ' Please relaunch with a new file name or delete the previous data.']);
end

%% fMRI/behavioral version of the task?
task_type = 'fMRI'; % 'behavioral'/'fMRI' (leave by default to avoid answering for every subject)

% depending on the version of the task (fMRI/behavioral),
% adapt the default options accordingly
switch task_type
    case 'behavioral' % behavioral version (Arthur protocol)
        IRM = 0; % will not include fMRI
        punishment_yn = 'yes'; % include punishment trials?
    case 'fMRI' % fMRI version (Nicolas protocol)
        IRM = 1; % will include fMRI => record TTL pulses + adapt key recordings accordingly
        punishment_yn = 'no'; % include punishment trials?
end

%% if physical effort task => requires the BioPac
% (mental effort task might also require the BioPac)
switch effort_type
    case 'physical'
        BioPac_yn = 'yes';
    case 'mental'
        BioPac_yn = 'no';
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

%% initialize screen
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = LGCM_ScreenConfiguration(IRM, testing_script);

%% stimulus related variables
[stim] = LGCM_stim_initialize(scr);

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

%% task parameters
n_trials = 90;
R_levels = [];
E_levels = [];
% create all possible combinations of these levels
warning('create all possible combinations of these levels');
% make as many mini-blocks containing all those combinations as there are
% trials
warning('make as many mini-blocks containing all those combinations as there are trials');
% randomize the size of each option across trials
warning('randomize the size of each option across trials');

%% timings
t_cross = 0.5;
t_choice = 3;
t_dispChoice = 1.5;
t_effort = 5; % stable total time, but duration of the effort will depend on the difficulty level at each trial

%% start recording fMRI TTL
if IRM == 1
    dummy_scan = 4; % number of TTL to wait before starting the task
    trigger_id = 84; % trigger value corresponding to the TTL code
    [T0, TTL] = LGCM_keyboard_check_start(dummy_scans, trigger_id);
end % fMRI check

%% perform the task
%% initial MVC measurement
[MVC_initial, onsets_MVC_initial] = LGCM_MVC_measurement(scr, session_effort_type, speed, stim, n_MVC_repeat);

%% Launch training trials
[onsets_training] = LGCM_training(scr, stim, speed,...
    audio_fbk_yn, sound,...
    session_effort_type, n_training_trials);

%% launch main task
for iTrial = 1:n_trials
    
    % re-initialize choice
    choice = 0;
    
    %% fixation cross
    Screen('FillRect',window,white,stim.cross_verticalLine); % vertical line
    Screen('FillRect',window,white,stim.cross_horizontalLine); % horizontal line
    [~,timeCrossNow] = Screen('Flip',window); % display the cross on screen
    onsets.cross(iTrial) = timeCrossNow;
    WaitSecs(t_cross);
    
    %% choice task
    
    [~,timeChoicePeriodNow] = Screen('Flip',window);
    onsets.choice_period(iTrial) = timeChoicePeriodNow;
    while timeNow <= timeChoicePeriodNow + t_choice
        [keyisdown, secs, keycode] = KbCheck;
        if keyisdown == 1 &&...
                keycode(key.left) == 1 && keycode(key.right) == 0
            timedown = GetSecs;
            choice = -1;
        elseif keyisdown == 1 &&...
                keycode(key.left) == 0 && keycode(key.right) == 1
            timedown = GetSecs;
            choice = 1;
        end
        onsets.choice(iTrial) = timedown;
        warning('here you can either leave it that way (no visual feedback at the moment they chose an option',...
            ' or show some visual feedback of the selected option already',...
            ' or require the subject to keep pressing the button until the end of the trial.');
    end % choice period
    
    %% display chosen option only
    
    [~,timeDispChoiceNow] = Screen('Flip',window);
    onsets.dispChoice(iTrial) = timeDispChoiceNow;
    WaitSecs(t_dispChoice);
    
    %% perform the effort
    % extract reward (or punishment) and effort levels according to choice
    % made
    switch choice
        case -1 % left option selected
            R_level_tmp = R_level_left;
            E_level_tmp = E_level_left;
        case 1 % right option selected
            R_level_tmp = R_level_right;
            E_level_tmp = E_level_right;
        case 0 % no choice made
            R_level_tmp = 0; % no reward (or punishment? we have to decide)
            warning('no reward or punishment when subject does not chose anything? We need to decide.');
            E_level_tmp = max(E_level_left, E_level_right);
    end
    
    [~,timeEffortPeriodNow] = Screen('Flip',window);
    onsets.startEffortPeriod(iTrial) = timeEffortPeriodNow;
    
    
    % update the effort display with time
    % default money size if no effort is being made
    [~,timeEffortPeriodNow] = Screen('Flip',window);
    % store money size and timing for each timestep
    
    % increase of money size if effort > threshold with time
    [~,timeEffortPeriodNow] = Screen('Flip',window);
    % store money size and timing for each timestep
    
    % decrease of money size if effort < threshold while was above
    [~,timeEffortPeriodNow] = Screen('Flip',window);
    % store money size and timing for each timestep
    
end % trial loop

%% Measure maximum power again at the end
[MVC_last, onsets_MVC_last] = LGCM_MVC_measurement(scr, session_effort_type, speed, stim, n_MVC_repeat);

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