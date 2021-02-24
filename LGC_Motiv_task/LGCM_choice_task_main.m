% main task of the experiment:
% 1) fixation cross
% 2) choice between 2 options
% 3) display chosen option
% 4) perform the effort corresponding to the chosen option
%
% The effort can be a physical effort of a fixed intensity but a 
% varying duration (in the physical version of the task) or a mental effort
% with a fixed number 
%
% developed by Arthur Barakat & Nicolas Clairis - 2020/2021
%
% See also LGCM_ScreenConfiguration.m, LGCM_MVC_measurement.m, etc.

%% Clear the workspace and the screen, instrreset resets the udp channels
sca; % close all PTB screens
close all; % close all windows
clearvars; % clear variables from memory
instrreset; % Disconnect and delete all instrument objects

%% define if you are currently testing the script (1)
% (no need to have correct timings and everything in PTB)
% or if this is the actual experiment => use optimal timings of the
% computer (0)
testing_script = 1;

%% working directories
% main_folder = ['D:' filesep 'Matlab_codes' filesep]; % Arthur's laptop path
main_folder = ['C:',filesep,'Users',filesep,'Loco',filesep,...
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
    sub_nber_str,...
    session_nm] = deal('');

% effort type?
% while ~ismember(effort_type_letter,{'p','m'})
%     effort_type_letter = input('effort type? For physical, press ''p'', for mental press ''m''.','s');
% end
% switch effort_type_letter
%     case 'p'
%         effort_type = 'physical';
%     case 'm'
%         effort_type = 'mental';
% end
effort_type = 'mental';

% subject
% while isempty(sub_initials)
%     sub_initials = input('Subject initials?','s');
% end
% while isempty(sub_nber) || length(sub_nber) ~= 3
%     sub_nber_str = input('Subject id number? (3 numbers)','s'); % string
% %     number of the subject
% end
sub_initials = 'NC';
sub_nber_str = '000';
% convert subject number string into a number
i_sub = str2double(sub_nber_str);
% Create subjectCodeName which is used as a file saving name
subjectCodeName = strcat(sub_initials,'_s',sub_nber_str);


% ask session number (especially for fMRI mapping with behavioral data)
% n_sess_max = 10; % maximal number of sessions
% while isempty(session_nm) ||...
%         str2double(session_nm) < 1 ||...
%         str2double(session_nm) > n_sess_max
%     session_nm = input('Session number? (write 0 if training session)','s');
% end
session_nm = '00';

if ismember(session_nm,{'0','00'})
    learning_done = false;
elseif str2double(session_nm) > 0
    learning_done = true;
else
    error('problem with session_nm specification');
end

% file name
file_nm = [subjectCodeName,'_session',session_nm,'_',effort_type,'_task'];
% verify the files do not already exist
if exist([savePath, file_nm,'.mat'],'file')
    error(['The file name ',file_nm,'.mat already exists.',...
        ' Please relaunch with a new file name or delete the previous data.']);
end

%% fMRI/behavioral version of the task?
task_type = 'behavioral'; % 'behavioral'/'fMRI' (leave by default to avoid answering for every subject)

% depending on the version of the task (fMRI/behavioral),
% adapt the default options accordingly
switch task_type
    case 'behavioral' % behavioral version (Arthur protocol)
        IRM = 0; % will not include fMRI
    case 'fMRI' % fMRI version (Nicolas protocol)
        IRM = 1; % will include fMRI => record TTL pulses + adapt key recordings accordingly
end

%% include punishment condition?
punishment_yn = 'yes'; % include punishment trials?

%% if physical effort task => requires the BioPac
% (mental effort task might also require the BioPac)
switch effort_type
    case 'physical'
        BioPac_yn = 'yes';
        warning(['check this out ',...
            'https://fr.mathworks.com/help/daq/getting-started-with-session-based-interface-using-ni-devices.html',...
            ' for NI card compatibility']);
    case 'mental'
        BioPac_yn = 'no';
end

%% Open a connection to the biopack datastream
% (for measuring physical force + maybe also other variables like Skin
% Conductance or heart rate)
if strcmp(BioPac_yn,'yes')
    warning('update what is necessary for grip based on NI device requirements');
%     % add biopac functions if missing so that they can be used for grip
%     % acquisition
%     addpath(BioPac_folder); % probably useless once you switch to NI card
%     
%     % start BioPac recording
%     [u_out] = BioPac_start();
end

%% initialize screen
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = LGCM_ScreenConfiguration(IRM, testing_script);
white = scr.colours.white;
% black = scr.colours.black;

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
% calibration
switch effort_type % in case you use different numbers for each effort type
    case 'mental'
        n_calibTrials = 5;
    case 'physical'
        n_calibTrials = 5;
end
% actual task
n_trials = 90;
n_R_levels = 4;
n_E_levels = 4;

% stimulus related variables for the display
[stim] = LGCM_stim_initialize(scr, n_R_levels, n_E_levels, pics_folder);

% determine reward and punishment trials
choice_opt = LGCM_choice_option_design(n_R_levels, n_E_levels, punishment_yn, n_trials);
if strcmp(punishment_yn,'yes')
    R_or_P = choice_opt.R_or_P;
end

%% timings
t_cross = 0.5;
t_choice = 5;
t_dispChoice = 1.5;
switch effort_type % in case you use different numbers for each effort type
    case 'mental'
        t_max_effort = 5; % time to perform the task
    case 'physical'
        t_max_effort = 5; % time to perform the task
end
calibTimes.instructions = 5;
calibTimes.effort_max = t_max_effort;
calibTimes.fbk = 2;

%% learning + calibration
switch effort_type
    case 'mental' % for mental effort, define all the number sequences in advance
        mentalE_prm = LGCM_mental_effort_parameters(i_sub);
        mentalE_prm.startAngle = 0; % for learning always start at zero
        
        %% LGCM_mental_learning
        learning_time_limit = false;
        % extract numbers to use for each learning phase
        [numberVector_learning] = LGCM_mental_numbers(2);
        
        %% adapt number of required correct answers according to if learning was done or not
        if ~learning_done
            n_max_learning_withInstructions = 6;%20;
            n_max_learning_withoutInstructions = 16;%30;
            % calibration max
            n_calibMax = 4;
        elseif learning_done % short retraining at the beginning of each block
            n_max_learning_withInstructions = 6;%10;
            n_max_learning_withoutInstructions = 6;%10;
        end % learning + training session
        
        % perform a first session with instructions (left/right vs odd/even
        % and lower/higher than 5 mapping)
        instructions_disp = 1;
        % display instructions
        [onset_endInstru_withInstru] = LGCM_mental_learning(scr, instructions_disp, mentalE_prm);
        % perform the learning with instructions
        [mentalE_perf_learning_withInstru] = LGCM_mental_effort_perf(scr, stim, key,...
            numberVector_learning(1,:),...
            mentalE_prm, n_max_learning_withInstructions, instructions_disp, learning_time_limit);
        
        % do a second learning session without instructions (= no mapping
        % between side and answer type)
        instructions_disp = 0;
        % display instructions
        [onset_endInstru_withoutInstru] = LGCM_mental_learning(scr, instructions_disp, mentalE_prm);
        % perform the learning without instructions as in the real task
        [mentalE_perf_learning_withoutInstru] = LGCM_mental_effort_perf(scr, stim, key,...
            numberVector_learning(2,:),...
            mentalE_prm, n_max_learning_withoutInstructions, instructions_disp, learning_time_limit);
        
        %% max performance measurement
        if ~learning_done
            % extract numbers to use for each calibration trial
            [numberVector_calib] = LGCM_mental_numbers(n_calibTrials);
            n_mental_max_perTrial = LGCM_mental_calib(scr, stim, key,...
                numberVector_calib, mentalE_prm, n_calibTrials,...
                n_calibMax, calibTimes);
        elseif learning_done
            % indicate max performance
            n_mental_max_perTrial = input('please write manually n_max of pairs done during training');            
        end % learning + training session
        
        %% randomize the order of the numbers appearing on screen
        mental_nbers_per_trial = LGCM_mental_numbers(n_trials);
        
        %% randomize the type of the first trial (odd/even or higher/lower
        % than 5)
        mental_taskType_trialStart = LGCM_mental_task_start(n_trials);
        
    case 'physical'% for physical effort, ask the MVC
        
        % take an initial MVC measurement
        
        % store global MVC
        if exist([],'file') % retrieve past MVC value if has been extracted already
            
            
        else % perform MVC measurement otherwise
            
        end
        % [MVC_initial, onsets_MVC_initial] = LGCM_MVC_measurement(scr, stim, session_effort_type, n_MVC_repeat);
        warning('MVC measurement function needs to be updated given last task changes');
        
end

% %%
% warning('For monetary displays: add 4th layer after RGB colour values for transparency');
% %%

%% initialize onsets for main task
[onsets.cross,...
    onsets.dispChoiceOptions,...
    onsets.choice,...
    onsets.effortPeriod,...
    onsets.keyReleaseMessage,...
    onsets.cross_after_buttonRelease,...
    onset.fbk, onset.fbk_win, onset.fbk_loss,...
    onset.fbk_fail] = deal(NaN(1,n_trials));

%% start recording fMRI TTL
if IRM == 1
    dummy_scan = 4; % number of TTL to wait before starting the task
    trigger_id = 84; % trigger value corresponding to the TTL code
    [T0, TTL] = LGCM_keyboard_check_start(dummy_scans, trigger_id, key);
end % fMRI check

%% launch main task
choice = zeros(1,n_trials);
[R_chosen, E_chosen] = deal(NaN(1, n_trials));
was_a_key_pressed_bf_trial = NaN(1,n_trials);
for iTrial = 1:n_trials
    
    %% fixation cross
    Screen('FillRect',window,white, stim.cross.verticalLine); % vertical line
    Screen('FillRect',window,white, stim.cross.horizontalLine); % horizontal line
    [~,timeCrossNow] = Screen('Flip',window); % display the cross on screen
    onsets.cross(iTrial) = timeCrossNow;
    WaitSecs(t_cross);
    
    %% check that no key is being pressed before the choice trial starts
    [was_a_key_pressed_bf_trial(iTrial),...
        onsets.keyReleaseMessage(iTrial)] = LGCM_check_keys_are_up(scr, key);
    
    % if a key was pressed before starting the trial => show the fixation
    % cross again with a similar amount of time
    if was_a_key_pressed_bf_trial(iTrial) == 1
        Screen('FillRect',window,white, stim.cross.verticalLine); % vertical line
        Screen('FillRect',window,white, stim.cross.horizontalLine); % horizontal line
        [~,timeCrossNow] = Screen('Flip',window); % display the cross on screen
        onsets.cross_after_buttonRelease(iTrial) = timeCrossNow;
        WaitSecs(t_cross);
    end
    
    %% choice phase
    [choice(iTrial),...
        onsets.dispChoiceOptions(iTrial),...
        onsets.choice(iTrial),...
        stoptask] = LGCM_choice_period(scr, stim,...
        choice_opt, iTrial, t_choice, key);
    
    %% check if escape was pressed => stop everything if so
    if stoptask == 1
        % save all the data in case you still want to analyze it
        save([savePath, file_nm,'_earlyEnd_tmp.mat'],'-struct','all');
        break;
    end

    %% display chosen option only
    [time_dispChoice,...
        R_chosen(iTrial),...
        E_chosen(iTrial)] = LGCM_choice_task_dispChosen(scr, stim, choice_opt,...
        choice(iTrial),...
        iTrial);
    onsets.dispChoice(iTrial) = time_dispChoice;
    WaitSecs(t_dispChoice);
    
    %% perform the effort if a choice was made, otherwise punish the
    % subject for not answering to the choice
    if choice(iTrial) == 0 % no choice was made => failure
        trial_failed = 1;
        
    elseif ismember(choice(iTrial), [-1,1]) % choice done => perform the corresponding effort
        
        
        %% perform the effort
        switch effort_type
            case 'physical'
                [trial_failed] = LGCM_physical_effort();
            case 'mental'
                [trial_failed, perf, onsetEffortPeriod(iTrial)] = LGCM_mental_effort(scr, stim,...
                    t_effort_max,...
                    R_chosen, R_or_P{iTrial}, E_chosen, n_E_levels,...
                    mental_nbers_per_trial(iTrial,:),...
                    switchPerc,...
                    mental_taskType_trialStart(iTrial));
                %%
                warning('update things here based on what has been coded already for training and calibration');
                %%
        end
        
        %     %% if the effort was not performed correctly,
        %     % repeat the trial without the reward until the effort is achieved
        %     % note the timings and performance for each repetition
        %     %
        %     % OR: new version: punish them with a loss
        %     if effortDoneTrial == 0
        % %         while ~effortDoneTrial
        % %
        % %         end % effort done loop
        %     end % if effort done
    end
    %%
    switch trial_failed
        case 0 % trial is a success
            
            % display monetary incentive
            
            % display feedback
            switch R_or_P{iTrial}
                case 'R'
%                     DrawFormattedText(window,'You won', xScreenCenter, (1/5)*yScreenCenter);
                    DrawFormattedText(window,'Vous avez obtenu', xScreenCenter, (1/5)*yScreenCenter);
                case 'P'
%                     DrawFormattedText(window,'You lost', xScreenCenter, (1/5)*yScreenCenter);
                    DrawFormattedText(window,'Vous avez perdu', xScreenCenter, (1/5)*yScreenCenter);
            end
            
            [~,onset.fbk(iTrial)] = Screen(window,'Flip');
            switch R_or_P{iTrial}
                case 'R'
                    onset.fbk_win(iTrial) = onset.fbk(iTrial);
                case 'P'
                    onset.fbk_loss(iTrial) = onset.fbk(iTrial);
            end
        case 1 % trial failed
            % display message
%             DrawFormattedText(window,'Too slow!', xScreenCenter, (1/5)*yScreenCenter);
            DrawFormattedText(window,'Trop lent!', xScreenCenter, (1/5)*yScreenCenter);
            % display money loss for failing
            
            [~,onset.fbk_fail(iTrial)] = Screen(window,'Flip');
    end
    
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
% timings
all.times.cross         = t_cross;
all.times.choice        = t_choice;
all.times.dispChosen    = t_dispChoice;
all.times.effort_max    = t_max_effort;

% store timings in all structure
all.onsets = onsets;

% Double save to finish: .mat and .csv format
save([savePath, file_nm,'.mat'],'-struct','all');
saveDataExcel(all, nbTrialPerPhase, file_nm);

save([savePath, file_nm,'_messyAllStuff.mat']);