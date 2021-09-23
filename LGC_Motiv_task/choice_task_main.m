% main task of the experiment:
% 1) fixation cross
% 2) choice between 2 options
% 3) display chosen option
% 4) perform the effort corresponding to the chosen option
%
% The effort can be a physical effort of a fixed intensity but a
% varying duration (in the physical version of the task) or a mental effort
% with a fixed maximal duration but a varying number of questions to solve
% according to difficulty (mental version of the task).
%
% developed by Arthur Barakat & Nicolas Clairis - 2020/2021
%
% See also ScreenConfiguration.m, physical_effort_perf.m,
% mental_effort_perf_Nback.m

%% Clear the workspace and the screen, instrreset resets the udp channels
ShowCursor;
sca; % close all PTB screens
close all; % close all windows
clearvars; % clear variables from memory
instrreset; % Disconnect and delete all instrument objects
clc;

%% define if you are currently testing the script (1)
% (no need to have correct timings and everything in PTB)
% or if this is the actual experiment => use optimal timings of the
% computer (0)
testing_script = 0;

%% langage to display instructions
langage = 'fr';

%% working directories
scripts_folderName = 'LGC_Motiv_task';
% define working directories
cd ..
main_folder                 = [pwd filesep];
main_task_folder            = [main_folder, scripts_folderName filesep];
results_folder              = [main_folder, 'LGC_Motiv_results' filesep];
% BioPac_folder               = [main_folder, 'BioPac_functions' filesep];
% pics_folder                 = [main_task_folder, 'Coin_PNG', filesep];
Matlab_DIY_functions_folder = [main_folder, 'Matlab_DIY_functions', filesep];

% add personal functions (needed for PTB opening at least)
addpath(genpath(main_task_folder));
% addpath(BioPac_folder);
addpath(Matlab_DIY_functions_folder);

% create results folder if no subject has been acquired yet
if ~exist(results_folder,'dir')
    mkdir(results_folder);
end
% go back to folder with scripts
cd(main_task_folder);

%% Define subject ID

% Insert the initials, the number of the participants
[iSubject, effort_type, session_nm] = deal([]);
while isempty(iSubject) ||...
        isempty(effort_type) || ~ismember(effort_type,{'p','m'}) ||...
        isempty(session_nm) || str2double(session_nm) < 0 % repeat until both are answered
    info = inputdlg({'Subject ID','Type d''effort (p/m)','Session number(0-4) (0 for initial calibration)'});
    [iSubject, effort_type, session_nm] = info{[1,2,3]};
end
session_nb = str2double(session_nm);
switch effort_type
    case 'p'
        effort_type = 'physical';
    case 'm'
        effort_type = 'mental';
end

% Create subjectCodeName which is used as a file saving name
subjectCodeName = strcat('_CID',iSubject);

% file name
file_nm = [subjectCodeName,'_session',session_nm,'_',effort_type,'_task'];
% verify the files do not already gexist
if exist([results_folder, file_nm,'.mat'],'file')
    error(['The file name ',file_nm,'.mat already exists.',...
        ' Please relaunch with a new file name or delete the previous data.']);
end

%% fMRI/behavioral version of the task?
IRM = 1;
% (0) does not include fMRI = training
% (1) include fMRI

%% include punishment condition?
punishment_yn = 'yes'; % include punishment trials?

%% task parameters
% initialize screen
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(IRM, testing_script);
white = scr.colours.white;
black = scr.colours.black;

% how many possible answers
n_buttonsChoice = 4;

% define relevant keys and dynamometer
switch effort_type
    case 'mental'
        key = relevant_key_definition('mental', IRM, n_buttonsChoice);
    case 'physical' % need dq output to record the handgrip data
        [key, dq] = relevant_key_definition('physical', IRM, n_buttonsChoice);
end

% initial calibration
switch effort_type
    case 'physical'
        n_calibTrials = 3;
    case 'mental'
        n_calibTrials = 3;
end

% calibration before/after end of each fMRI session
n_MaxPerfTrials = 2;

% actual task
n_R_levels = 3;
n_E_levels = 3;
nTrials = 48;

% extract money amount corresponding to each reward level for the
% computation of the gains
R_money = R_amounts(n_R_levels, punishment_yn);

% check trial number is ok based on the number of entered conditions
% you should have a pair number of trials so that you can define an equal
% amount of reward and punishment trials
if strcmp(punishment_yn,'yes') && mod(nTrials,2) ~= 0
    error(['you added punishments in the task, hence you need a pair number of total trials.',...
        ' Please fix this so that you can have the same number of trials in punishment and reward']);
end
%
% determine reward/punishment and effort level combinations for each trial
% choiceOptions = choice_option_design(n_R_levels, n_E_levels, punishment_yn, nTrials, R_money);
bestMatrix = getfield( load([main_task_folder,'DaBestDesignMat.mat'],'bestMatrix'),'bestMatrix');
choiceOptions = RP_moneyLevels(bestMatrix, R_money);
% display of different answers only during the training
confDispDuringChoice = false;

switch effort_type
    case 'physical'
        F_threshold = 55; % force should be maintained above this threshold (expressed in % of MVC)
        F_tolerance = 2.5; % tolerance allowed around the threshold (expressed in % of MVC)
        % need to define timings for each level of force
        [Ep_time_levels] = physical_effortLevels(n_E_levels);
    case 'mental'
        % define number of pairs to solve for each level of difficulty
        n_to_reach = mental_N_answersPerLevel(n_E_levels);
        % calibration: calibrate the maximal duration required for the
        % top effort
        n_calibMax = n_to_reach.(['E_level_',num2str(n_E_levels)]);
end

% stimulus related variables for the display
[stim] = stim_initialize(scr, n_E_levels, langage, R_money);
barTimeWaitRect = stim.barTimeWaitRect;

% define number of training conditions
switch punishment_yn
    case 'yes'
        trainingConditions = {'R','P','RP'};
    case 'no'
        trainingConditions = {'R'};
end
n_trainingConditions = length(trainingConditions);

%% load timings for each phase of the experiment
[~, calibTimes, ~, taskTimes, mainTimes] = timings_definition(trainingConditions, n_R_levels, n_E_levels, nTrials, effort_type);
t_endfMRI = mainTimes.endfMRI;
t_endSession = mainTimes.endSession;

%% calibration (before the MRI)
% calibration performance file name
calibPerf_file_nm = [results_folder,effort_type,'_calibPerf_sub',subjectCodeName,'.mat'];
if strcmp(effort_type,'mental')
    mentalE_prm_calib = mental_effort_parameters(iSubject);
    mentalE_prm_calib.startAngle = 0; % for learning always start at zero
    % no error threshold nor mapping of answers when errors are
    % made
    calib_errorLimits_Em.useOfErrorMapping = false;
    calib_errorLimits_Em.useOfErrorThreshold = true;
    calib_errorLimits_Em.errorThreshold = 20;
end
if session_nb == 0
    %% initial calibration Fmax and mental time (inside the scanner)
    switch effort_type
        case 'mental'
            %% max performance measurement
            % extract numbers to use for each calibration trial
            %             [numberVector_calib] = mental_numbers(n_calibTrials);
            [numberVector_calib] = mental_calibNumberVector(n_calibTrials, n_calibMax);
            
            %% alternatively, use fixed number of correct answers to provide for each effort
            % level
            % repeat calibration until the subject performance is better
            % than the requested time threshold
            calibSuccess = false;
            calibSession = 0;
            while calibSuccess == false
                calibSession = calibSession + 1;
                [t_min_calib, calibSessionSummary, calibSuccess] = mental_calibTime(scr, stim, key,...
                    numberVector_calib, mentalE_prm_calib, n_calibTrials, n_calibMax, calibTimes, calib_errorLimits_Em, langage);
                calibSummary.(['calibSession_',num2str(calibSession)]).calibSummary = calibSessionSummary;
                calibSummary.(['calibSession_',num2str(calibSession)]).calibSuccess = calibSuccess;
                calibSummary.(['calibSession_',num2str(calibSession)]).t_mental_max_perTrial = t_min_calib;
            end
            % save calib performance
            save(calibPerf_file_nm,'t_min_calib',...
                'calibSummary','calibSessionSummary','calibSuccess','calibSession',...
                'n_calibTrials','n_calibMax','numberVector_calib');
            
            
        case 'physical'% for physical effort, ask the MVC
            % record and store global MVC
            [calib_MVC, onsets_calib_MVC] = physical_effort_MVC(scr, stim, dq, n_calibTrials, calibTimes);
            MVC = calib_MVC.MVC; % expressed in Voltage
            % save calib performance
            save(calibPerf_file_nm,'MVC',...
                'onsets_calib_MVC','calib_MVC');
    end
    
elseif session_nb > 0
    %% actual task sessions in the scanner
    
    % load max performance from calibration
    switch effort_type
        case 'mental'
            t_min_calib = getfield(load(calibPerf_file_nm,'t_min_calib'),'t_min_calib');
            % increase time available based on calibration minimum time
            % performed to allow for more time
            t_max_effort = t_min_calib*taskTimes.t_min_scalingFactor; % allow more time then min performance
            taskTimes.max_effort = t_max_effort;
        case 'physical'
            MVC = getfield(load(calibPerf_file_nm,'MVC'),'MVC'); % expressed in Voltage
    end
    
    %% max perf measurement before start of each session
    % (not for training out of MRI)
    switch effort_type
        case 'mental'
            % perform max perf
            [numberVector_initialMaxPerf] = mental_calibNumberVector(n_MaxPerfTrials, n_calibMax);
            [t_min_initialMaxPerf, initialMaxPerfSessionSummary, initialMaxPerfSuccess] = mental_calibTime(scr, stim, key,...
                numberVector_initialMaxPerf, mentalE_prm_calib, n_MaxPerfTrials, n_calibMax, calibTimes, calib_errorLimits_Em, langage);
        case 'physical'
            % take an initial MVC measurement (even if it has been done in a
            % previous session, will allow us to keep track of the force level
            % of our participants)
            [initial_MVC, onsets_initial_MVC] = physical_effort_MVC(scr, stim, dq, n_MaxPerfTrials, calibTimes);
    end
    
    %% launch main task
    % define task parameters
    switch effort_type
        case 'mental'
            Ep_or_Em_vars.i_sub = iSubject;
            Ep_or_Em_vars.n_to_reach = n_to_reach;
            % for actual task: no display of mapping but consider 3
            % errors as a trial failure
            Ep_or_Em_vars.errorLimits.useOfErrorMapping = false;
            Ep_or_Em_vars.errorLimits.useOfErrorThreshold = false;
        case 'physical'
            Ep_or_Em_vars.MVC = MVC;
            Ep_or_Em_vars.dq = dq;
            Ep_or_Em_vars.Ep_time_levels = Ep_time_levels;
            Ep_or_Em_vars.F_threshold = F_threshold;
            Ep_or_Em_vars.F_tolerance = F_tolerance;
    end
    Ep_or_Em_vars.timeRemainingEndTrial_ONOFF = 0;
    
    %% launch physiological recording
    disp('Please start physiological recording and then press space.');
    [~, ~, keyCode] = KbCheck();
    while(keyCode(key.space) ~= 1)
        % wait until the key has been pressed
        [~, ~, keyCode] = KbCheck();
    end
    disp('OK - space was pressed, physio recording started');
    
    %% instruction that main task will start soon
    DrawFormattedText(window, stim.expWillStart.text,...
        stim.expWillStart.x, stim.expWillStart.y, scr.colours.white, scr.wrapat);
    [~, onsets.taskWillStart] = Screen(window, 'Flip');
    disp('Please press space and then launch fMRI (Be careful to respect this order for the T0...');
    [~, ~, keyCode] = KbCheck();
    while(keyCode(key.space) ~= 1)
        % wait until the key has been pressed
        [~, ~, keyCode] = KbCheck();
    end
    disp('OK - space was pressed');
    
    %% start recording fMRI TTL and wait for a given amount of TTL before
    % starting the task in order to calibrate all timings on T0
    if IRM == 1
        dummy_scans = 1; % number of TTL to wait before starting the task (dummy scans are already integrated in CIBM scanner)
        [T0, TTL] = keyboard_check_start(dummy_scans, key.trigger_id, key);
    end % fMRI check
    
    %% perform choice and performance task
    [perfSummary] = choice_and_perf(scr, stim, key,...
        effort_type, Ep_or_Em_vars, R_money,...
        'mainTask', nTrials, choiceOptions, confDispDuringChoice,...
        taskTimes,...
        results_folder, file_nm);
    
    %% add fixation cross to terminate the acquisition (to avoid weird fMRI behavior for last trial)
    Screen('FillRect',window, stim.cross.colour, stim.cross.verticalLine); % vertical line
    Screen('FillRect',window, stim.cross.colour, stim.cross.horizontalLine); % horizontal line
    [~,onsets.finalCross] = Screen('Flip',window); % display the cross on screen
    WaitSecs(taskTimes.finalCross);
    
    %% display feedback to prepare for last calibration and to inform for wait of end of fMRI acquisition
    DrawFormattedText(window, stim.endfMRIMessage.text,...
        stim.endfMRIMessage.x, stim.endfMRIMessage.y,...
        white, scr.wrapat);
    [~,onsets.endSessionFbk] = Screen(window,'Flip');
    WaitSecs(t_endfMRI);
    
    % indicate to the experimenter when to stop the fMRI acquisition
    disp('You can stop the fMRI acquisition now. When and only when fMRI acquisition has been stopped, please press space.');
    [~, ~, keyCode] = KbCheck();
    while(keyCode(key.space) ~= 1)
        % wait until the key has been pressed
        [~, ~, keyCode] = KbCheck();
    end
    disp('OK - space was pressed');
    
    %% get all TTL from the task
    if IRM == 1
        [TTL, keyLeft, keyRight,...
            keyLeftUnsure, keyLeftSure, keyRightUnsure, keyRightSure] = keyboard_check_end(TTL, key);
        % key storage of when left/right key have been pressed
        all.keys.keyLeft    = keyLeft;
        all.keys.keyRight   = keyRight;
        if key.n_buttonsChoice == 4
            all.keys.keyLeftUnsure  = keyLeftUnsure;
            all.keys.keyLeftSure    = keyLeftSure;
            all.keys.keyRightUnsure = keyRightUnsure;
            all.keys.keyRightSure   = keyRightSure;
        end
        
        % store T0 and TTL timings in onsets structure
        onsets.T0 = T0;
        onsets.TTL = TTL;
    end
    
    %% first save before last MVC
    save([results_folder, file_nm,'_messyAllStuff.mat']);
    
    %% Measure maximum power again at the end of each scan
    % add instructions
    DrawFormattedText(window, stim.postTaskMVCmeasurement.text,...
        stim.postTaskMVCmeasurement.x, stim.postTaskMVCmeasurement.y,...
        stim.postTaskMVCmeasurement.colour, scr.wrapat);
    Screen(window,'Flip');
    
    % MVC maximum
    switch effort_type
        case 'physical'
            [last_MVC, onsets_last_MVC] = physical_effort_MVC(scr, stim, dq, n_MaxPerfTrials, calibTimes);
        case 'mental'
            % extract numbers to use for each calibration trial
            %             [numberVector_endCalib] = mental_numbers(nFinalTrial);
            [numberVector_endCalib] = mental_calibNumberVector(n_MaxPerfTrials, n_calibMax);
            
            % use fixed number of correct answers to provide for each effort
            % level
            % repeat calibration until the subject performance is better
            % than the requested time threshold
            [t_min_finalMaxPerf, finalMaxPerf_SessionSummary, finalMaxPerf_calibSuccess] = mental_calibTime(scr, stim, key,...
                numberVector_endCalib, mentalE_prm_calib, n_MaxPerfTrials, n_calibMax, calibTimes, calib_errorLimits_Em, langage);
            calibEndSessionSummary.calibEndSession.calibSummary = finalMaxPerf_SessionSummary;
            calibEndSessionSummary.calibEndSession.calibSuccess = finalMaxPerf_calibSuccess;
            calibEndSessionSummary.calibEndSession.t_mental_max_perTrial = t_min_finalMaxPerf;
    end
    
    %% display feedback for the current session
    DrawFormattedText(window,...
        ['Felicitations! Cette session est maintenant terminee.',...
        'Vous avez obtenu: ',num2str(perfSummary.totalGain(nTrials)),...
        ' chf au cours de cette session.'],...
        stim.endSessionMessage.x, stim.endSessionMessage.y,...
        white, scr.wrapat);
    [~,onsets.endSessionFbk] = Screen(window,'Flip');
    WaitSecs(t_endSession);
    
    %% Save Data
    if IRM == 1
        % store all relevant data in final output variable
        all.choice_opt = choiceOptions; % choice options
        % store max before and after the session
        switch effort_type
            case 'physical'
                % initial calibration
                all.start_maxPerf.MVC = initial_MVC;
                all.start_maxPerf.onsets = onsets_initial_MVC;
                % last calibration
                all.end_maxPerf.MVC = last_MVC;
                all.end_maxPerf.onsets = onsets_last_MVC;
            case 'mental'
                if session_nb == 0
                    all.calibSummary = calibSummary;
                else
                    all.start_maxPerf.t_min = t_min_initialMaxPerf;
                    all.start_maxPerf.sessionSummary = initialMaxPerfSessionSummary;
                    all.start_maxPerf.success = initialMaxPerfSuccess;
                end
                % last max perf
                all.end_maxPerf.t_min = t_min_finalMaxPerf;
                all.end_maxPerf.SessionSummary = finalMaxPerf_SessionSummary;
                all.end_maxPerf.calibSuccess = finalMaxPerf_calibSuccess;
                all.end_maxPerf.calibEndSessionSummary = calibEndSessionSummary;
        end
        switch effort_type
            case 'physical'
                all.physicalPerf = perfSummary;
            case 'mental'
                if IRM == 0
                    all.mentalE_prm_learning = mentalE_prm_calib;
                end
                all.mentalE_perf = perfSummary;
        end
        % store timings in all structure
        all.calibTimes = calibTimes;
        all.taskTimes = taskTimes;
        all.onsets = onsets;
        
        
        % Double save to finish: .mat and .csv format
        save([results_folder, file_nm,'.mat'],'-struct','all');
    end
    
end % session number (calibration vs actual task)

%% STOP physiological recording
disp('Please stop physiological recording and then press space.');
[~, ~, keyCode] = KbCheck();
while(keyCode(key.space) ~= 1)
    % wait until the key has been pressed
    [~, ~, keyCode] = KbCheck();
end
disp('OK - space was pressed, physio recording stopped');

%% Clear the PTB screen
sca;

%% save all variables
save([results_folder, file_nm,'_messyAllStuff.mat']);