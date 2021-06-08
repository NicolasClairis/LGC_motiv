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
% See also ScreenConfiguration.m, physical_effort_perf.m,
% mental_effort_perf.m

%% Clear the workspace and the screen, instrreset resets the udp channels
sca; % close all PTB screens
close all; % close all windows
clearvars; % clear variables from memory
instrreset; % Disconnect and delete all instrument objects
clc;

%% define if you are currently testing the script (1)
% (no need to have correct timings and everything in PTB)
% or if this is the actual experiment => use optimal timings of the
% computer (0)
testing_script = 1;

%% working directories
% main_folder = ['D:' filesep 'Matlab_codes' filesep]; % Arthur's laptop path
% main_folder = ['C:',filesep,'Users',filesep,'Loco',filesep,...
%     'Documents',filesep,'GitHub',filesep,...
%     'LGC_motiv', filesep]; % Nico's laptop Github path
cd ..
main_folder                 = [pwd filesep]; % you have to be sure that you are in the correct path when you launch the script
main_task_folder            = [main_folder, 'LGC_Motiv_task' filesep];
results_folder              = [main_folder, 'LGC_Motiv_results' filesep];
% BioPac_folder               = [main_folder, 'BioPac_functions' filesep];
pics_folder                 = [main_task_folder, 'Coin_PNG', filesep];
Matlab_DIY_functions_folder = [main_folder, 'Matlab_DIY_functions', filesep];

% add personal functions (needed for PTB opening at least)
addpath(genpath(main_task_folder));
% addpath(BioPac_folder);
addpath(Matlab_DIY_functions_folder);

% create results folder if no subject has been acquired yet
if ~exist(results_folder,'dir')
    mkdir(results_folder);
end

%% task type (physical/mental effort) + subject + session identification
[effort_type_letter,...
    sub_initials,...
    sub_nber_str,...
    session_nm] = deal('');

% % effort type?
% while ~ismember(effort_type_letter,{'p','m'})
%     effort_type_letter = input('effort type? For physical, press ''p'', for mental press ''m''.','s');
% end
% switch effort_type_letter
%     case 'p'
%         effort_type = 'physical';
%     case 'm'
effort_type = 'mental';
% end

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

% calibration performance file name
calibPerf_file_nm = [results_folder,effort_type,'_calibPerf_sub',subjectCodeName,'.mat'];

% ask session number (especially for fMRI mapping with behavioral data)
% n_sess_max = 5; % maximal number of sessions
% while isempty(session_nm) ||...
%         str2double(session_nm) < 1 ||...
%         str2double(session_nm) > n_sess_max
%     session_nm = input('Session number? (write 0 if training session)','s');
% end
session_nm = '00';

if ismember(session_nm,{'0','00'}) || ~exist(calibPerf_file_nm,'file')
    learning_done = false;
elseif str2double(session_nm) > 0 && exist(calibPerf_file_nm,'file')
    learning_done = true;
else
    error('problem with session_nm specification or file loading');
end

% file name
file_nm = [subjectCodeName,'_session',session_nm,'_',effort_type,'_task'];
% verify the files do not already exist
if exist([results_folder, file_nm,'.mat'],'file')
    error(['The file name ',file_nm,'.mat already exists.',...
        ' Please relaunch with a new file name or delete the previous data.']);
end

%% fMRI/behavioral version of the task?
IRM = 0;
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

% define relevant keys and dynamometer
switch effort_type
    case 'mental'
        key = relevant_key_definition(effort_type, IRM);
    case 'physical' % need dq output to record the handgrip data
        [key, dq] = relevant_key_definition(effort_type, IRM);
end

% learning
switch effort_type
    case 'mental'
        % define number of pairs to solve for each level of difficulty
        
        % for learning: define the number of correct answers to provide
        % for the learning to be considered ok
        % adapt number of required correct answers according to if learning was done or not
        if ~learning_done
            n_maxLearning.learning_withInstructions = 6;%20;
            n_maxLearning.learning_withoutInstructions = 6;%30;
        elseif learning_done % short retraining at the beginning of each block to remind the mapping
            n_maxLearning.learning_withInstructions = 3;
            n_maxLearning.learning_withoutInstructions = 3;
        end % learning + training session
    case 'physical'
        n_learningForceRepeats = 3; % number of learning repetitions for each level of force
end

% calibration
switch effort_type % in case you use different numbers for each effort type
    case 'mental'
        n_calibTrials = 5;
    case 'physical'
        n_calibTrials = 5;
end

% max performance (beginning and end of each MRI session)
n_MaxPerfTrials = 3;

% actual task
% n_R_levels = 4;
% n_E_levels = 4;
% n_trials = 88;
n_R_levels = 3;
n_E_levels = 3;
nTrials = 44;

% extract money amount corresponding to each reward level for the
% computation of the gains
R_money = R_amounts(n_R_levels);

% check trial number is ok based on the number of entered conditions
% you should have a pair number of trials so that you can define an equal
% amount of reward and punishment trials
if strcmp(punishment_yn,'yes') && mod(nTrials,2) ~= 0
    error(['you added punishments in the task, hence you need a pair number of total trials.',...
        ' Please fix this so that you can have the same number of trials in punishment and reward']);
end
%
% determine reward/punishment and effort level combinations for each trial
choiceOptions = choice_option_design(n_R_levels, n_E_levels, punishment_yn, nTrials);

switch effort_type
    case 'physical'
        F_threshold = 50; % force should be maintained above this threshold (expressed in % of MVC)
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
[stim] = stim_initialize(scr, n_E_levels);
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
[trainingTimes, calibTimes, learningTimes, taskTimes, mainTimes] = timings_definition(trainingConditions, n_R_levels, n_E_levels, nTrials, effort_type);
t_endSession = mainTimes.endSession;

%% learning (for mental effort)
% for physical effort, you need to perform the MVC calibration before the
% learning phase
switch effort_type
    case 'mental' % for mental effort, define all the number sequences in advance
        
        %% mental_learning: to do at the beginning of each visit (but
        % not of every session in the scanner)
        if IRM == 0
            mentalE_prm_learning_and_calib = mental_effort_parameters(i_sub);
            mentalE_prm_learning_and_calib.startAngle = 0; % for learning always start at zero
            % no time limit for each trial: as long as needed until learning is
            % ok
            learning_time_limit = false;
            
            % perform 2 learning sessions, one with instructions and then one without
            % (left/right) vs (odd/even) and (lower/higher than 5) - mapping indicated the first time)
            % need to remind the mapping the second time
            learning_cols = {'col1','col2','all'};
            n_learningColours = length(learning_cols);
            learning_instructions = {'fullInstructions','noInstructions'}; %,'partialInstructions'
            n_learningInstructions = length(learning_instructions);
            % extract numbers to use for each learning phase
            [numberVector_learning] = mental_numbers(n_learningColours*n_learningInstructions);
            jLearningSession = 0;
            for iCol = 1:n_learningColours
                curr_learning_col = learning_cols{iCol};
                for iLearning_Instructions = 1:n_learningInstructions
                    curr_learning_instructions = learning_instructions{iLearning_Instructions};
                    
                    jLearningSession = jLearningSession + 1;
                    learning_sess_nm = ['learning_session',num2str(jLearningSession)];
                    % display instructions for the current learning type
                    [onsets.endLearningInstructions.(learning_sess_nm).(curr_learning_col).(curr_learning_instructions)] = mental_learning(scr,...
                        curr_learning_col, curr_learning_instructions, mentalE_prm_learning_and_calib);
                    
                    % perform the learning
                    [mentalE_learningPerf.(learning_sess_nm).(curr_learning_col).(curr_learning_instructions)] = mental_effort_perf(scr, stim, key,...
                        numberVector_learning(jLearningSession,:),...
                        mentalE_prm_learning_and_calib, n_maxLearning.learning_withInstructions,...
                        curr_learning_col, curr_learning_instructions, learning_time_limit);
                end % learning instructions loop
            end % learning colour loop
        end % MRI filter
end

%% calibration (before the MRI)
if IRM == 0 && ~learning_done
    switch effort_type
        case 'mental'
            %% max performance measurement
            % extract numbers to use for each calibration trial
            [numberVector_calib] = mental_numbers(n_calibTrials);
            
            %% alternatively, use fixed number of correct answers to provide for each effort
            % level
            % repeat calibration until the subject performance is better
            % than the requested time threshold
            calibSuccess = false;
            calibSession = 0;
            while calibSuccess == false
                calibSession = calibSession + 1;
                [t_min_calib, calibSessionSummary, calibSuccess] = mental_calibTime(scr, stim, key,...
                    numberVector_calib, mentalE_prm_learning_and_calib, n_calibTrials, n_calibMax, calibTimes);
                calibSummary.(['calibSession_',num2str(calibSession)]).calibSummary = calibSessionSummary;
                calibSummary.(['calibSession_',num2str(calibSession)]).calibSuccess = calibSuccess;
                calibSummary.(['calibSession_',num2str(calibSession)]).t_mental_max_perTrial = t_min_calib;
            end
            % save calib performance
            save(calibPerf_file_nm,'t_min_calib');
            
            
        case 'physical'% for physical effort, ask the MVC
            % record and store global MVC
            [initial_MVC, onsets_initial_MVC] = physical_effort_MVC(scr, dq, n_MVC_repeat, calibTimes);
            MVC = nanmax(initial_MVC); % expressed in Voltage
            save(calibPerf_file_nm,'MVC');
    end
    
else % load max performance from calibration
    switch effort_type
        case 'mental'
            t_min_calib = getfield(load(calibPerf_file_nm,'t_min_calib'),'t_min_calib');
        case 'physical'
            MVC = getfield(load(calibPerf_file_nm,'MVC'),'MVC'); % expressed in Voltage
    end
end

% mental effort: increase time available based on calibration max
% performance
switch effort_type
    case 'mental'
        t_training_max_effort = t_min_calib*trainingTimes.t_min_scalingFactor; % allow more time then min performance
        trainingTimes.max_effort = t_training_max_effort;
        t_max_effort = t_min_calib*taskTimes.t_min_scalingFactor; % allow more time then min performance
        taskTimes.max_effort = t_max_effort;
end

%% learning (for physical effort task)
% (needs to be done after the calibration)
switch effort_type
    case 'physical'
        % learning should be short to avoid exhausting the participant.
        % just perform X times each level of effort so that they get an idea of what it implies
        [learningPerfSummary, learningOnsets] = physical_learning(scr, stim, dq, n_E_levels, Ep_time_levels,...
            F_threshold, F_tolerance, MVC,...
            n_learningForceRepeats, learningTimes);
end

%% max perf measurement before start of each session
% (not for training out of MRI)
if IRM ~= 1 || learning_done
    switch effort_type
        case 'mental'
            % perform max perf
            [numberVector_initialMaxPerf] = mental_numbers(n_MaxPerfTrials);
            [t_min_initialMaxPerf, initialMaxPerfSessionSummary, initialMaxPerfSuccess] = mental_calibTime(scr, stim, key,...
                numberVector_initialMaxPerf, mentalE_prm_learning_and_calib, n_MaxPerfTrials, n_calibMax, calibTimes);
        case 'physical'
            % take an initial MVC measurement (even if it has been done in a
            % previous session, will allow us to keep track of the force level
            % of our participants)
            [initial_MVC, onsets_initial_MVC] = physical_effort_MVC(scr, stim, dq, n_MaxPerfTrials, calibTimes);
    end
end

%% short training should be added there (as similar as possible to the main task
% but no actual incentives involved? ie specify to the participants
% this will not be rewarded) => this should be included at the
% beginning of the training session
if IRM == 0
    for iTrainingCondition = 1:n_trainingConditions
        trainingCond_tmp = trainingConditions{iTrainingCondition};
        
        % define parameters for the training
        % reward/punishment and effort levels
        [trainingChoiceOptions_tmp, n_trainingTrials_tmp] = training_options(trainingCond_tmp, n_R_levels, n_E_levels);
        
        % start with reward training alone
        switch effort_type
            case 'physical'
                Ep_or_Em_vars.MVC = MVC;
                Ep_or_Em_vars.dq = dq;
                Ep_or_Em_vars.Ep_time_levels = Ep_time_levels;
                Ep_or_Em_vars.F_threshold = F_threshold;
                Ep_or_Em_vars.F_tolerance = F_tolerance;
            case 'mental'
                Ep_or_Em_vars.i_sub = i_sub;
                Ep_or_Em_vars.n_to_reach = n_to_reach;
        end
        [onsets_training.(trainingCond_tmp)] = choice_and_perf_trainingInstructions(scr, trainingCond_tmp, trainingTimes.instructions);
        [trainingSummary.(trainingCond_tmp)] = choice_and_perf(scr, stim, key, effort_type, Ep_or_Em_vars, R_money,...
            trainingCond_tmp, n_trainingTrials_tmp, trainingChoiceOptions_tmp, trainingTimes,...
            results_folder, file_nm);
    end % learning condition loop
    
    DrawFormattedText(window,'Bravo! Votre entraînement est terminé.',...
        'center','center',scr.colours.white, scr.wrapat);
    [~,onsets.EndTrainingMsg] = Screen('Flip',window); % display the cross on screen
    WaitSecs(trainingTimes.trainingEnd);
    sca; % close PTB if learning out of the MRI
end % training before MRI

%% launch main task
if IRM == 1 % || IRM == 0 %for piloting
    %% instruction that main task will start soon
    DrawFormattedText(window,...
        'L''expérimentateur va bientôt démarrer la tâche.',...
        'center', yScreenCenter*(5/3), scr.colours.white, scr.wrapat);
    [~, onsets.taskWillStart] = Screen(window, 'Flip');
    disp('Please press space and then launch fMRI (Be careful to respect this order for the T0...');
    [~, ~, keyCode] = KbCheck();
    while(keyCode(key.space) ~= 1)
        % wait until the key has been pressed
        [~, ~, keyCode] = KbCheck();
    end
    
    %% start recording fMRI TTL and wait for a given amount of TTL before
    % starting the task in order to calibrate all timings on T0
    if IRM == 1
        dummy_scans = 4; % number of TTL to wait before starting the task
        [T0, TTL] = keyboard_check_start(dummy_scans, key.trigger_id, key);
    end % fMRI check

    %% perform choice and performance task
    [perfSummary] = choice_and_perf(scr, stim, key,...
        effort_type, Ep_or_Em_vars, R_money,...
        'mainTask', nTrials, choiceOptions, taskTimes,...
        results_folder, file_nm);
    
    %% add fixation cross to terminate the acquisition (to avoid weird fMRI behavior for last trial)
    Screen('FillRect',window,white, stim.cross.verticalLine); % vertical line
    Screen('FillRect',window,white, stim.cross.horizontalLine); % horizontal line
    [~,onsets.finalCross] = Screen('Flip',window); % display the cross on screen
    WaitSecs(taskTimes.finalCross);
    
    %% display feedback for the current session
    DrawFormattedText(window,...
        ['Félicitations! Cette session est maintenant terminée.',...
        'Vous avez obtenu: ',num2str(perfSummary.totalGain(nTrials)),' chf au cours de cette session.'],white);
    [~,onsets.endSessionFbk] = Screen(window,'Flip');
    WaitSecs(t_endSession);
        
    % indicate to the experimenter when to stop the fMRI acquisition
    disp('You can stop the fMRI acquisition now. When and only when it is done, please press space.');
    [~, ~, keyCode] = KbCheck();
    while(keyCode(key.space) ~= 1)
        % wait until the key has been pressed
        [~, ~, keyCode] = KbCheck();
    end
    
end % MRI only (not for training out of fMRI)

%% get all TTL from the task
if IRM == 1
    switch session_effort_type
        case 'physical'
            TTL = keyboard_check_end(TTL);
        case 'mental'
            [TTL, keyLeft, keyRight] = keyboard_check_end(TTL);
            % key storage of when left/right key have been pressed
            all.keys.keyLeft    = keyLeft;
            all.keys.keyRight   = keyRight;
    end
    % store T0 and TTL timings in onsets structure
    onsets.T0 = T0;
    onsets.TTL = TTL;
end

%% first save before last MVC
save([results_folder, file_nm,'_messyAllStuff.mat']);

%% Measure maximum power again at the end of each scan
if IRM == 1
    % add instructions
    DrawFormattedText(window,...
        ['Pour finir cette session, nous allons vous demander ',...
        'd''essayer à nouveau de battre votre record.'],...
        'center', yScreenCenter*(5/3), scr.colours.white, scr.wrapat);
    Screen(window,'Flip');
    % MVC maximum
    nFinalTrial = 1;
    switch effort_type
        case 'physical'
            [MVC_last, onsets_MVC_last] = physical_effort_MVC(scr, stim, dq, nFinalTrial, calibTimes);
        case 'mental'
            % extract numbers to use for each calibration trial
            [numberVector_endCalib] = mental_numbers(nFinalTrial);
            
            % use fixed number of correct answers to provide for each effort
            % level
            % repeat calibration until the subject performance is better
            % than the requested time threshold
            [t_min_finalMaxPerf, finalMaxPerf_SessionSummary, finalMaxPerf_calibSuccess] = mental_calibTime(scr, stim, key,...
                numberVector_endCalib, mentalE_prm_learning_and_calib, nFinalTrial, n_calibMax, calibTimes);
            calibEndSessionSummary.calibEndSession.calibSummary = finalMaxPerf_SessionSummary;
            calibEndSessionSummary.calibEndSession.calibSuccess = finalMaxPerf_calibSuccess;
            calibEndSessionSummary.calibEndSession.t_mental_max_perTrial = t_min_finalMaxPerf;
    end
end

%% Clear the PTB screen
sca;

%% Save Data
if IRM == 1
    % store all relevant data in final output variable
    all.choice_opt = choiceOptions; % choice options
    % store max before and after the session
    switch effort_type
        case 'physical'
            all.start_maxPerf.MVC = initial_MVC;
            all.start_maxPerf.onsets = onsets_initial_MVC;
        case 'mental'
            if ~learning_done
                all.calibSummary = calibSummary;
            elseif learning_done
                all.start_maxPerf.t_min = t_min_initialMaxPerf;
                all.start_maxPerf.sessionSummary = initialMaxPerfSessionSummary;
                all.start_maxPerf.success = initialMaxPerfSuccess;
            end
            % last max perf
            all.end_maxPerf.t_min = t_min_finalMaxPerf;
            all.end_maxPerf.SessionSummary = finalMaxPerf_SessionSummary;
            all.end_maxPerf.calibSuccess = finalMaxPerf_calibSuccess;
    end
    % store performance
    all.calibSummary = calibSummary;
    switch effort_type
        case 'physical'
            all.physicalPerf = perfSummary;
            all.learning.perf = learningPerfSummary;
            all.learning.onsets = learningOnsets;
        case 'mental'
            if IRM == 0
                all.learningPerf = mentalE_learningPerf;
                all.mentalE_prm_learning = mentalE_prm_learning_and_calib;
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

save([results_folder, file_nm,'_messyAllStuff.mat']);