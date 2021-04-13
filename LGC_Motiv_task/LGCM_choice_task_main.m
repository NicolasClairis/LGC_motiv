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
% main_folder = ['C:',filesep,'Users',filesep,'Loco',filesep,...
%     'Documents',filesep,'GitHub',filesep,...
%     'LGC_motiv', filesep]; % Nico's laptop Github path
cd ..
main_folder                 = [pwd filesep]; % you have to be sure that you are in the correct path when you launch the script
main_task_folder            = [main_folder, 'LGC_Motiv_task' filesep];
savePath                    = [main_folder, 'LGC_Motiv_results' filesep];
% BioPac_folder               = [main_folder, 'BioPac_functions' filesep];
pics_folder                 = [main_task_folder, 'Coin_PNG', filesep];
Matlab_DIY_functions_folder = [main_folder, 'Matlab_DIY_functions', filesep];

% add personal functions (needed for PTB opening at least)
addpath(genpath(main_task_folder));
% addpath(BioPac_folder);
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

%% if physical effort task => requires opening the connection to the BioPac
switch effort_type
    case 'physical'
        % open communication with National Instruments card (which should
        % be connected to the Biopac module, itself connected to the
        % handgrip. Be careful to set everything properly in the channel 1
        % as input (grip) and output (to NI card)
        
        % you can use the following functions to check the name of the
        % device
        % d = daqlist;
        % d{1, "DeviceInfo"}
       
        % initialize NI input
        dq = daq("ni");
%         dq.Rate = 500; % define the acquisition rate
        NI_module_nm = "cDAQ1Mod1";
        NI_channel_output = "ai0";
        addinput(dq, NI_module_nm, NI_channel_output,"Voltage");
        % will need data = read(dq) function only to read the signal
end

%% initialize screen
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = LGCM_ScreenConfiguration(IRM, testing_script);
white = scr.colours.white;
black = scr.colours.black;

%% define bar size for the waiting time
barTimeWaitRect = [xScreenCenter*(1/2),...
    yScreenCenter*(3/4),...
    xScreenCenter*(3/2),...
    yScreenCenter];

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
n_trials = 88;
n_R_levels = 4;
n_E_levels = 4;
switch effort_type
    case 'physical'
        F_threshold = 50; % force should be maintained above this threshold
        F_tolerance = 5; % tolerance allowed around the threshold
        % need to define timings for each level of force
        switch n_E_levels
            case 3
                E_time_levels.level_1 = 1;
                E_time_levels.level_2 = 2;
                E_time_levels.level_3 = 3;
            case 4
                E_time_levels.level_1 = 0.5;
                E_time_levels.level_2 = 1.5;
                E_time_levels.level_3 = 2.5;
                E_time_levels.level_4 = 3.5;
        end
end

% stimulus related variables for the display
[stim] = LGCM_stim_initialize(scr, n_R_levels, n_E_levels, pics_folder);

% determine reward and punishment trials
choice_opt = LGCM_choice_option_design(n_R_levels, n_E_levels, punishment_yn, n_trials);
%%
warning('need to adapt choice options design to fit with different number of conditions (3 instead of 4) and trials');
%%

switch punishment_yn
    case 'yes'
        R_or_P = choice_opt.R_or_P;
    case 'no' % by default all trials are rewarding if no punishment
        R_or_P = repmat({'R'},1,n_trials);
    otherwise
        error('punishment_yn has not been initiliazed.');
end

% define the number of correct answers to provide for the learning to be
% considered ok
% adapt number of required correct answers according to if learning was done or not
if ~learning_done
    n_max.learning_withInstructions = 6;%20;
    n_max.learning_withoutInstructions = 6;%30;
    % calibration max
    n_calibMax = 4;
elseif learning_done % short retraining at the beginning of each block
    n_max.learning_withInstructions = 6;%10;
    n_max.learning_withoutInstructions = 6;%10;
end % learning + training session

%% timings
% calibration timings
calibTimes.instructions = 5;
switch effort_type % in case you use different numbers for each effort type
    case 'mental'
        calibTimes.effort_max = 20; % maximal time to perform the task (calibrated time should be shorter)
    case 'physical'
        calibTimes.effort_max = 5;% time to perform the task
end
calibTimes.fbk = 2;

% main task timings
t_cross = 0.5;
t_choice = 5;
t_dispChoice = 1.5;
switch effort_type
    case 'physical' % in case you use different numbers for each effort type
        t_max_effort = 5; % time to perform the task
        taskTimes.max_effort = t_max_effort;
    case 'mental'
        t_min_scalingFactor = 150/100; % multiply calibrated minimal time by this value
end
t_fbk = 1; % feedback display
taskTimes.cross = t_cross;
taskTimes.choice = t_choice;
taskTimes.dispChoice = t_dispChoice;
taskTimes.feedback = t_fbk;

%% learning + calibration
switch effort_type
    case 'mental' % for mental effort, define all the number sequences in advance
        mentalE_prm = LGCM_mental_effort_parameters(i_sub);
        mentalE_prm.startAngle = 0; % for learning always start at zero
        
        %% LGCM_mental_learning
        % no time limit for each trial: as long as needed until learning is
        % ok
        learning_time_limit = false;
        
        % perform 2 learning sessions, one with instructions and then one without
        % (left/right) vs (odd/even) and (lower/higher than 5) - mapping indicated the first time)
        % need to remind the mapping the second time
        learning_cols = {'col1','col2','all'}; n_learningColours = length(learning_cols);
        learning_instructions = {'fullInstructions','partialInstructions','noInstructions'};
        n_learningInstructions = length(learning_instructions);
        % extract numbers to use for each learning phase
        [numberVector_learning] = LGCM_mental_numbers(n_learningColours*n_learningInstructions);
        jLearningSession = 0;
        for iCol = 1:n_learningColours
            curr_learning_col = learning_cols{iCol};
            for iLearning_Instructions = 1:n_learningInstructions
                curr_learning_instructions = learning_instructions{iLearning_Instructions};
                
                jLearningSession = jLearningSession + 1;
                learning_sess_nm = ['learning_session',num2str(jLearningSession)];
                % display instructions for the current learning type
                [onsets.endLearningInstructions.(learning_sess_nm).(curr_learning_col).(curr_learning_instructions)] = LGCM_mental_learning(scr,...
                    curr_learning_col, curr_learning_instructions, mentalE_prm);
                
                % perform the learning
                [mentalE_learningPerf.(learning_sess_nm).(curr_learning_col).(curr_learning_instructions)] = LGCM_mental_effort_perf(scr, stim, key,...
                    numberVector_learning(jLearningSession,:),...
                    mentalE_prm, n_max.learning_withInstructions,...
                    curr_learning_col, curr_learning_instructions, learning_time_limit);
            end % learning instructions loop
        end % learning colour loop
        
        %% max performance measurement
        if ~learning_done
            % extract numbers to use for each calibration trial
            [numberVector_calib] = LGCM_mental_numbers(n_calibTrials);
            
            %% in case of using a fixed time and calibrating the maximal number of correct answers
%             n_mental_max_perTrial = LGCM_mental_calibNumbers(scr, stim, key,...
%                 numberVector_calib, mentalE_prm, n_calibTrials,...
%                 n_calibMax, calibTimes);

            %% alternatively, use fixed number of correct answers to provide for each effort
            % level, but calibrate the maximal duration required for the
            % top effort
            n_to_reach = LGCM_mental_N_answersPerLevel(n_E_levels);
            n_calibMax = n_to_reach.(['E_level_',num2str(n_E_levels)]);
            % repeat calibration until the subject performance is better
            % than the requested time threshold
            calibSuccess = false;
            calibSession = 0;
            while calibSuccess == false
                calibSession = calibSession + 1;
                [t_min_calib, calibSessionSummary, calibSuccess] = LGCM_mental_calibTime(scr, stim, key,...
                    numberVector_calib, mentalE_prm, n_calibTrials, n_calibMax, calibTimes);
                t_max_effort = t_min_calib*t_min_scalingFactor; % allow more time then min performance
                calibSummary.(['calibSession_',num2str(calibSession)]).calibSummary = calibSessionSummary;
                calibSummary.(['calibSession_',num2str(calibSession)]).calibSuccess = calibSuccess;
                calibSummary.(['calibSession_',num2str(calibSession)]).t_mental_max_perTrial = t_min_calib;
            end
        elseif learning_done
            % indicate max performance
            %  n_mental_max_perTrial = input('please write manually n_max of pairs done during training');
            t_max_effort = input('please write manually t_max_to_use based on calibration');
            warning('ideally replace manual input by loading the value with matlab');
        end % learning + training session
        taskTimes.max_effort = t_max_effort;
        
        % in case you calibrate number of correct answers during the
        % training
%         % define levels of difficulty based on the calibration
%         for iE_level = 1:n_E_levels
%             n_to_reach.(['E_level_',num2str(iE_level)]) = round((iE_level/(n_E_levels + 1))*n_mental_max_perTrial);
%             warning('you should check that all levels are effort can be differentiated or you have a big problem');
%         end
        
        %% initializing variables for the main task
        % randomize the order of the numbers appearing on screen
        mental_nbers_per_trial = LGCM_mental_numbers(n_trials);
        
        % randomize the type of the first trial (odd/even or higher/lower
        % than 5)
        mental_taskType_trialStart = LGCM_mental_task_start(n_trials);
        
        %% display message
        DrawFormattedText(window,...
            'L''expérimentateur va bientôt démarrer la tâche.',...
            'center', yScreenCenter*(5/3), scr.colours.black, scr.wrapat);
        [~, time_fbkPress] = Screen(window, 'Flip');
        KbWait;
        
    case 'physical'% for physical effort, ask the MVC
        
        % take an initial MVC measurement
        
        % store global MVC
        if exist([],'file') % retrieve past MVC value if has been extracted already
            
            
        else % perform MVC measurement otherwise
            
        end
        % [MVC_initial, onsets_MVC_initial] = LGCM_MVC_measurement(scr, stim, session_effort_type, n_MVC_repeat);
        warning('MVC measurement function needs to be updated given last task changes');
        
end

%% initialize onsets for main task
[onsets.cross,...
    onsets.dispChoiceOptions,...
    onsets.choice,...
    onsets.keyReleaseMessage,...
    onsets.cross_after_buttonRelease,...
    onsets.fbk, onsets.fbk_win, onsets.fbk_loss,...
    onsets.fbk_fail,...
    onsets.timeBarWait] = deal(NaN(1,n_trials));
% variables during effort period should record the data for each trial
[onsets.effortPeriod,...
    mentalE_perf] = deal(cell(1, n_trials));

%% start recording fMRI TTL
if IRM == 1
    dummy_scan = 4; % number of TTL to wait before starting the task
    trigger_id = 84; % trigger value corresponding to the TTL code
    [T0, TTL] = LGCM_keyboard_check_start(dummy_scans, trigger_id, key);
end % fMRI check

%% launch main task
choice = zeros(1,n_trials);
[R_chosen, E_chosen,...
    effortTime] = deal(NaN(1, n_trials));
was_a_key_pressed_bf_trial = NaN(1,n_trials);

time_limit = true; % time limit to reach level of force required
for iTrial = 1:n_trials
    
    %% fixation cross period
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
    
    %% choice period
    [choice(iTrial),...
        onsets.dispChoiceOptions(iTrial),...
        onsets.choice(iTrial),...
        stoptask] = LGCM_choice_period(scr, stim,...
        choice_opt, R_or_P{iTrial}, iTrial, t_choice, key);
    
    %% check if escape was pressed => stop everything if so
    if stoptask == 1
        % save all the data in case you still want to analyze it
        save([savePath, file_nm,'_earlyEnd_tmp.mat'],'-struct','all');
        break;
    end

    %% chosen option display period
    [time_dispChoice,...
        R_chosen(iTrial),...
        E_chosen(iTrial)] = LGCM_choice_task_dispChosen(scr, stim, choice_opt,...
        choice(iTrial), R_or_P{iTrial},...
        iTrial);
    onsets.dispChoice(iTrial) = time_dispChoice;
    WaitSecs(t_dispChoice);
    
    %% Effort period 
    % perform the effort if a choice was made, otherwise punish the
    % subject for not answering to the choice
    if choice(iTrial) == 0 % no choice was made => failure
        trial_was_successfull_tmp = 0;
        
    elseif ismember(choice(iTrial), [-1,1]) % choice done => perform the corresponding effort
        
        
        %% perform the effort
        tic;
        switch effort_type
            case 'physical'
                [physicalE_perf{iTrial},...
                    trial_was_successfull_tmp,...
                    onsets.effortPeriod{iTrial}] = LGCM_physical_effort(scr, stim, key,...
                    time_limit, t_max_effort);
                
                [physicalE_perf{iTrial}, trial_was_successfull_tmp, onsets.effortPeriod{iTrial}] = LGCM_physical_effort_perf(scr, stim,...
                    E_chosen,...
                    E_time_levels,...
                    F_threshold, F_tolerance,...
                    time_limit, t_max_effort);
            case 'mental'
                mentalE_prm.startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen(iTrial))]); % adapt start angle to current level of difficulty
                n_max_to_reach_tmp = n_to_reach.(['E_level_',num2str(E_chosen(iTrial))]);
                [mentalE_perf{iTrial},...
                    trial_was_successfull_tmp,...
                    onsets.effortPeriod{iTrial}] = LGCM_mental_effort_perf(scr, stim, key,...
                    mental_nbers_per_trial(iTrial,:),...
                    mentalE_prm, n_max_to_reach_tmp,...
                    'all', 'noInstructions', time_limit, t_max_effort);
        end % effort type loop
        effortTime(iTrial) = toc;
    end % choice made or not?
    
    %% Feedback period
    switch trial_was_successfull_tmp
        case 0 % trial failed
            % display message
%             DrawFormattedText(window,'Too slow!', xScreenCenter, (1/5)*yScreenCenter);
            DrawFormattedText(window,'Trop lent!',...
                'center', (1/5)*yScreenCenter,...
                white);
            % display money loss for failing
            
            [~,onsets.fbk_fail(iTrial)] = Screen(window,'Flip');
            
        case 1 % trial is a success
            % display feedback
            switch R_or_P{iTrial}
                case 'R'
%                     DrawFormattedText(window,'You won', xScreenCenter, (1/5)*yScreenCenter);
                    DrawFormattedText(window,'Vous avez obtenu',...
                        'center', (1/5)*yScreenCenter,...
                        white);
                case 'P'
%                     DrawFormattedText(window,'You lost', xScreenCenter, (1/5)*yScreenCenter);
                    DrawFormattedText(window,'Vous avez perdu',...
                        'center', (1/5)*yScreenCenter,...
                        white);
            end
            
            % display monetary incentive
            Screen('DrawTexture', window,...
            stim.reward.texture.(['reward_',num2str(R_chosen(iTrial))]),...
            [],...
            stim.chosenOption.reward.(['reward_',num2str(R_chosen(iTrial))]));
            
            [~,onsets.fbk(iTrial)] = Screen(window,'Flip');
            switch R_or_P{iTrial}
                case 'R'
                    onsets.fbk_win(iTrial) = onsets.fbk(iTrial);
                case 'P'
                    onsets.fbk_loss(iTrial) = onsets.fbk(iTrial);
            end
    end
    WaitSecs(t_fbk);
    
    %% Time waiting period
    % this period allows to de-confound effort and delay, ie even if lower
    % effort has been selected and performed quicker, a short waiting time
    % will force to wait
    if effortTime(iTrial) < t_max_effort
        tic;
        onsets.timeBarWait(iTrial) = GetSecs; % record start of time bar
        
        % show a dynamic waiting bar until the timing ends
        while toc < (t_max_effort - effortTime(iTrial))
            % update bar with time remaining
            percTimeAchieved = (toc + effortTime(iTrial))./t_max_effort;
            barTimeWaitRect_bis = barTimeWaitRect;
            % start on the right corner of the bar + percentage already
            % achieved and move to the left
            barTimeWaitRect_bis(3) = barTimeWaitRect(3) - percTimeAchieved*(barTimeWaitRect(3) - barTimeWaitRect(1));
            %
            DrawFormattedText(window,'Temps restant','center',yScreenCenter*(1/2));
            % draw one global fixed rectangle showing the total duration
            Screen('FrameRect',window, black, barTimeWaitRect);
            % draw one second rectangle updating dynamically showing the
            % time remaining
            Screen('FillRect',window, black, barTimeWaitRect_bis);
            
            Screen(window,'Flip');
        end % display until time catches up with maximum effort time
    end % if all time taken, no need for time penalty
    
    %% display number of trials done for the experimenter
    disp(['Trial ',num2str(iTrial),'/',num2str(n_trials),' done']);
end % trial loop

%% Measure maximum power again at the end
% [MVC_last, onsets_MVC_last] = LGCM_MVC_measurement(scr, session_effort_type, speed, stim, n_MVC_repeat);

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
all.choice_opt = choice_opt; % choice options
% choice made
all.choice = choice;
all.R_chosen = R_chosen;
all.E_chosen = E_chosen;
% effort informations
all.effortTime = effortTime;
if strcmp(effort_type,'physical')
    all.MVC(1) = MVC_initial.MVC; % store initial MVC
    all.MVC(2) = MVC_last.MVC; % store final MVC
    all.MVC_allValues_perTest(:,1) = MVC_initial.MVC_perCalibSession;
    all.MVC_allValues_perTest(:,2) = MVC_last.MVC_perCalibSession;
    all.MVC_initial = MVC_initial;
    all.MVC_last = MVC_last;
    onsets.initial_MVC = onsets_MVC_initial;
    onsets.final_MVC = onsets_MVC_final;
end
% store performance
all.calibSummary = calibSummary;
if strcmp(effort_type,'mental')
    all.learningPerf = mentalE_learningPerf;
    all.perf = mentalE_perf;
end
% store timings in all structure
all.calibTimes = calibTimes;
all.taskTimes = taskTimes;
all.onsets = onsets;


% Double save to finish: .mat and .csv format
save([savePath, file_nm,'.mat'],'-struct','all');

save([savePath, file_nm,'_messyAllStuff.mat']);