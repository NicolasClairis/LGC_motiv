% script for behavioral pilots

%% clean workspace before starting
sca;
clearvars;
close all;
instrreset; % Disconnect and delete all instrument objects
clc;
%% working directories
% launch within the folder where scripts are stored or will not work
cd ..
main_folder                 = [pwd filesep]; % you have to be sure that you are in the correct path when you launch the script
main_task_folder            = [main_folder, 'LGC_Motiv_task' filesep];
results_folder              = [main_folder, 'LGC_Motiv_results' filesep];
% BioPac_folder               = [main_folder, 'BioPac_functions' filesep];
% pics_folder                 = [main_task_folder, 'Coin_PNG', filesep];
Matlab_DIY_functions_folder = [main_folder, 'Matlab_DIY_functions', filesep];

% add personal functions (needed for PTB opening at least)
addpath(genpath(main_task_folder));
% addpath(BioPac_folder);
addpath(Matlab_DIY_functions_folder);

% create results folder if no -subject has been acquired yet
if ~exist(results_folder,'dir')
    mkdir(results_folder);
end
% go back to folder with scripts
cd(main_task_folder);

%% Define subject ID

% Insert the initials, the number of the participants
iSubject = [];
langue = 'f';
IRM = 0;
while isempty(iSubject) || length(iSubject) ~= 3
    % repeat until all questions are answered
    info = inputdlg({'Subject CID (XXX)','p/m'});
    [iSubject,p_or_m] = info{[1,2]};
    %     warning('when real experiment starts, remember to block in french and IRM = 1');
end
if ischar(IRM)
    IRM = str2double(IRM);
end
% Create subjectCodeName which is used as a file saving name
subjectCodeName = strcat('CID',iSubject);

if strcmp(p_or_m,'m') == 0 && strcmp(p_or_m,'p') == 0
    error('this letter has no definition');
end

file_nm_training_Em = ['training_data_Em_CID',num2str(iSubject)];
file_nm_training_Ep = ['training_data_Ep_CID',num2str(iSubject)];
file_nm = ['training_data_CID',num2str(iSubject)];
file_nm_IP = ['delta_IP_CID',num2str(iSubject)];
% convert subject CID into number (only if used to perform actual task)
if ischar(iSubject)
    iSubject = str2double(iSubject);
end
%% general parameters
% define subparts of the task to perform (on/off)
taskToPerform.physical.calib = 'off';
taskToPerform.physical.learning = 'off';
taskToPerform.physical.training = 'off';
switch IRM
    case 0
        taskToPerform.physical.task = 'off';
    case 1 % task will be done in the scanner after the training
        taskToPerform.physical.task = 'off';
end
taskToPerform.mental.learning = 'on';
taskToPerform.mental.calib = 'on';
taskToPerform.mental.training = 'on';
switch IRM
    case 0
        taskToPerform.mental.task = 'on';
    case 1 % task will be done in the scanner after the training
        taskToPerform.mental.task = 'on';
end
switch langue
    case 'f'
        langage = 'fr';
    case 'e'
        langage = 'engl'; % 'fr'/'engl' french or english?
    otherwise
        error('langage not recognised');
end
% initialize screen
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(IRM, 1);
white = scr.colours.white;
black = scr.colours.black;

% include punishment condition?
punishment_yn = 'yes'; % include punishment trials?

% number of reward and effort conditions
n_R_levels = 4;
n_E_levels = 4;

% prepare multiple versions of efforts for indifference point
E_right = [3];
E_left = [1 1];

%Number of repeats of the whole code
nbRepeat = 1;
nbEffortLvl = 1;

% Baseline Reward (CHF)
baselineR = 0.5;
baselineP = 0.5;

% Total amount of money to be given
totalGain = 0;
finalGain = 0;
% define number of trials per staircase procedure
n_trialsPerSession = 5;

% mapping between reward levels and actual monetary amounts
R_money = R_amounts(n_R_levels, punishment_yn);

% initialize visual stimuli to use in the experiment
[stim] = stim_initialize(scr, n_E_levels, langage, R_money);

% define number of training conditions
switch punishment_yn
    case 'yes'
        switch IRM
            case 0
                trainingConditions = {'RP'};
            case 1
                trainingConditions = {'RP'};
        end
    case 'no'
        trainingConditions = {'R'};
end
n_trainingConditions = length(trainingConditions);

% load timings for each phase of the experiment
[trainingTimes_Em, calibTimes_Em, learningTimes_Em, taskTimes_Em] = timings_definition(trainingConditions, n_R_levels, n_E_levels, n_trialsPerSession, 'mental');
[trainingTimes_Ep, calibTimes_Ep, learningTimes_Ep, taskTimes_Ep, mainTimes] = timings_definition(trainingConditions, n_R_levels, n_E_levels, n_trialsPerSession, 'physical');

% number of times we apply a staircase procedure. (2 mental (R or P) 2 physical (R or P))
if strcmp(taskToPerform.physical.task,'on') && strcmp(taskToPerform.mental.task,'on')
    n_sessions = 4;
else
    n_sessions = 2;
end
% number of buttons to answer
switch IRM
    case 0
        n_buttonsChoice = 2;
    case 1 % test buttons
        n_buttonsChoice = 4;
end

% mental calibration error management: no fail after 3 errors nor
% mapping display
% calibration happens in 2 cases:
% 1) if it hasn't been made yet (during training)
% 2) before and after each session of the main task
% therefore you need to set up the parameters for both cases, not just for
% the main calibration
if strcmp(taskToPerform.mental.calib,'on') || strcmp(taskToPerform.mental.task,'on')
    calib_errorLimits_Em.useOfErrorMapping = false;
    calib_errorLimits_Em.useOfErrorThreshold = true;
    calib_errorLimits_Em.errorThreshold = 2;
end
% time for end of session
t_endSession = trainingTimes_Ep.endSession;

%% physical parameters
if strcmp(taskToPerform.physical.calib,'on') ||...
        strcmp(taskToPerform.physical.learning,'on') ||...
        strcmp(taskToPerform.physical.training,'on') ||...
        strcmp(taskToPerform.physical.task,'on')
    % define relevant keys and dynamometer module
    [key_Ep, dq] = relevant_key_definition('physical', IRM, n_buttonsChoice);
    % define conditions
    n_MVC_repeat = 3;
    n_learningForceRepeats = 1; % number of learning repetitions for each level of difficulty (= each level of force)
    F_threshold = 55; % force should be maintained above this threshold (expressed in % of MVC)
    F_tolerance = 2.5; % tolerance allowed around the threshold (expressed in % of MVC)
    % need to define timings for each level of force
    [Ep_time_levels] = physical_effortLevels(n_E_levels);
end

%% mental parameters
if strcmp(taskToPerform.mental.calib,'on') ||...
        strcmp(taskToPerform.mental.learning,'on') ||...
        strcmp(taskToPerform.mental.training,'on') ||...
        strcmp(taskToPerform.mental.task,'on')
    % define relevant keys and dynamometer module
    key_Em = relevant_key_definition('mental', IRM, n_buttonsChoice);
    
    % calibration: calibrate the maximal duration required for the top effort
    n_calibTrials_Em = 3;
    
    % learning
    % perform 2 learning sessions, one with instructions and then one without
    % (left/right) and (lower/higher than 5) - mapping indicated the first time)
    % need to remind the mapping the second time
    learning_instructions = {'fullInstructions','noInstructions'}; %,'partialInstructions'
    n_learningInstructions = length(learning_instructions);
    % initial learning: careful to enter a pair number here
    n_maxLearning.learning_withInstructions = 15;
    n_maxLearning.learning_withoutInstructions = 15;
    
end

<<<<<<< Updated upstream
%% calibration mental
if strcmp(taskToPerform.mental.calib,'on')
    mentalE_prm_learning_and_calib = mental_effort_parameters();
    mentalE_prm_learning_and_calib.startAngle = 0; % for learning always start at zero
    n_calibMax = mentalE_prm_learning_and_calib.n_maxToReachCalib;
    % extract numbers to use for each calibration trial
%     [numberVector_calib] = mental_numbers(n_calibTrials_Em);
    [numberVector_calib] = mental_calibNumberVector(n_calibTrials_Em, n_calibMax);
    % error handling for extended learning
    extendedLearning_errorLimits.useOfErrorThreshold = false;
    extendedLearning_errorLimits.useOfErrorMapping = false;
    % alternatively, use fixed number of correct answers to provide for each effort
    % level
    % repeat calibration until the subject performance is better
    % than the requested time threshold
    calibSuccess = false;
    iCalibSession = 0;
    while calibSuccess == false
        iCalibSession = iCalibSession + 1;
        [n_mental_max_perTrial, calib_summary] = mental_calibNumbers(scr, stim, key_Em,...
            numberVector_calib, mentalE_prm_learning_and_calib, n_calibTrials_Em, calibTimes_Em, langage);
        calibSummary.(['calibSession_',num2str(iCalibSession)]).calibSummary = calib_summary;
        calibSummary.(['calibSession_',num2str(iCalibSession)]).n_mental_max_perTrial = n_mental_max_perTrial;
=======
%% run the code twice, each time for p or m conditions
for i_pm = 1:2
    switch p_or_m
        case 'p'
            %% physical preparation
            %% physical MVC
            if strcmp(taskToPerform.physical.calib,'on')
                waitSpace(langage, window, yScreenCenter, scr, key_Ep);
                
                [initial_MVC, onsets_initial_MVC] = physical_effort_MVC(scr, stim, dq, n_MVC_repeat, calibTimes_Ep);
                MVC = max(initial_MVC.MVC); % expressed in Voltage
            end
            
            % learning physical
            if strcmp(taskToPerform.physical.learning,'on')
                showTitlesInstruction(scr,stim,'learning',false);
                waitSpace(langage, window, yScreenCenter, scr, key_Ep);
                
                [learningPerfSummary_Ep, learningOnsets_Ep] = physical_learning(scr, stim, dq, n_E_levels, Ep_time_levels,...
                    F_threshold, F_tolerance, MVC,...
                    n_learningForceRepeats, learningTimes_Ep);
            end
            
            % training physical
            if strcmp(taskToPerform.physical.training,'on')
                showTitlesInstruction(scr,stim,'training',false);
                waitSpace(langage, window, yScreenCenter, scr, key_Ep);
                % initialize training parameters
                Ep_vars_training.MVC = MVC;
                Ep_vars_training.dq = dq;
                Ep_vars_training.Ep_time_levels = Ep_time_levels;
                Ep_vars_training.F_threshold = F_threshold;
                Ep_vars_training.F_tolerance = F_tolerance;
                Ep_vars_training.timeRemainingEndTrial_ONOFF = 0;
                
                for iTrainingCondition = 1:n_trainingConditions
                    trainingCond = trainingConditions{iTrainingCondition};
                    % display confidence mapping only for first training sessions and
                    % only for fMRI experiment
                    if strcmp(trainingCond,'RP') || IRM == 0 || n_buttonsChoice == 2
                        confidenceChoiceDisplay = false;
                    elseif ~strcmp(trainingCond,'RP') && IRM == 1 && n_buttonsChoice == 4
                        confidenceChoiceDisplay = true;
                    end
                    
                    % de4fine parameters for the training
                    % reward/punishment and effort levels
                    [trainingChoiceOptions_Ep_tmp, n_trainingTrials_Ep_tmp] = training_options(trainingCond, n_R_levels, n_E_levels, R_money);
                    
                    % start with reward training alone
                    [onsets_Ep_training.(trainingCond)] = choice_and_perf_trainingInstructions(scr, stim, trainingCond, trainingTimes_Ep.instructions);
                    [trainingSummary_Ep.(trainingCond)] = choice_and_perf(scr, stim, key_Ep, 'physical', Ep_vars_training, R_money,...
                        trainingCond, n_trainingTrials_Ep_tmp, trainingChoiceOptions_Ep_tmp, confidenceChoiceDisplay,...
                        trainingTimes_Ep,...
                        results_folder, file_nm_training_Ep);
                end % learning condition loop
                
                DrawFormattedText(window, stim.training.Ep.endMsg.text,...
                    'center','center',stim.training.Ep.endMsg.colour, scr.wrapat);
                [~,onsets.EndTrainingMsg] = Screen('Flip',window); % display the cross on screen
                WaitSecs(trainingTimes_Ep.trainingEnd);
            end
            
        case 'm'
            %% mental preparation
            %% learning mental
            if strcmp(taskToPerform.mental.learning,'on')
                
                %% first learn the mapping for answering left/right <5/>5 with and then without the display on the screen
                showTitlesInstruction(scr,stim,'learning',true)
                mentalE_prm_learning_and_calib = mental_effort_parameters();
                mentalE_prm_learning_and_calib.startAngle = 0; % for learning always start at zero
                % no time limit for each trial: as long as needed until learning is ok
                learning_useOfTimeLimit = false;
                
                % for learning display the mapping after 2 errors, avoid displaying
                learning_errorLimits.useOfErrorThreshold = false; % no error limit for the learning period
                learning_errorLimits.useOfErrorMapping = true;
                learning_errorLimits.errorMappingLimit = 2; % display mapping after this number of errors
                % extract numbers to use for each learning phase
                nMentalLearning_totalTrials = n_learningInstructions;
                [numberVector_learning] = mental_numbers(nMentalLearning_totalTrials);
                jLearningSession = 0;
                jMentalLearningTrial = 0;
                for iLearning_Instructions = 1:n_learningInstructions
                    jMentalLearningTrial = jMentalLearningTrial + 1;
                    curr_learning_instructions = learning_instructions{iLearning_Instructions};
                    
                    jLearningSession = jLearningSession + 1;
                    learning_sess_nm = ['learning_session',num2str(jLearningSession)];
                    % display instructions for the current learning type
                    [onsets.endLearningInstructions.(learning_sess_nm).(curr_learning_instructions)] = mental_learningInstructions(scr, stim,...
                        curr_learning_instructions, mentalE_prm_learning_and_calib);
                    
                    % perform the learning
                    [learningPerfSummary_Em.(learning_sess_nm).(curr_learning_instructions)] = mental_effort_perf(scr, stim, key_Em,...
                        numberVector_learning(jLearningSession,:),...
                        mentalE_prm_learning_and_calib, n_maxLearning.learning_withInstructions,...
                        curr_learning_instructions, learning_useOfTimeLimit, [], learning_errorLimits);
                    
                    
                    % for experimenter display how many trials have been performed
                    disp(['Mental learning trial ',num2str(jMentalLearningTrial),'/',num2str(nMentalLearning_totalTrials),' done']);
                end % learning instructions loop
                
                %% extended learning for each difficulty level (in N-back version now)
                mentalE_prm_extendedLearning = mentalE_prm_learning_and_calib;
                % always start at zero
                mentalE_prm_extendedLearning.startAngle = 0;
                n_extendedLearningTrials = 30; % number of repetitions
                % Nback version
                Nback_str = num2str(mentalE_prm_extendedLearning.Nback);
                learningVersion = ['extendedLearning_Nback',Nback_str];
                % time limits
                extendedLearning_useOfTimeLimit = true;
                extendedLearning_timeLimit = trainingTimes_Em.max_effort;
                
                % define conditions for the extended learning
                n_maxToReachForCalib = mentalE_prm_extendedLearning.n_maxToReachCalib;
                [numberVector_extendedLearning] = mental_numbers(n_extendedLearningTrials);
                % error handling for extended learning
                extendedLearning_errorLimits.useOfErrorThreshold = false;
                extendedLearning_errorLimits.useOfErrorMapping = false;
                % start at zero
                nMaxReachedUntilNowLearning = 0;
                
                % perform the extended training
                [onsets.endLearningInstructions.(['extendedLearning_session',num2str(1)]).extendedLearning] = mental_learningInstructions(scr, stim,...
                    learningVersion, mentalE_prm_learning_and_calib);
                for iExtendedLearningTrial = 1:n_extendedLearningTrials
                    mentalE_extendedLearningPerfSummary_tmp = mental_effort_perf_Nback(scr, stim, key_Em,...
                        numberVector_extendedLearning(iExtendedLearningTrial,:),...
                        mentalE_prm_extendedLearning, n_maxToReachForCalib,...
                        'noInstructions', extendedLearning_useOfTimeLimit, extendedLearning_timeLimit, extendedLearning_errorLimits, nMaxReachedUntilNowLearning);
                    learningPerfSummary_Em.extendedLearning.(['trial_',num2str(iExtendedLearningTrial)]) = mentalE_extendedLearningPerfSummary_tmp;
                    
                    % extract new best performance
                    nMaxReachedUntilNowLearning = max(nMaxReachedUntilNowLearning, mentalE_extendedLearningPerfSummary_tmp.n_correctAnswersProvided);
                    % small break between each answer
                    DrawFormattedText(window, stim.training.Em.endTrialMsg.text,'center',yScreenCenter/2,white);
                    DrawFormattedText(window,stim.training.Em.endTrialMsg_bis.text,'center','center',white);
                    [~,~,onsets.timeExtendedLearningFbk.(['trial_',num2str(iExtendedLearningTrial)])] = Screen(window,'Flip');
                    WaitSecs(learningTimes_Em.learning_rest);
                    disp(['Mental extended learning trial ',num2str(iExtendedLearningTrial),'/',num2str(n_extendedLearningTrials),' done']);
                end % trial loop
            end
            
            %% calibration mental
            if strcmp(taskToPerform.mental.calib,'on')
                mentalE_prm_learning_and_calib = mental_effort_parameters();
                mentalE_prm_learning_and_calib.startAngle = 0; % for learning always start at zero
                n_calibMax = mentalE_prm_learning_and_calib.n_maxToReachCalib;
                % extract numbers to use for each calibration trial
                %     [numberVector_calib] = mental_numbers(n_calibTrials_Em);
                [numberVector_calib] = mental_calibNumberVector(n_calibTrials_Em, n_calibMax);
                % error handling for extended learning
                extendedLearning_errorLimits.useOfErrorThreshold = false;
                extendedLearning_errorLimits.useOfErrorMapping = false;
                % alternatively, use fixed number of correct answers to provide for each effort
                % level
                % repeat calibration until the subject performance is better
                % than the requested time threshold
                calibSuccess = false;
                iCalibSession = 0;
                while calibSuccess == false
                    iCalibSession = iCalibSession + 1;
                    [n_mental_max_perTrial, calib_summary] = mental_calibNumbers(scr, stim, key_Em,...
                        numberVector_calib, mentalE_prm_learning_and_calib, n_calibTrials_Em, calibTimes_Em);
                    calibSummary.(['calibSession_',num2str(iCalibSession)]).calibSummary = calib_summary;
                    calibSummary.(['calibSession_',num2str(iCalibSession)]).n_mental_max_perTrial = n_mental_max_perTrial;
                end
            end
            
            %% training mental
            if strcmp(taskToPerform.mental.training,'on')
                showTitlesInstruction(scr,stim,'training',true)
                trainingTimes_Em.max_effort = t_min_calib*trainingTimes_Em.t_min_scalingFactor; % allow more time then min performance
                for iTrainingCondition = 1:n_trainingConditions
                    trainingCond = trainingConditions{iTrainingCondition};
                    % display confidence mapping only for first training sessions and
                    % only for fMRI experiment
                    if strcmp(trainingCond,'RP') || IRM == 0 || n_buttonsChoice == 2
                        confidenceChoiceDisplay = false;
                    elseif ~strcmp(trainingCond,'RP') && IRM == 1 && n_buttonsChoice == 4
                        confidenceChoiceDisplay = true;
                    end
                    
                    % define parameters for the training
                    % reward/punishment and effort levels
                    [trainingChoiceOptions_Em_tmp , n_trainingTrials_Em_tmp] = training_options(trainingCond, n_R_levels, n_E_levels,R_money);
                    
                    % start with reward training alone
                    Em_vars_training.i_sub = iSubject;
                    Em_vars_training.n_to_reach = n_to_reach;
                    error('n_to_reach needs updating');
                    % for training: no failures, no display of mapping
                    Em_vars_training.errorLimits.useOfErrorMapping = false;
                    Em_vars_training.errorLimits.useOfErrorThreshold = false;
                    Em_vars_training.timeRemainingEndTrial_ONOFF = 0;
                    [onsets_Em_training.(trainingCond)] = choice_and_perf_trainingInstructions(scr, stim, trainingCond, trainingTimes_Em.instructions);
                    [trainingSummary_Em.(trainingCond)] = choice_and_perf(scr, stim, key_Em, 'mental', Em_vars_training, R_money,...
                        trainingCond, n_trainingTrials_Em_tmp, trainingChoiceOptions_Em_tmp, confidenceChoiceDisplay,...
                        trainingTimes_Em,...
                        results_folder, file_nm_training_Em);
                end % learning condition loop
            end
    end
    %change m/p for next loop
    if strcmp(p_or_m,'p')
        p_or_m = 'm';
    elseif strcmp(p_or_m,'m')
        p_or_m = 'p';
>>>>>>> Stashed changes
    end
end
% reset m/p
if strcmp(p_or_m,'p')
    p_or_m = 'm';
elseif strcmp(p_or_m,'m')
    p_or_m = 'p';
end
%% actual task
if strcmp(taskToPerform.physical.task,'on') || strcmp(taskToPerform.mental.task,'on')
    
    %     waitSpace(langage, window, yScreenCenter, scr, key_Em);
    for iTimeLoop = 1:2
        DrawFormattedText(window,...
            stim.staircase.text,...
            stim.staircase.x,...
            stim.staircase.y,...
            stim.staircase.colour, scr.wrapat);
        if iTimeLoop == 1 % force them to read at first
            [~, onsets.trainingWillStart] = Screen(window, 'Flip');
            WaitSecs(trainingTimes_Ep.instructions);
        elseif iTimeLoop == 2 % after t_instructions seconds, they can manually start
            % display text: Press when you are ready to start
            DrawFormattedText(window, stim.pressWhenReady.text,...
                stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour);
            [~, onsets.trainingWillStart_bis] = Screen(window, 'Flip');
            KbWait;
        end
    end % loop over forced reading/manual pass loop
    
    % for physical effort
    if strcmp(taskToPerform.physical.task,'on')
        Ep_vars.MVC = MVC;
        Ep_vars.dq = dq;
        Ep_vars.Ep_time_levels = Ep_time_levels;
        Ep_vars.F_threshold = F_threshold;
        Ep_vars.F_tolerance = F_tolerance;
        Ep_vars.timeRemainingEndTrial_ONOFF = 0;
    end
    
    % keep track of which block was first
    switch p_or_m
        case 'm'
            % started with mental
            all.mentalfirst = 1;
        case 'p'
            % started with physical
            all.mentalfirst = 0;
    end
    
    
    % repeat twice all of it to have good measures
    for iRepeat=1:nbRepeat
        for i_pm = 1:2
            for iEffortLevel=1:nbEffortLvl
                
                
                % for mental effort timing
                if strcmp(taskToPerform.mental.task,'on')
                    taskTimes_Em.max_effort     = t_min_calib*taskTimes_Em.t_min_scalingFactor; % allow more time then min performance
                end
                
                % keep track of the current session for mental and physical (only important for saving)
                iMental = 1;
                iPhysical = 1;
                
                % number of runs
                for iSession = 1:n_sessions
                    
                    % define if they will see a punishment or reward trial
                    % always start with a reward, they have to win money first
                    session_nm = ['session_nb',num2str(iSession)];
                    if ismember(iSession,[1,2])
                        R_or_P = 'R';
                    elseif ismember(iSession,[3,4])
                        R_or_P = 'P';
                    end
                    
                    % instruction that main task will start soon
                    % for mental effort timing
                    if strcmp(taskToPerform.mental.task,'on')
                        [onsets.taskWillStart] = waitSpace(langage, window, yScreenCenter, scr, key_Em);
                    end
                    
                    % iSubject makes sure each new subject has pattern physical mental in a different order.
                    % iSession switches which code at each session
                    switch p_or_m
                        case 'p'
                            if strcmp(taskToPerform.physical.task,'on')
                                showTitlesInstruction(scr,stim,'task',false);
                                waitSpace(langage, window, yScreenCenter, scr, key_Ep);
                                
                                % run physical task
                                [perfSummary.physical.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iPhysical)]).(['Effort_lvl',(num2str(iEffortLevel))])] = choice_and_perf_staircase(scr, stim, key_Ep,...
                                    'physical', Ep_vars, R_money,...
                                    'mainTask',R_or_P,E_right(iEffortLevel),E_left(iEffortLevel), n_trialsPerSession, taskTimes_Ep,...
                                    results_folder, [file_nm,'_physical_repeat_nb',num2str(iRepeat),'_session_nb',session_nm,'_effort_lvl',num2str(iEffortLevel)]);
                                finalGain = perfSummary.physical.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iPhysical)]).(['Effort_lvl',(num2str(iEffortLevel))]).totalGain(end);
                                iPhysical = iPhysical +1;
                                
                            end
                            
                        case 'm'
                            if strcmp(taskToPerform.mental.task,'on')
                                showTitlesInstruction(scr,stim,'task',true);
                                
                                
                                % instructions
                                mentalE_prm_instruDisplay = mental_effort_parameters();
                                % Nback version
                                Nback_str = num2str(mentalE_prm_instruDisplay.Nback);
                                learningVersion = ['extendedLearning_Nback',Nback_str];
                                [onset_Press] = mental_learningInstructions(scr, stim, learningVersion, mentalE_prm_instruDisplay);
                                Em_vars.i_sub = iSubject;
                                Em_vars.n_to_reach = 0;
                                error('n_to_reach needs updating');
                                % run mental task
                                % for actual task: no display of mapping but consider 3
                                % errors as a trial failure
                                Em_vars.errorLimits.useOfErrorMapping = false;
                                Em_vars.errorLimits.useOfErrorThreshold = false;
                                Em_vars.errorLimits.errorThreshold = 3;
                                
                                Em_vars.timeRemainingEndTrial_ONOFF = 0;
                                [perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iMental)]).(['Effort_lvl',(num2str(iEffortLevel))])] = choice_and_perf_staircase(scr, stim, key_Em,...
                                    'mental', Em_vars, R_money,...
                                    'mainTask',R_or_P,E_right((iEffortLevel)),E_left((iEffortLevel)), n_trialsPerSession, taskTimes_Em,...
                                    results_folder, [file_nm,'_mental_repeat_nb',num2str(iRepeat),'_session_nb',session_nm,'_effort_lvl',num2str(iEffortLevel)]);
                                
                                finalGain = perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iMental)]).(['Effort_lvl',(num2str(iEffortLevel))]).totalGain(end);
                                iMental = iMental +1;
                                
                            end
                    end
                    % Display feedback only for mental effort but in the end for physical too
                    if strcmp(taskToPerform.physical.task,'off') && strcmp(taskToPerform.mental.task,'on')
                        if mod(iSubject+iSession,2) == 0
                            
                            % display feedback for the current session
                            finalGain_str = sprintf('0.2%',finalGain);
                            switch langage
                                case 'fr'
                                    DrawFormattedText(window,...
                                        ['Felicitations! Cette session est maintenant terminee.',...
                                        'Vous avez obtenu: ',finalGain_str,' chf au cours de cette session.'],...
                                        'center', yScreenCenter*(5/3), scr.colours.white, scr.wrapat);
                                case 'engl'
                                    DrawFormattedText(window,...
                                        ['Congratulations! This session is now completed.',...
                                        'You got: ',finalGain_str,' chf during this session.'],...
                                        'center', yScreenCenter*(5/3), scr.colours.white, scr.wrapat);
                            end
                            Screen(window,'Flip');
                            % give break after 4 IP
                            WaitSecs(t_endSession);
                        end
                    end
                    totalGain = totalGain + finalGain;
                    finalGain = 0;
                end % session loop
            end % end of all sessions for 1 effort lvl
            %change m/p for next loop
            if strcmp(p_or_m,'p')
                p_or_m = 'm';
            else
                p_or_m = 'p';
            end
            
        end
        % reset m/p
        if strcmp(p_or_m,'p')
            p_or_m = 'm';
        elseif strcmp(p_or_m,'m')
            p_or_m = 'p';
        end
    end
end

%% save the data

% learning performance
if strcmp(taskToPerform.mental.learning,'on')
    all.physical.learning = learningPerfSummary_Em;
end
if strcmp(taskToPerform.physical.learning,'on')
    all.mental.learning = learningPerfSummary_Ep;
end

% training performance
if strcmp(taskToPerform.physical.training,'on')
    all.physical.training.performance = trainingSummary_Ep;
    all.physical.training.instructions = onsets_Ep_training;
end
if strcmp(taskToPerform.mental.training,'on')
    all.mental.training.performance = trainingSummary_Em;
    all.mental.training.instructions = onsets_Em_training;
end

% actual performance in the main task sessions
% record physical main task data
if strcmp(taskToPerform.physical.task,'on')
    for iEffort= 1:nbEffortLvl % for all effort levels
        for iSession = 1:n_sessions/2 % for the physical sessions
            
            for iRepeat = 1:nbRepeat
                % save data in all and reformat it in a specific order
                all.physical.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).(['repeat_nb',num2str(iRepeat)]).perfSummary = perfSummary.physical.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',(num2str(iEffort))]);
                IP_variables.physicalDeltaIP = all.physical.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).(['repeat_nb',num2str(iRepeat)]).perfSummary.IP - baselineR;
            end
        end
    end
end
% record mental main task data
if strcmp(taskToPerform.mental.task,'on')
    for iEffort= 1:nbEffortLvl % for all effort levels
        for iSession = 1:n_sessions/2 % for the mental sessions
            for iRepeat = 1:nbRepeat
                % save data in all and reformat it in a specific order
                all.mental.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).(['repeat_nb',num2str(iRepeat)]).perfSummary = perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',(num2str(iEffort))]);
                IP_variables.mentalDeltaIP = all.mental.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).(['repeat_nb',num2str(iRepeat)]).perfSummary.IP - baselineR;
            end
        end
    end
end

IP_variables.baselineR = baselineR;
IP_variables.baselineP = baselineP;
IP_variables.training.p_or_m = p_or_m;
IP_variables.training.MVC = MVC;
% actually save the data
save([results_folder, file_nm,'.mat']);

% save delta_IP and baselineR
save([results_folder, file_nm_IP,'.mat'],'IP_variables');

%% Show a final screen if and only if they performed the task or nonsense since no amount involved
if strcmp(taskToPerform.physical.task,'on') || strcmp(taskToPerform.mental.task,'on')
    totalGain_str = sprintf('%0.2f',totalGain);
    % display feedback for the current session
    switch langage
        case 'fr'
            DrawFormattedText(window,...
                ['Felicitations! Cette experience est maintenant terminee.',...
                'Vous avez obtenu: ',totalGain_str,' chf au cours de cette session.'],...
                'center', 'center', scr.colours.white, scr.wrapat);
        case 'engl'
            DrawFormattedText(window,...
                ['Congratulations! This session is now completed.',...
                'You got: ',totalGain_str,' chf during this session.'],...
                'center', 'center', scr.colours.white, scr.wrapat);
    end
    Screen(window,'Flip');
    WaitSecs(t_endSession);
end
%% close PTB
ShowCursor;
sca;