% script for training participants before performing the main experiment

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
subResultFolder = [results_folder, subjectCodeName, filesep,'behavior',filesep];
if ~exist(subResultFolder,'dir')
    mkdir(subResultFolder);
end

if strcmp(p_or_m,'m') == 0 && strcmp(p_or_m,'p') == 0
    error('this letter has no definition');
end

file_nm_training_Em = ['training_data_Em_CID',num2str(iSubject)];
file_nm_training_Ep = ['training_data_Ep_CID',num2str(iSubject)];
file_nm = ['training_data_CID',num2str(iSubject)];
Ep_calib_filenm = [subResultFolder,subjectCodeName,'_physicalCalib.mat'];
Em_calib_filenm = [subResultFolder,subjectCodeName,'_mentalCalib.mat'];
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
taskToPerform.physical.task = 'off';
taskToPerform.mental.learning_1 = 'on';
taskToPerform.mental.calib = 'on';
taskToPerform.mental.learning_2 = 'on';
taskToPerform.mental.training = 'on';
taskToPerform.mental.task = 'on';
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
ShowCursor;
white = scr.colours.white;
black = scr.colours.black;

% include punishment condition?
punishment_yn = 'yes'; % include punishment trials?

% number of reward and effort conditions
n_R_levels = 4;
n_E_levels = 4;

% initialize visual stimuli to use in the experiment
[stim] = stim_initialize(scr, n_E_levels, langage);

% define number of training conditions
switch punishment_yn
    case 'yes'
        trainingConditions = {'RP'};
    case 'no'
        trainingConditions = {'R'};
end
trainingCondition = trainingConditions{1};

% define number of training trials (ie when choice + perf)
n_trainingTrials = 4;
% define number of trials per staircase procedure
n_trialsPerSession = 5;
% load timings for each phase of the experiment
[trainingTimes_Em, calibTimes_Em, learningTimes_Em, taskTimes_Em] = timings_definition(trainingConditions, n_trialsPerSession, n_trainingTrials, 'mental');
[trainingTimes_Ep, calibTimes_Ep, learningTimes_Ep, taskTimes_Ep] = timings_definition(trainingConditions, n_trialsPerSession, n_trainingTrials, 'physical');

% number of times we apply a staircase procedure. (2 mental (R or P) 2 physical (R or P))
if strcmp(taskToPerform.physical.task,'on') &&...
        strcmp(taskToPerform.mental.task,'on')
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

% mental calibration: consider calibration failed if more than
% errorThreshold errors have been made during the calibration => do again
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
    F_threshold = 55; % force should be maintained above this threshold (expressed in % of MVC)
    F_tolerance = 2.5; % tolerance allowed around the threshold (expressed in % of MVC)
    % need to define timings for each level of force
    [Ep_time_levels] = physical_effortLevels(n_E_levels);
end

%% mental parameters
if strcmp(taskToPerform.mental.calib,'on') ||...
        strcmp(taskToPerform.mental.learning_1,'on') ||...
        strcmp(taskToPerform.mental.learning_2,'on') ||...
        strcmp(taskToPerform.mental.training,'on') ||...
        strcmp(taskToPerform.mental.task,'on')
    % define relevant keys and dynamometer module
    key_Em = relevant_key_definition('mental', IRM, n_buttonsChoice);
end

%% run the code twice, each time for p or m conditions
for i_pm = 1:2
    switch p_or_m
        case 'p'
            %% physical MVC (calibrate the Fmax for the whole experiment)
            if strcmp(taskToPerform.physical.calib,'on')
                n_MVC_repeat = 3;
                waitSpace(langage, window, yScreenCenter, scr, key_Ep);
                [MVC_tmp, onsets_MVC] = physical_effort_MVC(scr, stim, dq, n_MVC_repeat, calibTimes_Ep);
                MVC = max(MVC_tmp.MVC); % expressed in Voltage
                save(Ep_calib_filenm,'MVC');
            elseif strcmp(taskToPerform.physical.calib,'off') &&...
                    ( strcmp(taskToPerform.physical.learning,'on') ||...
                    strcmp(taskToPerform.physical.training,'on') ||...
                    strcmp(taskToPerform.physical.task,'on') )
                MVC = getfield(load(Ep_calib_filenm,'MVC'),'MVC');
            end
            
            %% learning physical (learn each level of force)
            if strcmp(taskToPerform.physical.learning,'on')
                n_Ep_learningForceRepeats = 1; % number of learning repetitions for each level of difficulty (= each level of force)
                % introduce physical learning
                showTitlesInstruction(scr,stim,'learning',p_or_m);
                waitSpace(langage, window, yScreenCenter, scr, key_Ep);
                % perform physical learning
                [learningPerfSummary_Ep, learningOnsets_Ep] = physical_learning(scr, stim, dq, n_E_levels, Ep_time_levels,...
                    F_threshold, F_tolerance, MVC,...
                    n_Ep_learningForceRepeats, learningTimes_Ep);
            end
            
            %% training physical (choice + effort)
            if strcmp(taskToPerform.physical.training,'on')
                
                % introduce physical training
                showTitlesInstruction(scr,stim,'training',p_or_m);
                waitSpace(langage, window, yScreenCenter, scr, key_Ep);
                % define parameters for the training
                Ep_vars_training.MVC = MVC;
                Ep_vars_training.dq = dq;
                Ep_vars_training.Ep_time_levels = Ep_time_levels;
                Ep_vars_training.F_threshold = F_threshold;
                Ep_vars_training.F_tolerance = F_tolerance;
                Ep_vars_training.timeRemainingEndTrial_ONOFF = 0;
                % mapping between reward levels and fake monetary amounts
                R_money = R_amounts(n_R_levels, punishment_yn);
                % reward/punishment and effort levels
                [trainingChoiceOptions_Ep_tmp] = training_options(trainingCondition, n_R_levels, n_E_levels, R_money, n_trainingTrials);
                % perform physical training in 2 phases:
                % 1) with confidence mapping;
                % 2) without confidence mapping
                trainingConfConditions = {'withConfMapping','withoutConfMapping'};
                n_trainingConfConditions = length(trainingConfConditions); % with/without confidence mapping
                if n_trainingConfConditions == 2
                    n_trialsPerTrainingCondition = n_trainingTrials/2;
                    if floor(n_trialsPerTrainingCondition) < (n_trainingTrials/2)
                        error('number of training trials should be pair. Please fix your training design matrix.');
                    end
                end
                % perform the training
                for iTrainingCondition = 1:n_trainingConfConditions
                    trainingConfCond = trainingConfConditions{iTrainingCondition};
                    % display confidence mapping only for first training sessions and
                    % only for fMRI experiment
                    if strcmp(trainingConfCond,'withoutConfMapping') || n_buttonsChoice == 2
                        confidenceChoiceDisplay = false;
                    elseif strcmp(trainingConfCond,'withConfMapping') && n_buttonsChoice == 4
                        confidenceChoiceDisplay = true;
                    end
                    % extract trials to use
                    trainingTrials_idx = (1:n_trialsPerTrainingCondition) + n_trialsPerTrainingCondition*(iTrainingCondition - 1);
                    % training instruction
                    [onsets_Ep_training.(['session',num2str(iTrainingCondition),'_',trainingConfCond])] = choice_and_perf_trainingInstructions(scr, stim, trainingCondition, trainingTimes_Ep.instructions);
                    % perform the training
                    [trainingSummary_Ep.(['session',num2str(iTrainingCondition),'_',trainingConfCond])] = choice_and_perf(scr, stim, key_Ep, 'physical', Ep_vars_training,...
                        trainingCondition, n_trialsPerTrainingCondition, trainingChoiceOptions_Ep_tmp, confidenceChoiceDisplay,...
                        trainingTimes_Ep,...
                        subResultFolder, file_nm_training_Ep);
                end % learning condition loop
                
                DrawFormattedText(window, stim.training.Ep.endMsg.text,...
                    'center','center',stim.training.Ep.endMsg.colour, scr.wrapat);
                [~,onsets.EndTrainingMsg] = Screen('Flip',window); % display the cross on screen
                WaitSecs(trainingTimes_Ep.trainingEnd);
            end
            
        case 'm'
            %% learning mental: 0-back and 2-back as a calibration
            if strcmp(taskToPerform.mental.learning_1,'on')
                
                %% learning the mapping for answering left/right <5/>5 with and then without the display on the screen (0-back)
                % learning parameters
                % perform 2 short learning sessions: one with mapping
                % (left/right) - (<5/>5) and one without to associate
                % left/right answer with a given side on the screen
                learning1_instructions = {'fullInstructions','noInstructions'}; %,'partialInstructions'
                n_learningInstructions = length(learning1_instructions);
                % initial learning: careful to enter a pair number here
                n_maxLearning.learning_withInstructions = 15;
                n_maxLearning.learning_withoutInstructions = 15;
                mentalE_prm_learning = mental_effort_parameters();
                mentalE_prm_learning.startAngle = 0; % for learning always start at zero
                % no time limit for each trial: as long as needed until learning is ok
                learning_useOfTimeLimit = false;
                % display instructions for learning
                showTitlesInstruction(scr,stim,'learning',p_or_m)
                
                % for learning display the mapping after 'errorMappingLimit' number of errors
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
                    curr_learning_instructions = learning1_instructions{iLearning_Instructions};
                    
                    jLearningSession = jLearningSession + 1;
                    learning_sess_nm = ['learning_0back_session',num2str(jLearningSession)];
                    % display instructions for the current learning type
                    [onsets.endLearningInstructions.(learning_sess_nm).(curr_learning_instructions)] = mental_learningInstructions(scr, stim,...
                        curr_learning_instructions, mentalE_prm_learning);
                    
                    % perform the learning
                    [learningPerfSummary_Em.(learning_sess_nm).(curr_learning_instructions)] = mental_effort_perf(scr, stim, key_Em,...
                        numberVector_learning(jLearningSession,:),...
                        mentalE_prm_learning, n_maxLearning.learning_withInstructions,...
                        curr_learning_instructions, learning_useOfTimeLimit, [], learning_errorLimits);
                    
                    
                    % for experimenter display how many trials have been performed
                    disp(['Mental 0-back learning trial ',num2str(jMentalLearningTrial),'/',num2str(nMentalLearning_totalTrials),' done']);
                end % learning instructions loop
                
                %% learning (1) by repeating the calibration many times before actual calibration
                % define number of trials to perform
                n_learning1calibLikeTrials = 2;
                
                mentalE_prm_learning1calibLike = mental_effort_parameters();
                % always start at zero
                mentalE_prm_learning1calibLike.startAngle = 0;
                % Nback version
                Nback_str = num2str(mentalE_prm_learning1calibLike.Nback);
                learningVersion = ['learning_Nback',Nback_str];
                % time limits
                learning1calibLike_useOfTimeLimit = true;
                learning1calibLike_timeLimit = trainingTimes_Em.max_effort;
                
                % define conditions for the learning
                n_maxToReachForCalib = mentalE_prm_learning1calibLike.n_maxToReachCalib;
                [numberVector_learning1calibLike] = mental_numbers(n_learning1calibLikeTrials);
                % error handling for learning
                learning1calibLike_errorLimits.useOfErrorThreshold = false;
                learning1calibLike_errorLimits.useOfErrorMapping = false;
                % start at zero so that they see the orange bar improving
                % from trial to trial
                nMaxReachedUntilNowLearning = 0;
                % use a minimal amount of correct answers to give or the
                % trial will be repeated
                n_Em_learning1calibLike_MinToReach = 6;
                
                % perform the learning session
                [onsets.endLearningInstructions.learning1calibLike_session] = mental_learningInstructions(scr, stim,...
                    learningVersion, mentalE_prm_learning1calibLike);
                n_maxReachedDuringLearning = NaN(1,n_learning1calibLikeTrials);
                for iLearning1Trial = 1:n_learning1calibLikeTrials
                    mentalE_learning1calibLikePerfSummary_tmp = mental_effort_perf_Nback(scr, stim, key_Em,...
                        numberVector_learning1calibLike(iLearning1Trial,:),...
                        mentalE_prm_learning1calibLike, n_maxToReachForCalib,...
                        'noInstructions', learning1calibLike_useOfTimeLimit, learning1calibLike_timeLimit,...
                        learning1calibLike_errorLimits, nMaxReachedUntilNowLearning, n_Em_learning1calibLike_MinToReach);
                    learningPerfSummary_Em.learning1calibLike.(['trial_',num2str(iLearning1Trial)]) = mentalE_learning1calibLikePerfSummary_tmp;
                    n_maxReachedDuringLearning(iLearning1Trial) = mentalE_learning1calibLikePerfSummary_tmp.n_correctAnswersForDisplay;
                    
                    % extract new best performance
                    nMaxReachedUntilNowLearning = max(nMaxReachedUntilNowLearning, n_maxReachedDuringLearning(iLearning1Trial));
                    % small break between each answer
                    DrawFormattedText(window, stim.training.Em.endTrialMsg.text,'center',yScreenCenter/2,white);
                    DrawFormattedText(window,stim.training.Em.endTrialMsg_bis.text,'center','center',white);
                    [~,~,onsets.timeLearningFbk.(['trial_',num2str(iLearning1Trial)])] = Screen(window,'Flip');
                    WaitSecs(learningTimes_Em.learning_rest);
                    disp(['Mental learning calibration-like trial ',num2str(iLearning1Trial),'/',num2str(n_learning1calibLikeTrials),' done']);
                end % trial loop
                
                learning1done = 0;
                n_learning1bonusTrialsToLearn = 10;
                jLearningTrial = n_learning1calibLikeTrials;
                while learning1done == 0
                    if n_learning1calibLikeTrials > 10 &&...
                            ((n_maxReachedDuringLearning(end) < n_Em_learning1calibLike_MinToReach) ||...
                            (n_maxReachedDuringLearning(end-1) < n_Em_learning1calibLike_MinToReach) ||...
                            (n_maxReachedDuringLearning(end-2) < n_Em_learning1calibLike_MinToReach) ||...
                            (n_maxReachedDuringLearning(end-3) < n_Em_learning1calibLike_MinToReach) ||...
                            (n_maxReachedDuringLearning(end-4) < n_Em_learning1calibLike_MinToReach) )
                        disp('performance was too low in one of the last trials. We will redo ',...
                            num2str(n_learning1bonusTrialsToLearn),' more trials to compensate.');
                        [numberVector_learning1_bonus] = mental_numbers(n_learning1bonusTrialsToLearn);
                        for iLearning1Trial_bonus = 1:n_learning1bonusTrialsToLearn
                            jLearningTrial = jLearningTrial + 1;
                            mentalE_learning1calibLikePerfSummary_tmp = mental_effort_perf_Nback(scr, stim, key_Em,...
                                numberVector_learning1_bonus(iLearning1Trial_bonus,:),...
                                mentalE_prm_learning1calibLike, n_maxToReachForCalib,...
                                'noInstructions', learning1calibLike_useOfTimeLimit, learning1calibLike_timeLimit,...
                                learning1calibLike_errorLimits, nMaxReachedUntilNowLearning, n_Em_learning1calibLike_MinToReach);
                            learningPerfSummary_Em.learning1calibLike.(['trial_',num2str(jLearningTrial)]) = mentalE_learning1calibLikePerfSummary_tmp;
                            n_maxReachedDuringLearning(jLearningTrial) = mentalE_learning1calibLikePerfSummary_tmp.n_correctAnswersForDisplay;
                            
                            % extract new best performance
                            nMaxReachedUntilNowLearning = max(nMaxReachedUntilNowLearning, n_maxReachedDuringLearning(jLearningTrial));
                            % small break between each answer
                            DrawFormattedText(window, stim.training.Em.endTrialMsg.text,'center',yScreenCenter/2,white);
                            DrawFormattedText(window,stim.training.Em.endTrialMsg_bis.text,'center','center',white);
                            [~,~,onsets.timeLearningFbk.(['trial_',num2str(jLearningTrial)])] = Screen(window,'Flip');
                            WaitSecs(learningTimes_Em.learning_rest);
                            disp(['Mental learning calibration-like BONUS trial ',num2str(iLearning1Trial_bonus),'/',num2str(n_learning1bonusTrialsToLearn),' done']);
                        end % trial loop
                    else
                        learning1done = 1;
                    end
                end
            end % learning (1)
            
            %% calibration mental
            if strcmp(taskToPerform.mental.calib,'on')
                % number of calibration trials
                n_calibTrials_Em = 3;
                
                % mental calibration parameters
                mentalE_prm_calib = mental_effort_parameters();
                mentalE_prm_calib.startAngle = 0; % for learning always start at zero
                n_calibMax = mentalE_prm_calib.n_maxToReachCalib;
                % extract numbers to use for each calibration trial
                %     [numberVector_calib] = mental_numbers(n_calibTrials_Em);
                [numberVector_calib] = mental_calibNumberVector(n_calibTrials_Em, n_calibMax);
                % perform the calibration
                [NMP, calib_summary] = mental_calibNumbers(scr, stim, key_Em,...
                    numberVector_calib, mentalE_prm_calib, n_calibTrials_Em, calibTimes_Em, langage);
                calibSummary.calibSummary = calib_summary;
                calibSummary.n_mental_max_perTrial = NMP;
                % record number of maximal performance (NMP)
                save(Em_calib_filenm,'NMP');
            elseif strcmp(taskToPerform.mental.calib,'off') &&...
                    ( strcmp(taskToPerform.mental.learning_1,'on') ||...
                    strcmp(taskToPerform.mental.learning_2,'on') ||...
                    strcmp(taskToPerform.mental.training,'on') ||...
                    strcmp(taskToPerform.mental.task,'on') )
                NMP = getfield(load(Em_calib_filenm,'NMP'),'NMP');
            end % calibration
            
            %% learning (2) for each difficulty level
            if strcmp(taskToPerform.mental.learning_2,'on')
                % introduce physical learning
                showTitlesInstruction(scr,stim,'learning','m');
                waitSpace(langage, window, yScreenCenter, scr, key_Em);
                
                % define all difficulty levels based on calibration
                [n_to_reach] = mental_N_answersPerLevel(n_E_levels, NMP);
                % define number of learning trials
                n_Em_learningForceRepeats = 1; % number of learning repetitions for each level of difficulty (= each level of force)
                % timings
                Em_learningTimings.time_limit = true;
                Em_learningTimings.t_max = learningTimes_Em.max_effort;
                Em_learningTimings.learning_rest = learningTimes_Em.learning_rest;
                % perform all the difficulty levels
                [learning2PerfSummary_Em, onsets] = mental_learning(scr, stim, key_Em, n_E_levels, n_to_reach, n_Em_learningForceRepeats, Em_learningTimings);
            end % learning (2)
            
            %% training mental
            if strcmp(taskToPerform.mental.training,'on')
                % define parameters for the training
                [Em_vars_training.n_to_reach] = mental_N_answersPerLevel(n_E_levels, NMP);
                Em_vars_training.errorLimits.useOfErrorMapping = false;
                Em_vars_training.errorLimits.useOfErrorThreshold = false;
                Em_vars_training.timeRemainingEndTrial_ONOFF = 0;
                % mapping between reward levels and fake monetary amounts
                R_money = R_amounts(n_R_levels, punishment_yn);
                % reward/punishment and effort levels
                [trainingChoiceOptions_Em_tmp] = training_options(trainingCondition, n_R_levels, n_E_levels, R_money, n_trainingTrials);
                % perform physical training in 2 phases:
                % 1) with confidence mapping;
                % 2) without confidence mapping
                trainingConfConditions = {'withConfMapping','withoutConfMapping'};
                n_trainingConfConditions = length(trainingConfConditions); % with/without confidence mapping
                if n_trainingConfConditions == 2
                    n_trialsPerTrainingCondition = n_trainingTrials/2;
                    if floor(n_trialsPerTrainingCondition) < (n_trainingTrials/2)
                        error('number of training trials should be pair. Please fix your training design matrix.');
                    end
                end
                % perform the training
                for iTrainingCondition = 1:n_trainingConfConditions
                    trainingConfCond = trainingConfConditions{iTrainingCondition};
                    % display confidence mapping only for first training sessions and
                    % only for fMRI experiment
                    if strcmp(trainingConfCond,'withoutConfMapping') || n_buttonsChoice == 2
                        confidenceChoiceDisplay = false;
                    elseif strcmp(trainingConfCond,'withConfMapping') && n_buttonsChoice == 4
                        confidenceChoiceDisplay = true;
                    end
                    % extract trials to use
                    trainingTrials_idx = (1:n_trialsPerTrainingCondition) + n_trialsPerTrainingCondition*(iTrainingCondition - 1);
                    % training instruction
                    [onsets_Em_training.(['session',num2str(iTrainingCondition),'_',trainingConfCond])] = choice_and_perf_trainingInstructions(scr, stim, trainingCondition, trainingTimes_Em.instructions);
                    % select effort level
                    [trainingSummary_Em.(trainingConfCond)] = choice_and_perf(scr, stim, key_Em, 'mental', Em_vars_training,...
                        trainingCondition, n_trialsPerTrainingCondition, trainingChoiceOptions_Em_tmp, confidenceChoiceDisplay,...
                        trainingTimes_Em,...
                        subResultFolder, file_nm_training_Em);
                end % training condition loop
            end
    end
    
    %% change m/p for next loop
    if strcmp(p_or_m,'p')
        p_or_m = 'm';
    elseif strcmp(p_or_m,'m')
        p_or_m = 'p';
    end
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

%% indifference point measurement
if strcmp(taskToPerform.physical.task,'on') || strcmp(taskToPerform.mental.task,'on')
    
    % Number of repeats of the whole code (how many times will you measure
    % the IP?)
    nbRepeat = 2;
    % for how many levels of effort will you measure the IP?
    E_right = 2;
    E_left  = 0;
    nbEffortLvl = 1;
    
    % Baseline Reward (CHF)
    baselineR = 0.5;
    baselineP = 0.5;
    
    % Total amount of money to be given
    totalGain = 0;
    sessionFinalGain = 0;
    
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
    
    
    for iEffortLevel = 1:nbEffortLvl
        % keep track of the current session for mental and physical (only important for saving)
        iMental = 1;
        iPhysical = 1;
        
        % number of runs
        for iSession = 1:n_sessions
            % session gains
            sessionFinalGain = 0;
            
            % define if they will see a punishment or reward trial
            % always start with a reward, they have to win money first
            session_nm = ['session_nb',num2str(iSession)];
            R_or_P = 'R'; % compute IP only for rewards
            
            for i_pm = 1:2
                p_or_m
                
                
                % perform physical or mental IP according to inputs
                switch p_or_m
                    case 'p'
                        if strcmp(taskToPerform.physical.task,'on')
                            % instructions
                            showTitlesInstruction(scr,stim,'task',p_or_m);
                            waitSpace(langage, window, yScreenCenter, scr, key_Ep);
                            
                            % run physical task
                            perf_Ep_IP_tmp = choice_and_perf_staircase(scr, stim, key_Ep,...
                                'physical', Ep_vars,...
                                'mainTask',R_or_P,E_right(iEffortLevel),E_left(iEffortLevel), n_trialsPerSession, taskTimes_Ep,...
                                subResultFolder, [file_nm,'_physical_session_nb',session_nm,'_effort_lvl',num2str(iEffortLevel)]);
                            perfSummary.physical.(['session_nb',num2str(iPhysical)]).(['Effort_lvl',(num2str(iEffortLevel))]) = perf_Ep_IP_tmp;
                            sessionFinalGain = sessionFinalGain + perf_Ep_IP_tmp.totalGain(end);
                            iPhysical = iPhysical + 1;
                        end
                        
                    case 'm'
                        if strcmp(taskToPerform.mental.task,'on')
                            % instructions
                            showTitlesInstruction(scr,stim,'task',p_or_m);
                            waitSpace(langage, window, yScreenCenter, scr, key_Em);
                            
                            mentalE_prm_instruDisplay = mental_effort_parameters();
                            % Nback version
                            Nback_str = num2str(mentalE_prm_instruDisplay.Nback);
                            learningVersion = ['learning_Nback',Nback_str];
                            mental_learningInstructions(scr, stim, learningVersion, mentalE_prm_instruDisplay);
                            [Em_vars.n_to_reach] = mental_N_answersPerLevel(n_E_levels, NMP);
                            % run mental task
                            % for actual task: no display of mapping but consider 3
                            % errors as a trial failure
                            Em_vars.errorLimits.useOfErrorMapping = false;
                            Em_vars.errorLimits.useOfErrorThreshold = false;
                            Em_vars.timeRemainingEndTrial_ONOFF = 0;
                            perf_Em_IP_tmp = choice_and_perf_staircase(scr, stim, key_Em,...
                                'mental', Em_vars,...
                                'mainTask',R_or_P,E_right((iEffortLevel)),E_left(iEffortLevel), n_trialsPerSession, taskTimes_Em,...
                                subResultFolder, [file_nm,'_mental_session_nb',session_nm,'_effort_lvl',num2str(iEffortLevel)]);
                            perfSummary.mental.(['session_nb',num2str(iMental)]).(['Effort_lvl',(num2str(iEffortLevel))]) = perf_Em_IP_tmp;
                            sessionFinalGain = sessionFinalGain + perf_Em_IP_tmp.totalGain(end);
                            iMental = iMental +1;
                        end
                end
                if strcmp(p_or_m,'p')
                    p_or_m = 'm';
                else
                    p_or_m = 'p';
                end
            end % session loop
            % display feedback for the current session
            finalGain_str = sprintf('%0.2f',sessionFinalGain);
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
            totalGain = totalGain + sessionFinalGain;
            
            
        end % end of all sessions for 1 effort lvl
        %change m/p for next loop
        
        
    end
    
end

%% save the data

% calibration
if strcmp(taskToPerform.physical.calib,'on')
    all.physical.MVC = MVC;
end
if strcmp(taskToPerform.mental.calib,'on')
    all.physical.NMP = NMP;
end
% learning performance
if strcmp(taskToPerform.mental.learning_1,'on')
    all.mental.learning_1 = learningPerfSummary_Em;
    all.mental.learning_1.n_minToReach = n_Em_learning1calibLike_MinToReach;
end
if strcmp(taskToPerform.mental.learning_2,'on')
    all.mental.learning_2 = learning2PerfSummary_Em;
end
if strcmp(taskToPerform.physical.learning,'on')
    all.mental.learning = learningPerfSummary_Ep;
end

% training performance
if strcmp(taskToPerform.physical.training,'on')
    all.physical.training.performance = trainingSummary_Ep;
    all.physical.training.onsets = onsets_Ep_training;
end
if strcmp(taskToPerform.mental.training,'on')
    all.mental.training.performance = trainingSummary_Em;
    all.mental.training.onsets = onsets_Em_training;
end

% actual performance in the main task sessions
% record physical main task data
if strcmp(taskToPerform.physical.task,'on')
    for iEffort= 1:nbEffortLvl % for all effort levels
        for iSession = 1:n_sessions/2 % for the physical sessions
            % save data in all and reformat it in a specific order
            all.physical.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).perfSummary = perfSummary.physical.(['session_nb',num2str(iSession)]).(['Effort_lvl',(num2str(iEffort))]);
            IP_variables.physicalDeltaIP = all.physical.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).perfSummary.IP - baselineR;
        end
    end
end
% record mental main task data
if strcmp(taskToPerform.mental.task,'on')
    for iEffort= 1:nbEffortLvl % for all effort levels
        for iSession = 1:n_sessions/2 % for the mental sessions
            % save data in all and reformat it in a specific order
            all.mental.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).perfSummary = perfSummary.mental.(['session_nb',num2str(iSession)]).(['Effort_lvl',(num2str(iEffort))]);
            IP_variables.mentalDeltaIP = all.mental.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).perfSummary.IP - baselineR;
        end
    end
end

IP_variables.baselineR = baselineR;
IP_variables.baselineP = baselineP;
IP_variables.training.p_or_m = p_or_m;
if strcmp(taskToPerform.physical.task,'on')
    IP_variables.calibration.MVC = MVC;
end
if strcmp(taskToPerform.mental.task,'on')
    IP_variables.calibration.NMP = NMP;
end
% actually save the data
save([subResultFolder, file_nm,'.mat']);

% save delta_IP and baselineR
save([subResultFolder, file_nm_IP,'.mat'],'IP_variables');

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