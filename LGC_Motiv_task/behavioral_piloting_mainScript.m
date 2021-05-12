% script for behavioral pilots

%% working directories
% launch within the folder where scripts are stored or will not work
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


%% define subject ID

%% general parameters
% initialize screen
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(0, 1);
white = scr.colours.white;
black = scr.colours.black;

% define relevant keys and dynamometer module
[key, dq] = relevant_key_definition('physical', 0);

% include punishment condition?
punishment_yn = 'yes'; % include punishment trials?

% number of reward and effort conditions
n_R_levels = 3;
n_E_levels = 3;
n_trialsPerSession = 44;
% mapping between reward levels and actual monetary amounts
R_money = R_amounts(n_R_levels);

% initialize visual stimuli to use in the experiment
[stim] = stim_initialize(scr, n_R_levels, n_E_levels, pics_folder);

% define number of training conditions
switch punishment_yn
    case 'yes'
        trainingConditions = {'R','P','RP'};
    case 'no'
        trainingConditions = {'R'};
end
n_trainingConditions = length(trainingConditions);

% load timings for each phase of the experiment
[trainingTimes_Em, calibTimes_Em, learningTimes_Em, taskTimes_Em] = timings_definition(scr, trainingConditions, n_R_levels, n_E_levels, n_trialsPerSession, 'mental');
[trainingTimes_Ep, calibTimes_Ep, learningTimes_Ep, taskTimes_Ep] = timings_definition(scr, trainingConditions, n_R_levels, n_E_levels, n_trialsPerSession, 'physical');

%% physical parameters
n_MVC_repeat = 3;
n_learningForceRepeats = 3; % number of learning repetitions for each level of difficulty (= each level of force)
F_threshold = 50; % force should be maintained above this threshold (expressed in % of MVC)
F_tolerance = 2.5; % tolerance allowed around the threshold (expressed in % of MVC)
% need to define timings for each level of force
[Ep_time_levels] = physical_effortLevels(n_E_levels);

% calibration
n_calibTrials_Ep = 5;

%% mental parameters
% define number of pairs to solve for each level of difficulty
n_to_reach = mental_N_answersPerLevel(n_E_levels);

% calibration: calibrate the maximal duration required for the
% top effort
n_calibMax = n_to_reach.(['E_level_',num2str(n_E_levels)]);
n_calibTrials_Em = 5;

% learning
% perform 2 learning sessions, one with instructions and then one without
% (left/right) vs (odd/even) and (lower/higher than 5) - mapping indicated the first time)
% need to remind the mapping the second time
learning_cols = {'col1','col2','all'};
n_learningColours = length(learning_cols);
learning_instructions = {'fullInstructions','noInstructions'}; %,'partialInstructions'
n_learningInstructions = length(learning_instructions);
% initial learning
n_maxLearning.learning_withInstructions = 15;
n_maxLearning.learning_withoutInstructions = 15;

%% physical preparation
    %% physical MVC
[initial_MVC, onsets_initial_MVC] = physical_effort_MVC(scr, dq, n_MVC_repeat, calibTimes_Ep);
MVC = nanmax(initial_MVC.MVC); % expressed in Voltage

    %% learning physical
[learningPerfSummary, learningOnsets] = physical_learning(scr, stim, dq, n_E_levels, Ep_time_levels,...
    F_threshold, F_tolerance, MVC,...
    n_learningForceRepeats, learningTimes_Ep);

    %% training physical
for iTrainingCondition = 1:n_trainingConditions
    trainingCond = trainingConditions{iTrainingCondition};
    
    % define parameters for the training
    % reward/punishment and effort levels
    [trainingChoiceOptions_tmp, n_trainingTrials_tmp, R_or_P_training_tmp] = training_options(trainingCond, n_R_levels, n_E_levels);
    
    % start with reward training alone
    Ep_vars.MVC = MVC;
    Ep_vars.dq = dq;
    [trainingSummary.(trainingCond)] = choice_and_perf_training(scr, stim, key, 'physical', Ep_vars, R_money,...
        trainingCond, R_or_P_training_tmp, n_trainingTrials_tmp, trainingChoiceOptions_tmp, trainingTimes_Ep);
end % learning condition loop

DrawFormattedText(window,'Bravo! Votre entraînement est terminé.',...
    'center','center',scr.colours.black, scr.wrapat);
[~,onsets.EndTrainingMsg] = Screen('Flip',window); % display the cross on screen
WaitSecs(trainingTimes.trainingEnd);

%% mental preparation
    %% learning mental
mentalE_prm_learning_and_calib = mental_effort_parameters(i_sub);
mentalE_prm_learning_and_calib.startAngle = 0; % for learning always start at zero
% no time limit for each trial: as long as needed until learning is
% ok
learning_time_limit = false;
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

    %% calibration mental
% extract numbers to use for each calibration trial
[numberVector_calib] = mental_numbers(n_calibTrials);

% alternatively, use fixed number of correct answers to provide for each effort
% level
% repeat calibration until the subject performance is better
% than the requested time threshold
calibSuccess = false;
calibSession = 0;
while calibSuccess == false
    calibSession = calibSession + 1;
    [t_min_calib, calibSessionSummary, calibSuccess] = mental_calibTime(scr, stim, key,...
        numberVector_calib, mentalE_prm_learning_and_calib, n_calibTrials_Em, n_calibMax, calibTimes_Em);
    calibSummary.(['calibSession_',num2str(calibSession)]).calibSummary = calibSessionSummary;
    calibSummary.(['calibSession_',num2str(calibSession)]).calibSuccess = calibSuccess;
    calibSummary.(['calibSession_',num2str(calibSession)]).t_mental_max_perTrial = t_min_calib;
end
    %% training mental


%% actual task
% each block 1) MVC 2) task 3) MVC

% physical 1

% mental 1

% physical 2

% mental 2

%% close all
sca;
