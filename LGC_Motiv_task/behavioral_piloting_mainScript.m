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
% go back to folder with scripts
cd(main_task_folder);

%% define subject ID
iSubject = 0;
file_nm_training_Em = ['pilot_data_Em_sub',num2str(iSubject)];
file_nm_training_Ep = ['pilot_data_Ep_sub',num2str(iSubject)];
file_nm = ['pilot_data_sub',num2str(iSubject)];

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

n_sessions = 4; % 4 blocks in total (2 mental and 2 physical)
%% physical parameters
n_MVC_repeat = 3;
n_learningForceRepeats = 3; % number of learning repetitions for each level of difficulty (= each level of force)
F_threshold = 50; % force should be maintained above this threshold (expressed in % of MVC)
F_tolerance = 2.5; % tolerance allowed around the threshold (expressed in % of MVC)
% need to define timings for each level of force
[Ep_time_levels] = physical_effortLevels(n_E_levels);

% calibration
n_calibTrials_Ep = 3;

%% mental parameters
% define number of pairs to solve for each level of difficulty
n_to_reach = mental_N_answersPerLevel(n_E_levels);

% calibration: calibrate the maximal duration required for the
% top effort
n_calibMax = n_to_reach.(['E_level_',num2str(n_E_levels)]);
n_calibTrials_Em = 3;
% pre and post-task calibration (lower number of trials)
n_calibTrials_Em_bis = 3;

% learning
% perform 2 learning sessions, one with instructions and then one without
% (left/right) vs (odd/even) and (lower/higher than 5) - mapping indicated the first time)
% need to remind the mapping the second time
learning_cols = {'col1','col2','all'};
n_learningColours = length(learning_cols);
learning_instructions = {'fullInstructions','noInstructions'}; %,'partialInstructions'
n_learningInstructions = length(learning_instructions);
% initial learning: careful to enter a pair number here
n_maxLearning.learning_withInstructions = 8;
n_maxLearning.learning_withoutInstructions = 8;
warning('left few training trials for Arthur, but need to increase for actual subjects');

%% physical preparation
    %% physical MVC
[initial_MVC, onsets_initial_MVC] = physical_effort_MVC(scr, dq, n_MVC_repeat, calibTimes_Ep);
MVC = nanmax(initial_MVC.MVC); % expressed in Voltage
% 
%     %% learning physical
% [learningPerfSummary_Ep, learningOnsets_Ep] = physical_learning(scr, stim, dq, n_E_levels, Ep_time_levels,...
%     F_threshold, F_tolerance, MVC,...
%     n_learningForceRepeats, learningTimes_Ep);
% 
%     %% training physical
% for iTrainingCondition = 1:n_trainingConditions
%     trainingCond = trainingConditions{iTrainingCondition};
%     
%     % define parameters for the training
%     % reward/punishment and effort levels
%     [trainingChoiceOptions_Ep_tmp, n_trainingTrials_Ep_tmp] = training_options(trainingCond, n_R_levels, n_E_levels);
%     
    % start with reward training alone
    Ep_vars.MVC = MVC;
    Ep_vars.dq = dq;
    Ep_vars.Ep_time_levels = Ep_time_levels;
    Ep_vars.F_threshold = F_threshold;
    Ep_vars.F_tolerance = F_tolerance;
%     [onsets_Ep_training.(trainingCond)] = choice_and_perf_trainingInstructions(scr, trainingCond, trainingTimes_Ep.instructions);
%     [trainingSummary_Ep.(trainingCond)] = choice_and_perf(scr, stim, key, 'physical', Ep_vars, R_money,...
%         trainingCond, n_trainingTrials_Ep_tmp, trainingChoiceOptions_Ep_tmp, trainingTimes_Ep,...
%         results_folder, file_nm_training_Ep);
% end % learning condition loop
% 
% DrawFormattedText(window,'Bravo! Votre entra�nement physique est termin�.',...
%     'center','center',scr.colours.black, scr.wrapat);
% [~,onsets.EndTrainingMsg] = Screen('Flip',window); % display the cross on screen
% WaitSecs(trainingTimes_Ep.trainingEnd);

%% mental preparation
%     %% learning mental
mentalE_prm_learning_and_calib = mental_effort_parameters(iSubject);
mentalE_prm_learning_and_calib.startAngle = 0; % for learning always start at zero
% % no time limit for each trial: as long as needed until learning is
% % ok
% learning_time_limit = false;
% % extract numbers to use for each learning phase
% [numberVector_learning] = mental_numbers(n_learningColours*n_learningInstructions);
% jLearningSession = 0;
% for iCol = 1:n_learningColours
%     curr_learning_col = learning_cols{iCol};
%     for iLearning_Instructions = 1:n_learningInstructions
%         curr_learning_instructions = learning_instructions{iLearning_Instructions};
%         
%         jLearningSession = jLearningSession + 1;
%         learning_sess_nm = ['learning_session',num2str(jLearningSession)];
%         % display instructions for the current learning type
%         [onsets.endLearningInstructions.(learning_sess_nm).(curr_learning_col).(curr_learning_instructions)] = mental_learning(scr,...
%             curr_learning_col, curr_learning_instructions, mentalE_prm_learning_and_calib);
%         
%         % perform the learning
%         [learningPerfSummary_Em.(learning_sess_nm).(curr_learning_col).(curr_learning_instructions)] = mental_effort_perf(scr, stim, key,...
%             numberVector_learning(jLearningSession,:),...
%             mentalE_prm_learning_and_calib, n_maxLearning.learning_withInstructions,...
%             curr_learning_col, curr_learning_instructions, learning_time_limit);
%     end % learning instructions loop
% end % learning colour loop

    %% calibration mental
% extract numbers to use for each calibration trial
[numberVector_calib] = mental_numbers(n_calibTrials_Em);

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

%     %% training mental
% trainingTimes_Em.max_effort = t_min_calib*trainingTimes_Em.t_min_scalingFactor; % allow more time then min performance
% for iTrainingCondition = 1:n_trainingConditions
%     trainingCond = trainingConditions{iTrainingCondition};
%     
%     % define parameters for the training
%     % reward/punishment and effort levels
%     [trainingChoiceOptions_Em_tmp, n_trainingTrials_Em_tmp] = training_options(trainingCond, n_R_levels, n_E_levels);
%     
%     % start with reward training alone
%     Em_vars.i_sub = iSubject;
%     Em_vars.n_to_reach = n_to_reach;
%     [onsets_Em_training.(trainingCond)] = choice_and_perf_trainingInstructions(scr, trainingCond, trainingTimes_Em.instructions);
%     [trainingSummary_Em.(trainingCond)] = choice_and_perf(scr, stim, key, 'mental', Em_vars, R_money,...
%         trainingCond, n_trainingTrials_Em_tmp, trainingChoiceOptions_Em_tmp, trainingTimes_Em,...
%         results_folder, file_nm_training_Em);
% end % learning condition loop

%% actual task
% for mental effort timing
taskTimes_Em.max_effort     = t_min_calib*taskTimes_Em.t_min_scalingFactor; % allow more time then min performance

% instruction that main task will start soon
DrawFormattedText(window,...
    'L''exp�rimentateur va bient�t d�marrer la t�che.',...
    'center', yScreenCenter*(5/3), scr.colours.black, scr.wrapat);
[~, onsets.taskWillStart] = Screen(window, 'Flip');
disp('Please press space.');
[~, ~, keyCode] = KbCheck();
while(keyCode(key.space) ~= 1)
    % wait until the key has been pressed
    [~, ~, keyCode] = KbCheck();
end

% each block 1) MVC 2) task 3) MVC
for iSession = 1:n_sessions
    % define session number
    if ismember(iSession,[1,2])
        session_nm = 'session1';
    elseif ismember(iSession,[3,4])
        session_nm = 'session2';
    end
    
    % define trials (here is where you might want to replace the function
    % by a fixed matrix design
    [choiceOptions_tmp] = choice_option_design(n_R_levels, n_E_levels, punishment_yn, n_trialsPerSession);
    
    if ( (mod(iSubject,2) == 0) && ismember(iSession,[1,3]) ) ||...
            ( (mod(iSubject,2) ~= 0) && ismember(iSession,[2,4]) )% physical task
        
        % pre-task MVC
        [MVC_preTask.(session_nm), onsets_preTask_MVC.physical.(session_nm)] = physical_effort_MVC(scr, dq, n_MVC_repeat, calibTimes_Ep);

        % task
        [perfSummary.physical.(session_nm)] = choice_and_perf(scr, stim, key,...
            'physical', Ep_vars, R_money,...
            'mainTask', n_trialsPerSession, choiceOptions_tmp, taskTimes_Ep,...
            results_folder, [file_nm,'_physical_',session_nm]);
        choiceOptions.physical.(session_nm) = choiceOptions_tmp;

        % post-task MVC
        [MVC_postTask.(session_nm), onsets_postTask_MVC.physical.(session_nm)] = physical_effort_MVC(scr, dq, n_MVC_repeat, calibTimes_Ep);
    elseif ( (mod(iSubject,2) == 0) && ismember(iSession,[2,4]) ) ||...
            ( (mod(iSubject,2) ~= 0) && ismember(iSession,[1,3]) )% mental task
        % pre-task max perf
        % extract numbers to use for each calibration trial
        [numberVector_calib_tmp] = mental_numbers(n_calibTrials_Em_bis);
        
        % alternatively, use fixed number of correct answers to provide for each effort
        % level
        % repeat calibration until the subject performance is better
        % than the requested time threshold
        [t_min_calib_preTask.(session_nm), calibSessionSummary_preTask.(session_nm), calibSuccess_preTask.(session_nm)] = mental_calibTime(scr, stim, key,...
            numberVector_calib_tmp, mentalE_prm_learning_and_calib, n_calibTrials_Em_bis, n_calibMax, calibTimes_Em);
        % task
        [perfSummary.mental.(session_nm)] = choice_and_perf(scr, stim, key,...
            'mental', Em_vars, R_money,...
            'mainTask', n_trialsPerSession, choiceOptions_tmp, taskTimes_Em,...
            results_folder, [file_nm,'_mental_',session_nm]);
        choiceOptions.mental.(session_nm) = choiceOptions_tmp;

        % post-task max perf
        [numberVector_calib_tmp_bis] = mental_numbers(n_calibTrials_Em_bis);
        
        % alternatively, use fixed number of correct answers to provide for each effort
        % level
        % repeat calibration until the subject performance is better
        % than the requested time threshold
        [t_min_calib_postTask.(session_nm), calibSessionSummary_postTask.(session_nm), calibSuccess_postTask.(session_nm)] = mental_calibTime(scr, stim, key,...
            numberVector_calib_tmp_bis, mentalE_prm_learning_and_calib, n_calibTrials_Em_bis, n_calibMax, calibTimes_Em);
    end % nature of the task
    % display feedback for the current session
    DrawFormattedText(window,...
        ['F�licitations! Cette session est maintenant termin�e.',...
        'Vous avez obtenu: ',num2str(perfSummary.totalGain(nTrials)),' chf au cours de cette session.']);
    
end % session loop

%% save the data
% global variables
all.physical.global.MVC = MVC;
all.mental.global.t_min_calib = t_min_calib;
% learning performance
all.physical.learning = learningPerfSummary_Ep;
all.mental.learning = learningPerfSummary_Em;
% training performance
all.physical.training.performance = trainingSummary_Ep;
all.physical.training.instructions = onsets_Ep_training;
all.mental.training.performance = trainingSummary_Em;
all.mental.training.instructions = onsets_Em_training;
% best performance pre/post-task
all.physical.preTask.MVC = MVC_preTask;
all.physical.preTask.onsets_MVC = onsets_preTask_MVC.physical;
all.physical.postTask.MVC = MVC_postTask;
all.physical.postTask.onsets_MVC = onsets_postTask_MVC.physical;
all.mental.preTask.t_min_calib = t_min_calib_preTask;
all.mental.postTask.t_min_calib = t_min_calib_postTask;
all.mental.preTask.calibSessionSummary = calibSessionSummary_preTask;
all.mental.postTask.calibSessionSummary = calibSessionSummary_postTask;
all.mental.preTask.calibSuccess = calibSuccess_preTask;
all.mental.postTask.calibSuccess = calibSuccess_postTask;
% actual performance in the main task sessions
% record physical main task data
all.physical.session1.choiceOptions = choiceOptions.physical.session1;
all.physical.session1.perfSummary = perfSummary.physical.session1;
all.physical.session2.choiceOptions = choiceOptions.physical.session2;
all.physical.session2.perfSummary = perfSummary.physical.session2;
% record mental main task data
all.mental.session1.choiceOptions = choiceOptions.mental.session1;
all.mental.session1.perfSummary = perfSummary.mental.session1;
all.mental.session2.choiceOptions = choiceOptions.mental.session2;
all.mental.session2.perfSummary = perfSummary.mental.session2;

% actually save the data
save([results_folder, file_nm,'.mat'],'all');

%% close PTB
sca;
