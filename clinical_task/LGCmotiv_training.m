function [] = LGCmotiv_training(scr, iSubject, subResultFolder, n_visit, p_or_m, keys, n_buttonsChoice, taskToPerform)
% [] = LGCmotiv_training(window, visit_nm, n_buttonsChoice, taskToPerform)
% LGCmotiv_training will perform the LGC motivation training for visits 1,
% 2 and 3 of the clinical trial depending on the inputs.
%
% INPUTS
% scr: structure with information on window where stimuli will be displayed
%
% n_visit: visit number (1/2/3)
%
% n_buttonsChoice: number of buttons to use for the choice
%
% taskToPerform: structure indicating which subtasks have to be performed
% with on/off criterium
%
% Developed by Nicolas Clairis and Arthur Barakat - 2020-2022

% check visit number ok
if ~ismember(n_visit,[1,2,3])
    error('problem with visit number');
end
visit_nm = ['_v',num2str(n_visit)];

%% define output file names
file_nm_training_Em = ['training_data_Em_CID',num2str(iSubject),visit_nm];
file_nm_training_Ep = ['training_data_Ep_CID',num2str(iSubject),visit_nm];
file_nm = ['training_data_CID',num2str(iSubject),visit_nm];
fullFileNm = [subResultFolder, file_nm,'.mat'];
if exist(fullFileNm,'file')
    error([fullFileNm,' file already exists. If the training script was launched and crashed, ',...
        ' please consider renaming the files already saved to avoid losing data.']);
end
subjectCodeName = ['CID',iSubject];
Ep_visitCalib_filenm = [subResultFolder,subjectCodeName,'_physicalCalib',visit_nm,'.mat'];
Em_visitCalib_filenm = [subResultFolder,subjectCodeName,'_mentalCalib',visit_nm,'.mat'];
% baseline max
Ep_calib_filenm = [subResultFolder,subjectCodeName,'_physicalCalib_v1.mat'];
Em_calib_filenm = [subResultFolder,subjectCodeName,'_mentalCalib_v1.mat'];
% IP file name
file_nm_IP = ['delta_IP_CID',num2str(iSubject),visit_nm];

%% general parameters
window = scr.window;
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
        trainingRP_P_or_R = 'RP';
    case 'no'
        trainingRP_P_or_R = 'R';
end

% define number of training trials (ie when choice + perf)
n_trainingTrials = 4;
% define number of trials per staircase procedure
n_trialsPerSession = 5;
% load timings for each phase of the experiment
[trainingTimes_Em, calibTimes_Em, learningTimes_Em, taskTimes_Em] = timings_definition({trainingRP_P_or_R}, n_trialsPerSession, n_trainingTrials, 'mental');
[trainingTimes_Ep, calibTimes_Ep, learningTimes_Ep, taskTimes_Ep] = timings_definition({trainingRP_P_or_R}, n_trialsPerSession, n_trainingTrials, 'physical');


n_IPsessions = 2;

% time for end of session
t_endSession = trainingTimes_Ep.endSession;

%% physical parameters
if strcmp(taskToPerform.physical.calib,'on') ||...
        strcmp(taskToPerform.physical.learning,'on') ||...
        strcmp(taskToPerform.physical.training,'on') ||...
        strcmp(taskToPerform.physical.task,'on')
    % define relevant keys and dynamometer module
    [dq] = grip_initialization();
    % define conditions
    F_threshold = 55; % force should be maintained above this threshold (expressed in % of MVC)
    F_tolerance = 2.5; % tolerance allowed around the threshold (expressed in % of MVC)
    % need to define timings for each level of force
    [Ep_time_levels] = physical_effortLevels(n_E_levels);
else
    dq = [];
end

%% run the code twice, each time for p or m conditions
for i_pm = 1:2
    switch p_or_m
        case 'p'
            %% physical MVC (calibrate the Fmax for the whole experiment)
            if strcmp(taskToPerform.physical.calib,'on')
                n_MVC_repeat = 3;
                [MVC_tmp, onsets_MVC] = physical_effort_MVC(scr, stim, dq, n_MVC_repeat, calibTimes_Ep, 'MVC', keys);
                MVC = mean(MVC_tmp.MVC); % expressed in Voltage
                save(Ep_visitCalib_filenm,'MVC');
            end
            if n_visit == 1
                MVC_v1 = MVC;
            elseif n_visit > 1
                MVC_v1 = getfield(load(Ep_calib_filenm,'MVC'),'MVC');
            end
            
            %% learning physical (learn each level of force)
            if strcmp(taskToPerform.physical.learning,'on')
                % introduce physical learning
                showTitlesInstruction(scr,stim,'learning',p_or_m, keys);

                n_Ep_learningForceRepeats = 5; % number of learning repetitions for each level of difficulty (= each level of force)
                % perform physical learning
                [learningPerfSummary_Ep, learningOnsets_Ep] = physical_learning(scr, stim, dq, n_E_levels, Ep_time_levels,...
                    F_threshold, F_tolerance, MVC_v1,...
                    n_Ep_learningForceRepeats, learningTimes_Ep);
                
                %% temporary save of data
                save([subResultFolder, file_nm,'.mat']);
            end % physical learning
            
            %% training physical (choice + effort)
            if strcmp(taskToPerform.physical.training,'on')
                % introduce physical training
                showTitlesInstruction(scr,stim,'training',p_or_m, keys);
                
                % define parameters for the training
                Ep_vars_training.MVC = MVC_v1;
                Ep_vars_training.dq = dq;
                Ep_vars_training.Ep_time_levels = Ep_time_levels;
                Ep_vars_training.F_threshold = F_threshold;
                Ep_vars_training.F_tolerance = F_tolerance;
                Ep_vars_training.timeRemainingEndTrial_ONOFF = 0;
                % mapping between reward levels and fake monetary amounts
                R_money = R_amounts(n_R_levels, punishment_yn);
                % reward/punishment and effort levels
                [trainingChoiceOptions_Ep_tmp] = training_options(trainingRP_P_or_R, n_R_levels, n_E_levels, R_money, n_trainingTrials);
                % perform physical training in 2 phases:
                % 1) with confidence mapping;
                % 2) without confidence mapping
                trainingConfConditions = {'RP_withConfMapping','RP_withoutConfMapping'};
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
                    if strcmp(trainingConfCond,'RP_withoutConfMapping') || n_buttonsChoice == 2
                        confidenceChoiceDisplay = false;
                    elseif strcmp(trainingConfCond,'RP_withConfMapping') && n_buttonsChoice == 4
                        confidenceChoiceDisplay = true;
                    end
                    % extract trials to use
                    % training instruction
                    [onsets_Ep_training.(['session',num2str(iTrainingCondition),'_',trainingConfCond])] = choice_and_perf_trainingInstructions(scr, stim, trainingConfCond, trainingTimes_Ep.instructions);
                    % perform the training
                    [trainingSummary_Ep.(['session',num2str(iTrainingCondition),'_',trainingConfCond])] = choice_and_perf(scr, stim, keys, 'physical', Ep_vars_training,...
                        trainingRP_P_or_R, n_trialsPerTrainingCondition, trainingChoiceOptions_Ep_tmp, confidenceChoiceDisplay,...
                        trainingTimes_Ep,...
                        subResultFolder, file_nm_training_Ep);
                end % learning condition loop
                
                DrawFormattedText(window, stim.training.Ep.endMsg.text,...
                    'center','center',stim.training.Ep.endMsg.colour, scr.wrapat);
                [~,onsets.EndTrainingMsg] = Screen('Flip',window); % display the cross on screen
                WaitSecs(trainingTimes_Ep.trainingEnd);
                
                %% temporary save of data
                save([subResultFolder, file_nm,'.mat']);
            end % training physical
            
        case 'm'
            %% learning mental: 0-back and 2-back as a calibration
            if strcmp(taskToPerform.mental.learning_1,'on')
                % display instructions for learning
                showTitlesInstruction(scr,stim,'learning',p_or_m, keys);
                
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
                n_maxLearning.learning_2back = 6;
                mentalE_prm_learning = mental_effort_parameters();
                mentalE_prm_learning.startAngle = 0; % for learning always start at zero
                % no time limit for each trial: as long as needed until learning is ok
                learning_useOfTimeLimit = false;
                
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
                        curr_learning_instructions, mentalE_prm_learning); % inform about answer to give
                    
                    % perform the learning
                    if iLearning_Instructions == 1
                        n_maxToReachLearning_tmp = n_maxLearning.learning_withInstructions;
                    elseif iLearning_Instructions == 2
                        n_maxToReachLearning_tmp = n_maxLearning.learning_withoutInstructions;
                    end
                    [learningPerfSummary_Em.(learning_sess_nm).(curr_learning_instructions)] = mental_effort_perf(scr, stim, keys,...
                        numberVector_learning(jLearningSession,:),...
                        mentalE_prm_learning, n_maxToReachLearning_tmp,...
                        curr_learning_instructions, learning_useOfTimeLimit, [], learning_errorLimits);
                    
                    
                    % for experimenter display how many trials have been performed
                    disp(['Mental 0-back learning (1) trial ',num2str(jMentalLearningTrial),'/',num2str(nMentalLearning_totalTrials),' done']);
                end % learning instructions loop
                
                %% learning (1) of the 2-back without temporal pressure
                n_learning1_2back = 2;
                mentalE_prm_learning1_2back = mental_effort_parameters();
                % always start at zero
                mentalE_prm_learning1_2back.startAngle = 0;
                % Nback version
                Nback_str = num2str(mentalE_prm_learning1_2back.Nback);
                learningVersion = ['learning_Nback',Nback_str];
                % time limits
                learning1_2back_useOfTimeLimit = false;
                
                % define conditions for the learning
                [numberVector_learning1_2back] = mental_numbers(n_learning1_2back);
                % error handling for learning
                learning1_2back_errorLimits.useOfErrorThreshold = false;
                learning1_2back_errorLimits.useOfErrorMapping = false;
                
                % perform the learning session
                [onsets.endLearningInstructions.learning1_2back_session] = mental_learningInstructions(scr, stim,...
                    learningVersion, mentalE_prm_learning1_2back); % inform about 2-back
                 for iLearning_2backTrial = 1:n_learning1_2back
                    mentalE_learning1_2backPerfSummary_tmp = mental_effort_perf_Nback(scr, stim, keys,...
                        numberVector_learning1_2back(iLearning_2backTrial,:),...
                        mentalE_prm_learning1_2back, n_maxLearning.learning_2back,...
                        'noInstructions', learning1_2back_useOfTimeLimit, [],...
                        learning1_2back_errorLimits, [], []);
                    learningPerfSummary_Em.learning1_2back.(['trial_',num2str(iLearning_2backTrial)]) = mentalE_learning1_2backPerfSummary_tmp;
                    
                    % small break between each answer
                    DrawFormattedText(window, stim.training.Em.endTrialMsg_good.text,'center',yScreenCenter/2,white);
                    DrawFormattedText(window,stim.training.Em.endTrialMsg_bis.text,'center','center',white);
                    [~,~,onsets.timeLearning2backFbk.(['trial_',num2str(iLearning_2backTrial)])] = Screen(window,'Flip');
                    WaitSecs(learningTimes_Em.learning_rest);
                    disp(['Mental learning (1) 2-back trial ',num2str(iLearning_2backTrial),'/',num2str(n_learning1_2back),' done']);
                end % trial loop
                
                %% learning (1) by repeating the calibration many times before actual calibration
                % define number of trials to perform
                n_learning1calibLikeTrials = 30;
                
                mentalE_prm_learning1calibLike = mental_effort_parameters();
                % always start at zero
                mentalE_prm_learning1calibLike.startAngle = 0;
                % Nback version
                Nback_str = num2str(mentalE_prm_learning1calibLike.Nback);
                learningVersion = ['learning_Nback',Nback_str,'_bis'];
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
                    mentalE_learning1calibLikePerfSummary_tmp = mental_effort_perf_Nback(scr, stim, keys,...
                        numberVector_learning1calibLike(iLearning1Trial,:),...
                        mentalE_prm_learning1calibLike, n_maxToReachForCalib,...
                        'noInstructions', learning1calibLike_useOfTimeLimit, learning1calibLike_timeLimit,...
                        learning1calibLike_errorLimits, nMaxReachedUntilNowLearning, n_Em_learning1calibLike_MinToReach);
                    learningPerfSummary_Em.learning1calibLike.(['trial_',num2str(iLearning1Trial)]) = mentalE_learning1calibLikePerfSummary_tmp;
                    n_maxReachedDuringLearning(iLearning1Trial) = mentalE_learning1calibLikePerfSummary_tmp.n_correctAnswersForDisplay;
                    
                    % extract new best performance
                    nMaxReachedUntilNowLearning = max(nMaxReachedUntilNowLearning, n_maxReachedDuringLearning(iLearning1Trial));
                    % small break between each answer, when they press a
                    % button break is over
                    if n_maxReachedDuringLearning(iLearning1Trial) >= n_maxToReachForCalib % congrats only when perf is good
                        DrawFormattedText(window, stim.training.Em.endTrialMsg_good.text,'center',yScreenCenter/2,white);
                    else
                        DrawFormattedText(window, stim.training.Em.endTrialMsg_bad.text,'center',yScreenCenter/2,white);
                    end
                    DrawFormattedText(window,stim.training.Em.endTrialMsg_bis.text,'center','center',white);
                    DrawFormattedText(window, stim.pressWhenReady.text,...
                        stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour, wrapat);
                    [~,~,onsets.timeLearningFbk.(['trial_',num2str(iLearning1Trial)])] = Screen(window,'Flip');
                    KbQueueWait(0,3); % wait for button press and button release before moving on
                    disp(['Mental learning (1) calibration-like trial ',num2str(iLearning1Trial),'/',num2str(n_learning1calibLikeTrials),' done']);
                end % trial loop
                
                learning1done = 0;
                n_learning1bonusTrialsToLearn = 5; % how many trials to use as a learning penalty
                n_lastTrialsToCheck = 5; % how many trials to check
                n_trialsCorrectThreshold = 4; % if less (<) than this number of trials was correct in the n_lastTrialsToCheck trials, redo more trials
                jLearning1Trial = n_learning1calibLikeTrials;
                iBlockRepeats = 0;
                while learning1done == 0
                    learningPerf_lastTrials = zeros(1,n_lastTrialsToCheck);
                    for iLastTrial = 1:n_lastTrialsToCheck
                        learningPerf_lastTrials(iLastTrial) = n_maxReachedDuringLearning(end+1-iLastTrial) >= n_Em_learning1calibLike_MinToReach;
                    end
                    if n_learning1calibLikeTrials >= n_lastTrialsToCheck &&...
                            ( sum(learningPerf_lastTrials) < n_trialsCorrectThreshold)
                        disp(['performance was too low in one of the last trials. We will redo ',...
                            num2str(n_learning1bonusTrialsToLearn),' more trials to compensate.']);
                        iBlockRepeats = iBlockRepeats + 1;
                        disp(['starting now the ',num2str(iBlockRepeats),' bonus block of the mental learning.']);
                        [numberVector_learning1_bonus] = mental_numbers(n_learning1bonusTrialsToLearn);
                        for iLearning1Trial_bonus = 1:n_learning1bonusTrialsToLearn
                            jLearning1Trial = jLearning1Trial + 1;
                            mentalE_learning1calibLikePerfSummary_tmp = mental_effort_perf_Nback(scr, stim, keys,...
                                numberVector_learning1_bonus(iLearning1Trial_bonus,:),...
                                mentalE_prm_learning1calibLike, n_maxToReachForCalib,...
                                'noInstructions', learning1calibLike_useOfTimeLimit, learning1calibLike_timeLimit,...
                                learning1calibLike_errorLimits, nMaxReachedUntilNowLearning, n_Em_learning1calibLike_MinToReach);
                            learningPerfSummary_Em.learning1calibLike.(['trial_',num2str(jLearning1Trial)]) = mentalE_learning1calibLikePerfSummary_tmp;
                            n_maxReachedDuringLearning(jLearning1Trial) = mentalE_learning1calibLikePerfSummary_tmp.n_correctAnswersForDisplay;

                            % extract new best performance
                            nMaxReachedUntilNowLearning = max(nMaxReachedUntilNowLearning, n_maxReachedDuringLearning(jLearning1Trial));
                            % small break between each answer, when they press a
                            % button break is over
                            if n_maxReachedDuringLearning(jLearning1Trial) >= n_maxToReachForCalib % congrats only when perf is good
                                DrawFormattedText(window, stim.training.Em.endTrialMsg_good.text,'center',yScreenCenter/2,white);
                            else
                                DrawFormattedText(window, stim.training.Em.endTrialMsg_bad.text,'center',yScreenCenter/2,white);
                            end
                            DrawFormattedText(window,stim.training.Em.endTrialMsg_bis.text,'center','center',white);
                            DrawFormattedText(window, stim.pressWhenReady.text,...
                                stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour, wrapat);
                            [~,~,onsets.timeLearningFbk.(['trial_',num2str(jLearning1Trial)])] = Screen(window,'Flip');
                            KbQueueWait(0,3); % wait for button press and button release before moving on
                            disp(['Mental learning (1) calibration-like BONUS trial ',num2str(iLearning1Trial_bonus),'/',num2str(n_learning1bonusTrialsToLearn),' done']);
                        end % trial loop
                    else
                        learning1done = 1;
                    end
                end
                
                %% temporary save of data
                save([subResultFolder, file_nm,'.mat']);
            end % mental learning (1)
            
            %% calibration mental
            if strcmp(taskToPerform.mental.calib,'on')
                % number of calibration trials
                n_calibTrials_Em = 3;
                
                % mental calibration parameters
                mentalE_prm_calib = mental_effort_parameters();
                mentalE_prm_calib.calib_or_maxPerf = 'calib';
                mentalE_prm_calib.startAngle = 0; % for learning always start at zero
                n_calibMax = mentalE_prm_calib.n_maxToReachCalib;
                % extract numbers to use for each calibration trial
                %     [numberVector_calib] = mental_numbers(n_calibTrials_Em);
                [numberVector_calib] = mental_calibNumberVector(n_calibTrials_Em, n_calibMax);
                % perform the calibration
                [NMP, calib_summary] = mental_calibNumbers(scr, stim, keys,...
                    numberVector_calib, mentalE_prm_calib, n_calibTrials_Em, calibTimes_Em, langage);
                calibSummary.calibSummary = calib_summary;
                calibSummary.n_mental_max_perTrial = NMP;
                % record number of maximal performance (NMP)
                save(Em_visitCalib_filenm,'NMP');
            elseif strcmp(taskToPerform.mental.calib,'off') &&...
                    ( strcmp(taskToPerform.mental.learning_1,'on') ||...
                    strcmp(taskToPerform.mental.learning_2,'on') ||...
                    strcmp(taskToPerform.mental.training,'on') ||...
                    strcmp(taskToPerform.mental.task,'on') )
                NMP = getfield(load(Em_calib_filenm,'NMP'),'NMP');
            end % calibration
            if n_visit == 1
                NMP_v1 = NMP;
            elseif n_visit > 1
                NMP_v1 = getfield(load(Em_calib_filenm,'NMP'),'NMP');
            end
            
            %% learning (2) for each difficulty level
            if strcmp(taskToPerform.mental.learning_2,'on')
                % introduce mental learning
                showTitlesInstruction(scr,stim,'learning','m', keys);
                
                % define all difficulty levels based on calibration
                [n_to_reach] = mental_N_answersPerLevel(n_E_levels, NMP_v1);
                % define number of learning trials
                n_Em_learningForceRepeats = 5; % number of learning repetitions for each level of difficulty (= each level of force)
                % timings
                Em_learningTimings = learningTimes_Em;
                Em_learningTimings.time_limit = false;
                % perform all the difficulty levels
                [learning2PerfSummary_Em, onsets] = mental_learning(scr, stim, keys, n_E_levels, n_to_reach, n_Em_learningForceRepeats, Em_learningTimings);
                
                %% temporary save of data
                save([subResultFolder, file_nm,'.mat']);
            end % learning (2)
            
            %% training mental
            if strcmp(taskToPerform.mental.training,'on')
                % show title before instructions
                showTitlesInstruction(scr, stim, 'training', 'm', keys);

                % define parameters for the training
                [Em_vars_training.n_to_reach] = mental_N_answersPerLevel(n_E_levels, NMP_v1);
                Em_vars_training.errorLimits.useOfErrorMapping = false;
                Em_vars_training.errorLimits.useOfErrorThreshold = false;
                Em_vars_training.timeRemainingEndTrial_ONOFF = 0;
                % mapping between reward levels and fake monetary amounts
                R_money = R_amounts(n_R_levels, punishment_yn);
                % reward/punishment and effort levels
                [trainingChoiceOptions_Em_tmp] = training_options(trainingRP_P_or_R, n_R_levels, n_E_levels, R_money, n_trainingTrials);
                % perform physical training in 2 phases:
                % 1) with confidence mapping;
                % 2) without confidence mapping
                trainingConfConditions = {'RP_withConfMapping','RP_withoutConfMapping'};
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
                    if strcmp(trainingConfCond,'RP_withoutConfMapping') || n_buttonsChoice == 2
                        confidenceChoiceDisplay = false;
                    elseif strcmp(trainingConfCond,'RP_withConfMapping') && n_buttonsChoice == 4
                        confidenceChoiceDisplay = true;
                    end
                    % extract trials to use
                    trainingTrials_idx = (1:n_trialsPerTrainingCondition) + n_trialsPerTrainingCondition*(iTrainingCondition - 1);
                    % training instruction
                    [onsets_Em_training.(['session',num2str(iTrainingCondition),'_',trainingConfCond])] = choice_and_perf_trainingInstructions(scr, stim, trainingConfCond, trainingTimes_Em.instructions);
                    % select effort level
                    [trainingSummary_Em.(trainingConfCond)] = choice_and_perf(scr, stim, keys, 'mental', Em_vars_training,...
                        trainingRP_P_or_R, n_trialsPerTrainingCondition, trainingChoiceOptions_Em_tmp, confidenceChoiceDisplay,...
                        trainingTimes_Em,...
                        subResultFolder, file_nm_training_Em);
                end % training condition loop
                
                %% temporary save of data
                save([subResultFolder, file_nm,'.mat']);
            end % mental training
    end % physical/mental loop
    
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
    
    % Baseline Reward (CHF)
    baselineR = 0.5;
    baselineP = 0.5;
    
    % Total amount of money to be given
    totalGain = 0;
    
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
            KbQueueWait(0,3);
        end
    end % loop over forced reading/manual pass loop
    
    % for physical effort
    if strcmp(taskToPerform.physical.task,'on')
        Ep_vars.MVC = MVC_v1;
        Ep_vars.dq = dq;
        Ep_vars.Ep_time_levels = Ep_time_levels;
        Ep_vars.F_threshold = F_threshold;
        Ep_vars.F_tolerance = F_tolerance;
        Ep_vars.timeRemainingEndTrial_ONOFF = 0;
    end

    % keep track of the current session for mental and physical (only important for saving)
    iMental = 1;
    iPhysical = 1;

    % number of runs
    for iSession = 1:n_IPsessions
        % session gains
        sessionFinalGain = 0;

        % define if they will see a punishment or reward trial
        % always start with a reward, they have to win money first
        session_nm = ['session_nb',num2str(iSession)];
        R_or_P = 'R'; % compute IP only for rewards

        for i_pm = 1:2

            % perform physical or mental IP according to inputs
            switch p_or_m
                case 'p'
                    if strcmp(taskToPerform.physical.task,'on')
                        % instructions
                        showTitlesInstruction(scr,stim,'task',p_or_m, keys);

                        % run physical task
                        perf_Ep_IP_tmp = choice_and_perf_staircase(scr, stim, keys,...
                            'physical', Ep_vars,...
                            R_or_P,E_right(iEffortLevel),E_left(iEffortLevel), n_trialsPerSession, taskTimes_Ep,...
                            subResultFolder, [file_nm,'_physical_session_nb',session_nm,'_effort_lvl',num2str(iEffortLevel)]);
                        perfSummary.physical.(['session_nb',num2str(iPhysical)]).(['Effort_lvl',(num2str(iEffortLevel))]) = perf_Ep_IP_tmp;
                        sessionFinalGain = sessionFinalGain + perf_Ep_IP_tmp.totalGain(end);
                        iPhysical = iPhysical + 1;
                    end

                case 'm'
                    if strcmp(taskToPerform.mental.task,'on')
                        % instructions
                        showTitlesInstruction(scr,stim,'task',p_or_m, keys);

                        mentalE_prm_instruDisplay = mental_effort_parameters();
                        % Nback version
                        [Em_vars.n_to_reach] = mental_N_answersPerLevel(n_E_levels, NMP_v1);
                        % run mental task
                        % for actual task: no display of mapping but consider 3
                        % errors as a trial failure
                        Em_vars.errorLimits.useOfErrorMapping = false;
                        Em_vars.errorLimits.useOfErrorThreshold = false;
                        Em_vars.timeRemainingEndTrial_ONOFF = 0;
                        perf_Em_IP_tmp = choice_and_perf_staircase(scr, stim, keys,...
                            'mental', Em_vars,...
                            R_or_P,E_right((iEffortLevel)),E_left(iEffortLevel), n_trialsPerSession, taskTimes_Em,...
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
        switch langage
            case 'fr'
                DrawFormattedText(window,...
                    'Felicitations! Cette session est maintenant terminee.',...
                    'center', yScreenCenter*(5/3), scr.colours.white, scr.wrapat);
            case 'engl'
                DrawFormattedText(window,...
                    'Congratulations! This session is now completed.',...
                    'center', yScreenCenter*(5/3), scr.colours.white, scr.wrapat);
        end
        Screen(window,'Flip');
        % give break after 4 IP
        WaitSecs(t_endSession);
        totalGain = totalGain + sessionFinalGain;


    end % session loop
    %change m/p for next loop
    
end % mental/physical on?

%% save the data

% calibration
if strcmp(taskToPerform.physical.calib,'on')
    all.physical.MVC = MVC;
    all.physical.onsets_MVC = onsets_MVC;
end
if strcmp(taskToPerform.mental.calib,'on')
    all.mental.NMP = NMP;
    all.mental.calib_summary = calib_summary;
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
    all.physical.learning = learningPerfSummary_Ep;
    all.physical.learning = learningOnsets_Ep;
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
        for iSession = 1:n_IPsessions % for the physical sessions
            % save data in all and reformat it in a specific order
            all.physical.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).perfSummary = perfSummary.physical.(['session_nb',num2str(iSession)]).(['Effort_lvl',(num2str(iEffort))]);
            IP_variables.physicalDeltaIP_perSession(iSession) = all.physical.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).perfSummary.IP - baselineR;
        end
    end
    IP_variables.physicalDeltaIP = mean(IP_variables.physicalDeltaIP_perSession);
end
% record mental main task data
if strcmp(taskToPerform.mental.task,'on')
    for iEffort= 1:nbEffortLvl % for all effort levels
        for iSession = 1:n_IPsessions % for the mental sessions
            % save data in all and reformat it in a specific order
            all.mental.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).perfSummary = perfSummary.mental.(['session_nb',num2str(iSession)]).(['Effort_lvl',(num2str(iEffort))]);
            IP_variables.mentalDeltaIP_perSession(iSession) = all.mental.(['EffortLvl_',num2str(iEffort)]).(['session_nb',num2str(iSession)]).perfSummary.IP - baselineR;
        end
    end
    IP_variables.mentalDeltaIP = mean(IP_variables.mentalDeltaIP_perSession);
end

if strcmp(taskToPerform.physical.task,'on') || strcmp(taskToPerform.mental.task,'on')
    IP_variables.baselineR = baselineR;
    IP_variables.baselineP = baselineP;
    IP_variables.totalGain = totalGain;
end
IP_variables.training.p_or_m = p_or_m;
if strcmp(taskToPerform.physical.task,'on')
    IP_variables.calibration.MVC = MVC;
end
if strcmp(taskToPerform.mental.task,'on')
    IP_variables.calibration.NMP = NMP;
end
% actually save the final data
save([subResultFolder, file_nm,'.mat']);

% save delta_IP and baselineR
save([subResultFolder, file_nm_IP,'.mat'],'IP_variables');

end % function