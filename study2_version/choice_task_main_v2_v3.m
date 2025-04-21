function[] = choice_task_main_v2_v3(scr, iSubject, subResultsFolder, n_visit, iRun, effort_type,...
            keys, n_buttonsChoice, trainingTaskToPerform)
% main task of the experiment:
% 1) fixation cross
% 2) choice between 2 options
% 3) display chosen option
% 4) perform the effort corresponding to the chosen option
%
% The effort can be a physical effort of a fixed intensity but a
% varying duration (in the physical version of the task or a mental effort
% with a fixed maximal duration but a varying number of questions to solve
% according to difficulty (mental version of the task).
%
% developed by Arthur Barakat & Nicolas Clairis - 2020/2021
% Adapted by Nicolas Clairis may 2022 for clinical trial
%
% See also ScreenConfiguration.m, physical_effort_perf.m,
% mental_effort_perf_Nback.m

%% langage to display instructions
langage = 'fr';

%% Define subject ID
session_nm = num2str(iRun);

% Create subjectCodeName which is used as a file saving name
subjectCodeName = strcat('CID',iSubject);

% calibration performance file name
calibPerf_file_nm = [subResultsFolder, subjectCodeName,'_',effort_type,'Calib.mat'];

% file name for main session
file_nm = [subjectCodeName,'_session',session_nm,'_',effort_type,'_task'];
% verify the files do not already exist
if exist([subResultsFolder, file_nm,'.mat'],'file')
    error(['The file name ',file_nm,'.mat already exists.',...
        ' Please relaunch with a new file name or delete the previous data.']);
end

%% fMRI/behavioral version of the task?
IRM = 0; % (0)non-fMRI/(1) fMRI exp

%% include punishment condition?
punishment_yn = 'yes'; % include punishment trials?

%% task parameters
white = scr.colours.white;

% calibration before/after end of each fMRI session
n_MaxPerfTrials = 2;

% actual task
n_R_levels = 4;
n_E_levels = 4;
nTrials = 54;

% extract money amount corresponding to each reward level for the
% computation of the gains
file_nm_IP = ['delta_IP_CID',num2str(iSubject)];
if ~exist([subResultsFolder, file_nm_IP,'.mat'],'file')
    error(['couldn''t find ',file_nm_IP,' file. Please verify the IP was done.']);
end
IPdata = getfield(load([subResultsFolder, file_nm_IP,'.mat'],'IP_variables'),'IP_variables');
[R_money] = R_amounts_IP(n_R_levels, punishment_yn, IPdata, effort_type);

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

% define thresholds for each task
switch effort_type
    case 'physical'
        % initialize the grip
        dq = grip_initialization();
        F_threshold = 55; % force should be maintained above this threshold (expressed in % of MVC)
        F_tolerance = 2.5; % tolerance allowed around the threshold (expressed in % of MVC)
    case 'mental'
        mentalE_prm_calib = mental_effort_parameters();
        mentalE_prm_calib.calib_or_maxPerf = 'maxPerf';
        mentalE_prm_calib.startAngle = 0; % for learning always start at zero
        % no error threshold nor mapping of answers when errors are
        % made
        % calibration and max performances: define the maximal number to reach
        n_calibMax = mentalE_prm_calib.n_maxToReachCalib;
end

% stimulus related variables for the display
[stim] = stim_initialize(scr, n_E_levels, langage);

%% load timings for each phase of the experiment
trainingConditions = {'RP'};
trainingTrials = 3;
[~, calibTimes, ~, taskTimes, mainTimes] = timings_definition(trainingConditions, nTrials, trainingTrials, effort_type);
t_endfMRI = mainTimes.endfMRI;
t_endSession = mainTimes.endSession;

%% load max performance from calibration and define effort difficulty
% levels accordingly
switch effort_type
    case 'mental'
        % load calibration maximal performance
        NmaxCalib = getfield(load(calibPerf_file_nm,'NMP'),'NMP');
        % define number of pairs to solve for each level of difficulty
        n_to_reach = mental_N_answersPerLevel(n_E_levels, NmaxCalib);
    case 'physical'
        MVC = getfield(load(calibPerf_file_nm,'MVC'),'MVC'); % expressed in Voltage
        % define effort levels for each level of difficulty
        [Ep_time_levels] = physical_effortLevels(n_E_levels);
end

%% max perf measurement before start of each fMRI session
switch effort_type
    case 'mental'
        % perform max perf
        [numberVector_initialMaxPerf] = mental_calibNumberVector(n_MaxPerfTrials, n_calibMax);
        [n_initialMaxPerf, initialMaxPerfSessionSummary] = mental_calibNumbers(scr, stim, keys,...
            numberVector_initialMaxPerf, mentalE_prm_calib, n_MaxPerfTrials, calibTimes, langage);
    case 'physical'
        % take an initial MVC measurement (even if it has been done in a
        % previous session, will allow us to keep track of the force level
        % of our participants)
        [initial_MVC, onsets_initial_MVC] = physical_effort_MVC(scr, stim, dq, n_MaxPerfTrials, calibTimes, 'maxPerf');
end

%% launch main task
% define task parameters
switch effort_type
    case 'mental'
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

%% instruction that main task will start soon
DrawFormattedText(window,stim.expWillStartBis.text,...
stim.expWillStartBis.x,stim.expWillStartBis.y,scr.colours.white, scr.wrapat);
[~, onsets.taskWillStart] = Screen(window, 'Flip');
WaitSecs(3);

%% perform choice and performance task
[perfSummary] = choice_and_perf(scr, stim, keys,...
    effort_type, Ep_or_Em_vars,...
    'mainTask', nTrials, choiceOptions, confDispDuringChoice,...
    taskTimes,...
    subResultsFolder, file_nm);

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


%% first save before last MVC
save([subResultsFolder, file_nm,'_messyAllStuff.mat']);

%% Measure maximum power again at the end of each scan
% add instructions
DrawFormattedText(window, stim.postTaskMVCmeasurement.text,...
    stim.postTaskMVCmeasurement.x, stim.postTaskMVCmeasurement.y,...
    stim.postTaskMVCmeasurement.colour, scr.wrapat);
Screen(window,'Flip');

% MVC maximum
switch effort_type
    case 'physical'
        [last_MVC, onsets_last_MVC] = physical_effort_MVC(scr, stim, dq, n_MaxPerfTrials, calibTimes, 'maxPerf');
    case 'mental'
        % extract numbers to use for each calibration trial
        [numberVector_endCalib] = mental_calibNumberVector(n_MaxPerfTrials, n_calibMax);
        % last max performance measurement
        [n_finalMaxPerf, finalMaxPerf_SessionSummary] = mental_calibNumbers(scr, stim, keys,...
            numberVector_endCalib, mentalE_prm_calib, n_MaxPerfTrials, calibTimes, langage);
end

%% display feedback for the current session
DrawFormattedText(window, stim.endSessionMessage.text,...
    stim.endSessionMessage.x, stim.endSessionMessage.y,...
    white, scr.wrapat);
[~,onsets.endSessionFbk] = Screen(window,'Flip');
WaitSecs(t_endSession);

%% Save Data
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
        all.physicalPerf = perfSummary;
    case 'mental'
        if iRun == 0
            all.calibSummary = calibSummary;
        else
            all.start_maxPerf.n_maxPerf = n_initialMaxPerf;
            all.start_maxPerf.sessionSummary = initialMaxPerfSessionSummary;
        end
        % last max perf
        all.end_maxPerf.n_maxPerf = n_finalMaxPerf;
        all.end_maxPerf.SessionSummary = finalMaxPerf_SessionSummary;
        if IRM == 0
            all.mentalE_prm_learning = mentalE_prm_calib;
        end
        all.mentalE_perf = perfSummary;
end
% store timings in all structure
all.calibTimes = calibTimes;
all.taskTimes = taskTimes;
all.onsets = onsets;

% save to finish
save([subResultsFolder, file_nm,'.mat'],'-struct','all');

%% save all variables
save([subResultsFolder, file_nm,'_messyAllStuff.mat']);
end % function