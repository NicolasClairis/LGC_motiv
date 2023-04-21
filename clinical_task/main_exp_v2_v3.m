% main script for visits 2 and 3 of the clinical trial
% 1) re-do the training of the task to remind them how it works
% 2) perform the reward/effort tradeoff incentivized task
% 3) train for the "social" task
% 4) perform the "social" task


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
main_task_folder            = [main_folder, 'clinical_task' filesep];
results_folder              = [main_folder, 'INRB_results' filesep];
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
while isempty(iSubject) || length(iSubject) ~= 3 ||...
        ~ismember(p_or_m_first,{'p','m'}) ||...
        ~ismember(str2double(visit),[1,2,3])
    % repeat until all questions are answered
    info = inputdlg({'Subject CID (XXX)','p/m first?','visit?'});
    [iSubject, p_or_m_first, visit] = info{[1,2,3]};
    %     warning('when real experiment starts, remember to block in french and IRM = 1');
end
% Create subjectCodeName which is used as a file saving name
subjectCodeName = strcat('CID',iSubject);
subResultsFolder = [results_folder, subjectCodeName, filesep,'behavior',filesep];
if ~exist(subResultsFolder,'dir')
    mkdir(subResultsFolder);
end
% convert subject CID into number (only if used to perform actual task)
if ischar(iSubject)
    iSubject = str2double(iSubject);
end

n_visit = str2double(visit);

% other main parameters set by default
IRMdisp = 0; % 0 = display on main screen, 1=display on secondary screen
IRMbuttons = 1; % 0 = use keyboard, 1 = use keypad
IRM = 0; % is exp in MRI (ie need to check the fMRI TTL) then set IRM=1 or out of MRI (IRM=0)?
testing_script = 0; % 0 = ignore problems with computer resources,
% 1=use all computer resources (important for proper timing measurements)
n_buttonsChoice = 4; % set to 4 if you want to measure confidence

% initialize the langage for the experiment
langue = 'f';
switch langue
    case 'f'
        langage = 'fr';
    case 'e'
        langage = 'engl'; % 'fr'/'engl' french or english?
    otherwise
        error('langage not recognised');
end

%% in case of bug (or for debugging) allow to select which scripts to launch
taskToPerform.motivTraining = 'on';
% define which subtasks of the training you want to perform
if strcmp(taskToPerform.motivTraining,'on')
    trainingTaskToPerform.physical.calib = 'on';
    trainingTaskToPerform.physical.learning = 'on';
    trainingTaskToPerform.physical.training = 'on';
    trainingTaskToPerform.physical.task = 'on';
    trainingTaskToPerform.mental.learning_1 = 'on';
    trainingTaskToPerform.mental.calib = 'on';
    trainingTaskToPerform.mental.learning_2 = 'on';
    trainingTaskToPerform.mental.training = 'on';
    trainingTaskToPerform.mental.task = 'on';
end
taskToPerform.motivTask = 'on';
taskToPerform.socialTraining = 'on';
taskToPerform.socialTask = 'on';

%% initialize psychtoolbox once for all scripts
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(IRMdisp, testing_script);
white = scr.colours.white;
black = scr.colours.black;
xSize = xScreenCenter*2;
ySize = yScreenCenter*2;

%% initiliaze stimulus for pressing when ready
pressWhenReadyText = 'Appuyez quand vous etes pret(e) a poursuivre.';
[~,~,textSizePressWhenReady] = DrawFormattedText(window, stim.pressWhenReady.text, 'center', 'center', white);
pressWhenReadyX = xScreenCenter - (textSizePressWhenReady(3) - textSizePressWhenReady(1))/2;
pressWhenReadyY = ySize*(15/16) - (textSizePressWhenReady(4) - textSizePressWhenReady(2))/2;

%% launch check of button presses
keys = relevant_key_definition(IRM, n_buttonsChoice);
keyboard_check_start(keys, IRM);

%% re-training effort tasks
switch taskToPerform.motivTraining
    case 'on'
        % introduce the training
        DrawFormattedText(window,'Re-entraînement','center','center',white);
        DrawFormattedText(window,pressWhenReadyText,pressWhenReadyX,pressWhenReadyY,white,wrapat);
        Screen(window,'Flip');
        KbQueueWait(0,3);

        % launch the training
        LGCmotiv_training(scr, iSubject, subResultsFolder, n_visit, p_or_m_first,...
            keys, n_buttonsChoice, trainingTaskToPerform)

        % finalize the training
        DrawFormattedText(window,'Re-entraînement terminé','center','center',white);
        DrawFormattedText(window,pressWhenReadyText,pressWhenReadyX,pressWhenReadyY,white,wrapat);
        Screen(window,'Flip');
        KbQueueWait(0,3);
end

%% main effort task
switch taskToPerform.motivTask
    case 'on'
        % introduce the task
        DrawFormattedText(window,'Tâche principale','center','center',white);
        DrawFormattedText(window,pressWhenReadyText,pressWhenReadyX,pressWhenReadyY,white,wrapat);
        Screen(window,'Flip');
        KbQueueWait(0,3);

        for iRun = 1:4
            % define task type
            if (ismember(iRun,[1,3]) && strcmp(p_or_m_first,'p')) ||...
                    (ismember(iRun,[2,4]) && strcmp(p_or_m_first,'m'))
                runEffort_type = 'physical';
            elseif (ismember(iRun,[1,3]) && strcmp(p_or_m_first,'m')) ||...
                    (ismember(iRun,[2,4]) && strcmp(p_or_m_first,'p'))
                runEffort_type = 'mental';
            else
                runEffort_type = input(['Unexpected situation: run different from 1-4, ',...
                    'please define if run is ''physical'' or ''mental''.'],'s');
            end % task type

            % perform the task
            choice_task_main_v2_v3(scr, iSubject, subResultsFolder, n_visit, iRun, runEffort_type,...
                keys, n_buttonsChoice, trainingTaskToPerform);
        end % run loop

        % finalize the task
        DrawFormattedText(window,'Bravo! Vous avez terminé la tâche principale.','center','center',white,wrapat);
        DrawFormattedText(window,pressWhenReadyText,pressWhenReadyX,pressWhenReadyY,white,wrapat);
        Screen(window,'Flip');
        KbQueueWait(0,3);
end

%% social training
switch taskToPerform.socialTraining
    case 'on'
        % introduce the training
        DrawFormattedText(window,'Entraînement deuxième tâche','center','center',white);
        DrawFormattedText(window,pressWhenReadyText,pressWhenReadyX,pressWhenReadyY,white,wrapat);
        Screen(window,'Flip');
        KbQueueWait(0,3);



        % finalize the training
        DrawFormattedText(window,'Entraînement terminé','center','center',white);
        DrawFormattedText(window,pressWhenReadyText,pressWhenReadyX,pressWhenReadyY,white,wrapat);
        Screen(window,'Flip');
        KbQueueWait(0,3);
end

%% social task
switch taskToPerform.socialTask
    case 'on'
        % introduce the task
        DrawFormattedText(window,'Tâche avec visages','center','center',white);
        DrawFormattedText(window,pressWhenReadyText,pressWhenReadyX,pressWhenReadyY,white,wrapat);
        Screen(window,'Flip');
        KbQueueWait(0,3);


        % finalize the task
        DrawFormattedText(window,'Bravo! Vous avez terminé l''expérience!',...
            'center','center',white,wrapat);
        DrawFormattedText(window,pressWhenReadyText,pressWhenReadyX,pressWhenReadyY,white,wrapat);
        Screen(window,'Flip');
        KbQueueWait(0,3);
end

%% stop checking buttons: release key buffer
KbQueueStop;
KbQueueRelease;

%% close PTB
ShowCursor;
sca;