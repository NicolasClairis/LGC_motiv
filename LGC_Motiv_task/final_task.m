function[] = final_task()
% final_task will launch the final task to be performed at the end of the
% experiment. Its purpose it to recompute a new indifference point at the
% end of the experiment to see whether participants are impacted by fatigue
% in their final subjective estimation of how much effort they would be
% willing to make.
%
% See also main_TRAINING and choice_task_main for the previous parts of the
% experiment.
%
% Designed by N.Clairis - june 2023

%% start where the function is located
scriptPath = [fileparts(matlab.desktop.editor.getActiveFilename),filesep];
cd(scriptPath);

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

%% subject number?
% physical/mental order?
iSubject = [];
while isempty(iSubject) || length(iSubject) ~= 3
    % repeat until all questions are answered
    info = inputdlg({'Subject CID (XXX)','p/m'});
    [iSubject,p_or_m] = info{[1,2]};
    %     warning('when real experiment starts, remember to block in french and IRM = 1');
end
% Create subjectCodeName which is used as a file saving name
subjectCodeName = strcat('CID',iSubject);
subResultFolder = [results_folder, subjectCodeName, filesep,...
    'behavior',filesep];

if ~ismember(p_or_m,{'p','m'})
    error('this letter has no definition');
end

%% final file name
finalTask_fileName = ['finalTask_data_',subjectCodeName,'.mat'];

%% load indifference point (IP)
file_nm_IP = ['delta_IP_CID',num2str(iSubject)];
full_IP_filename = [subResultFolder,file_nm_IP,'.mat'];
if exist(full_IP_filename,'file')
    IP_variables = getfield(load(full_IP_filename,'IP_variables'),'IP_variables');
else
    error(['Could not find ',full_IP_filename,' please check path and file name is ok.']);
end

%% timings
t_instru = 2;
timings.cross.mainTask = 0.5;
% precise if the choice and the performance periods will have a time
% constraint
choiceTimeParameters.timeLimit = false;
t_dispChosen    = 2; % keep same timing as in main task
% final time
t_endSession = 3;

%% initialize screen
IRMdisp = 0; % defines the screen parameters (0 for training screen, 1 for fMRI screen)
testing_script = 1;
[scr, ~, ~,...
    window] = ScreenConfiguration(IRMdisp, testing_script);
white = scr.colours.white;
% black = scr.colours.black;

%% general parameters of the task
IRMbuttons = 1; % defines the buttons to use (1 = same as in fMRI)
% initialize visual stimuli to use in the experiment
langage = 'fr';
n_E_levels = 4;
[stim] = stim_initialize(scr, n_E_levels, langage);

% number of buttons to answer
switch IRMbuttons
    case 0
        n_buttonsChoice = 4;
    case 1 % test buttons
        n_buttonsChoice = 4;
end
% define relevant keys and dynamometer module
[key] = relevant_key_definition('mental', IRMbuttons, n_buttonsChoice);
keyboard_check_start(key, 0);
% confidence feedback visual display
confidenceDispChosen.display = true;

%% perform the task
% effort
% baseline level of effort difficulty
baselineE = 2;
% number of repetitions for the fixed effort option
E_lowE_nRepeats = 1;

% Baseline Reward (CHF)
baselineR = 0.5;
R_or_P = 'R';

% number of trial repetitions until considering task is over if the
% participant doesn't reach the breakpoint before
nTrials = 10; % max should be at 1024

% general introduction
DrawFormattedText(window,['Dans cette derniere partie, nous allons vous poser quelques questions ',....
    'sur vos preferences en general pour chaque type d''effort. Cette fois, ',...
    'vous n''aurez plus besoin d''executer les efforts mais essayez de repondre ',...
    'honnetement.'],...
    'center','center',white, scr.wrapat);
% display text: Press when you are ready to start
DrawFormattedText(window, stim.pressWhenReady.text,...
    stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour);
[~, onsets.finalTask_general_instructions] = Screen(window, 'Flip');
KbQueueWait(0,3);

% define order of mental/physical tasks
switch p_or_m
    case 'p'
        taskOrder = {'p','m'};
    case 'm'
        taskOrder = {'m','p'};
end
for iTask = 1:length(taskOrder)
    task_nm = taskOrder{iTask};
    E_nm = ['E',task_nm];

    %% introduction of the task
    switch task_nm
        case 'm'
            DrawFormattedText(window,'pour les efforts MENTAUX',...
                'center','center',white);
        case 'p'
            DrawFormattedText(window,'pour les efforts PHYSIQUES',...
                'center','center',white);
    end
    % display text: Press when you are ready to start
    DrawFormattedText(window, stim.pressWhenReady.text,...
        stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour);
    [~, onsets.(E_nm).finalTask_instructions] = Screen(window, 'Flip');
    WaitSecs(t_instru);

    %% extract information regarding fixed reward and effort values
    % reward
    switch task_nm
        case 'm'
            deltaIP = IP_variables.mentalDeltaIP;
        case 'p'
            deltaIP = IP_variables.physicalDeltaIP;
    end
    high_R.(E_nm) = round(baselineR + 10*deltaIP,2);
    R_left = baselineR;
    R_right = high_R.(E_nm);
    % effort
    E_left = baselineE;
    E_right = baselineE;
    
    % time variables
    [onsets.(E_nm).preChoiceCross,...
        onsets.(E_nm).dispChoiceOptions,...
        onsets.(E_nm).choice,...
        onsets.(E_nm).preChoiceCross_keyReleaseMessage,...
        onsets.(E_nm).preChoiceCross_after_buttonRelease,...
        onsets.(E_nm).dispChosen,...
        dur.(E_nm).preChoiceCross,...
        dur.(E_nm).dispChoiceOptions,...
        dur.(E_nm).preChoiceCross_keyReleaseMessage,...
        dur.(E_nm).preChoiceCross_after_buttonRelease,...
        dur.(E_nm).dispChosen,...
        was_a_key_pressed_bf_trial.(E_nm)] = deal(NaN(1,nTrials));
    % main variables
    [RT.(E_nm),...
        confidence.(E_nm),...
        R_chosen.(E_nm),...
        E_chosen.(E_nm),...
        E_chosen_repeats.(E_nm),...
        high_E_nRepeats.(E_nm)] = deal(NaN(1,nTrials));
    choice_LR.(E_nm) = zeros(1,nTrials);
    breakPointReached.(E_nm) = 0;
    iTrial = 0;
    while (iTrial < nTrials) && breakPointReached.(E_nm) == 0
        % increase trial number
        iTrial = iTrial + 1;
        
        %% fixation cross
        Screen('FillRect',window, white, stim.cross.verticalLine); % vertical line
        Screen('FillRect',window, white, stim.cross.horizontalLine); % horizontal line
        [~,onsets.(E_nm).preChoiceCross(iTrial)] = Screen('Flip',window); % display the cross on screen
        WaitSecs(timings.cross.mainTask);
        dur.(E_nm).preChoiceCross(iTrial) = GetSecs - onsets.(E_nm).preChoiceCross(iTrial);
        
        %% check that no key is being pressed before the choice trial starts
        [was_a_key_pressed_bf_trial.(E_nm)(iTrial),...
            onsets.(E_nm).keyReleaseMessage(iTrial),...
            dur.(E_nm).preChoiceCross_keyReleaseMessage(iTrial)] = check_keys_are_up(scr, stim, key);
        
        % if a key was pressed before starting the trial => show the fixation
        % cross again with a similar amount of time
        if was_a_key_pressed_bf_trial.(E_nm)(iTrial) == 1
            Screen('FillRect',window,white, stim.cross.verticalLine); % vertical line
            Screen('FillRect',window,white, stim.cross.horizontalLine); % horizontal line
            [~,onsets.(E_nm).cross_after_buttonRelease(iTrial)] = Screen('Flip',window); % display the cross on screen
            WaitSecs(1);
            dur.(E_nm).preChoiceCross_after_buttonRelease(iTrial) = GetSecs - onsets.(E_nm).preChoiceCross_after_buttonRelease(iTrial);
        end
        
        %% choice period
        % define effort difficulty
        E_left_nRepeats = E_lowE_nRepeats;
        if iTrial == 1
            high_E_nRepeats.(E_nm)(1) = 2;
        else
            high_E_nRepeats.(E_nm)(iTrial) = high_E_nRepeats.(E_nm)(iTrial-1)*2;
        end
        E_right_nRepeats = high_E_nRepeats.(E_nm)(iTrial);
        
        % keep choice period until a choice is done
        while choice_LR.(E_nm)(iTrial) == 0
            [choice_LR.(E_nm)(iTrial),...
                onsets.(E_nm).dispChoiceOptions(iTrial),...
                onsets.(E_nm).choice(iTrial)] = final_task_choice_period(scr, stim,...
                R_left, R_right, E_left, E_right,...
                E_left_nRepeats, E_right_nRepeats,...
                R_or_P,...
                choiceTimeParameters, key, task_nm);
        end % keep performing the trial until a choice is made

        % store information relative to choice made
        R_chosen.(E_nm)(iTrial) = R_left.*(choice_LR.(E_nm)(iTrial)<0) +...
            R_right.*(choice_LR.(E_nm)(iTrial)>0);
        E_chosen.(E_nm)(iTrial) = E_left.*(choice_LR.(E_nm)(iTrial)<0) +...
            E_right.*(choice_LR.(E_nm)(iTrial)>0);
        E_chosen_repeats.(E_nm)(iTrial) = E_left_nRepeats.*(choice_LR.(E_nm)(iTrial)<0) +...
            E_right_nRepeats.*(choice_LR.(E_nm)(iTrial)>0);
        confidence.(E_nm)(iTrial) = abs(choice_LR.(E_nm)(iTrial)) == 2;
        RT.(E_nm) = onsets.(E_nm).choice(iTrial) - onsets.(E_nm).dispChoiceOptions(iTrial);

        %% display chosen option
        confidenceDispChosen.lowOrHigh = abs(choice_LR.(E_nm)(iTrial));
        [time_dispChosen] = final_task_dispChosen(scr, stim, choice_LR.(E_nm)(iTrial),...
            R_chosen.(E_nm)(iTrial), E_chosen.(E_nm)(iTrial), E_chosen_repeats.(E_nm)(iTrial),...
            R_or_P, confidenceDispChosen);
        onsets.(E_nm).dispChosen(iTrial) = time_dispChosen;
        WaitSecs(t_dispChosen);
        dur.(E_nm).dispChosen(iTrial) = GetSecs - onsets.(E_nm).dispChosen(iTrial);
        
        %% check if break point has been reached
        if R_chosen.(E_nm)(iTrial) == baselineR
            breakPointReached.(E_nm) = 1;
            breakPointValue.(E_nm) = high_E_nRepeats.(E_nm)(iTrial);
        end
    end % trial loop
end % task loop

%% save the data
save([subResultFolder, finalTask_fileName],...
    'choice_LR',...
    'breakPointReached',...
    'confidence',...
    'onsets','dur','RT',...
    'R_chosen','E_chosen','E_chosen_repeats',...
    'high_E_nRepeats','breakPointValue');
    
%% end message
DrawFormattedText(window,...
    'Felicitations! L''experience est maintenant terminee.',...
    'center', 'center', white, scr.wrapat);
Screen(window,'Flip');
WaitSecs(t_endSession);

%% releyse buffer for key presses
KbQueueStop;
KbQueueRelease;

%% close PTB
ShowCursor;
sca;

end % function