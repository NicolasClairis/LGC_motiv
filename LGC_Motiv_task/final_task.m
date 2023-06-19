function[] = final_task()


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
IRMdisp = 0; % defines the screen parameters (0 for training screen, 1 for fMRI screen)
IRMbuttons = 1; % defines the buttons to use (1 = same as in fMRI)
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
fileName = ['finalTask_data_',subjectCodeName,'.mat'];

%% load indifference point (IP)
file_nm_IP = ['delta_IP_CID',num2str(iSubject)];
IP_variables = getfield(load([subResultFolder,...
    file_nm_IP,'.mat'],'IP_variables'),'IP_variables');

%% initialize screen
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(IRMdisp, testing_script);
window = scr.window;
white = scr.colours.white;
black = scr.colours.black;
% confidence feedback visual display
confidenceDispChosen.display = true;
% no display of the confidence mapping
confidenceChoiceDisp = false;
%% timings
timings.cross.mainTask = 0.5;
% precise if the choice and the performance periods will have a time
% constraint
choiceTimeParameters.timeLimit = false;
t_dispChoice    = timings.dispChoice;


% initialize visual stimuli to use in the experiment
langage = 'fr';
[stim] = stim_initialize(scr, n_E_levels, langage);

% number of buttons to answer
switch IRMbuttons
    case 0
        n_buttonsChoice = 4;
    case 1 % test buttons
        n_buttonsChoice = 4;
end

%% perform the task
% for how many levels of effort will you measure the IP?
E_right = 2;
E_left  = 2;
E_left_nRepeats = 1;
E_right_nRepeats = 1;

% Baseline Reward (CHF)
baselineR = 0.5;
baselineE = 2;

% number of repetitions
nTrials = 10; % max should be at 1024

% final time
t_endSession = 3;

switch p_or_m
    case 'p'
        taskOrder = {'p','m'};
    case 'm'
        taskOrder = {'m','p'};
end
for iTask = 1:length(taskOrder)
    task_nm = taskOrder{iTask};
    switch task_nm
        case 'm'
            high_R.(E_nm) = baselineR + 10*IP_variables.mentalDeltaIP;
        case 'p'
            high_R.(E_nm) = baselineR + 10*IP_variables.physicalDeltaIP;
    end
    E_nm = ['E',task_nm];
    
    % time variables
    [onsets.(E_nm).preChoiceCross,...
        onsets.(E_nm).dispChoiceOptions,...
        onsets.(E_nm).choice,...
        onsets.(E_nm).preChoiceCross_keyReleaseMessage,...
        onsets.(E_nm).preChoiceCross_after_buttonRelease,...
        dur.(E_nm).preChoiceCross,...
        dur.(E_nm).dispChoiceOptions,...
        dur.(E_nm).preChoiceCross_keyReleaseMessage,...
        dur.(E_nm).preChoiceCross_after_buttonRelease,...
        was_a_key_pressed_bf_trial.(E_nm)] = deal(NaN(1,nTrials));
    % main variables
    [RT.(E_nm),...
        confidence.(E_nm),...
        R_chosen.(E_nm),...
        E_chosen.(E_nm),...
        high_E_nRepeats.(E_nm)] = deal(NaN(1,nTrials));
    choice_LR.(E_nm) = zeros(1,nTrials);
    breakPointReached.(E_nm) = 0;
    while (iTrial < nTrials) && breakPointReached.(E_nm) == 0
        
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
        if iTrial == 1
            high_E_nRepeats.(E_nm)(1) = 2;
        else
            high_E_nRepeats.(E_nm)(iTrial) = high_E_nRepeats.(E_nm)(iTrial-1)*2;
        end
        
        % keep choice period until a choice is done
        while choice_LR.(E_nm)(iTrial) == 0
            [choice_LR.(E_nm)(iTrial),...
                onsets.(E_nm).dispChoiceOptions(iTrial),...
                onsets.(E_nm).choice(iTrial)] = finalTask_choice_period(scr, stim,...
                R_left, R_right, E_left, E_right, R_or_P,...
                timeParameter, key, confidenceDisp);
        end % keep performing the trial until a choice is made
        
        %% display chosen option
        [time_dispChoice] = final_task_dispChosen(scr, stim, choice(iTrial),...
            R_chosen(iTrial), E_chosen(iTrial), R_or_P, confidenceDispChosen);
        onsets.(E_nm).dispChoice(iTrial) = time_dispChoice;
        WaitSecs(t_dispChoice);
        dur.(E_nm).dispChoice(iTrial) = GetSecs - onsets.(E_nm).dispChoice(iTrial);
        
        %% check if break point has been reached
        if choice_LR.(E_nm)(iTrial) == 1
            breakPointReached.(E_nm) = 1;
        end
    end % trial loop
end % task loop

%% save the data
save([subResultFolder, file_nm_IP,'.mat'],...
    'choice_LR',...
    'breakPointReached',...
    'confidence',...
    'onsets','dur','RT');
    
%% end message
DrawFormattedText(window,...
    'Felicitations! L''experience est maintenant terminee.',...
    'center', 'center', scr.colours.white, scr.wrapat);
Screen(window,'Flip');
WaitSecs(t_endSession);
waitSpace;

%% releyse buffer for key presses
KbQueueStop;
KbQueueRelease;

%% close PTB
ShowCursor;
sca;

end % function