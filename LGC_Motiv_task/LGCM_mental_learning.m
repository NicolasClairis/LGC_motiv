function[learning_summary] = LGCM_mental_learning(scr, stim,...
    n_max_to_reach,...
    mental_n_col, sideQuestion, switchPerc, instructions_disp)
% [learning_summary] = LGCM_mental_learning(scr, stim,...
%     n_max_to_reach,...
%     mental_n_col, sideQuestion, switchPerc, instructions_disp)
% LGCM_mental_learning will train the participant to know how to perform
% the mental effort task with task switching correctly.
%
% INPUTS
% scr: structure with screen parameters
%
% stim: stucture with stimuli parameters
%
% n_max_to_reach: number of subsequent correct answers to
% provide until we consider the learning is ok
%
% mental_n_col: structure with the colour to use for the font
%   .oddEven: colour to use for odd or even question
%   .lowHigh: colour to use for lower/higher than 5 question
%
% sideQuestion: structure which tells where each answer is expected to be
% .oE.pair: -1 means left button corresponds to the pair answer and
% .oE.impair = 1 means right button corresponds to the impair answer
% same logic applies to .hL.low and .hL.high fields
%
% switchPerc: percentage of switches required (based on total number of
% subsequent correct answers you want)
%
% instructions_disp:
% (0) no reminder of what the question is nor of where you should answer
% (1) display instructions: ask if odd/even (or lower/higher than 5) and
% display also on the screen the relevant answer to each question
%
% OUTPUTS
% learning_summary: structure with all relevant information of the learning
% phase:
%   .nTrials: number of trials it took to reach a correct
% performance
%   .rt: reaction time (rt) for each question
%   .taskType: task type (0: odd/even; 1: lower/higher than 5)


%% extract and define main variables of interest
window = scr.window;
% angle values
startAngle = 0;
endAngle = 360;
totalAngleDistance = endAngle - startAngle;

t_instructions = 5; % minimal amount to see the instructions

% extract numbers to use for each learning phase
[numberVector] = LGCM_mental_numbers(1);

% determine number of switches to implement in a given sequence
% with instructions
n_switch = switchPerc*n_max_to_reach;
if n_switch ~= round(n_switch)
    error('something is wrong with your script, number of switch has to be integer');
end

%% define moments when the task switches: should be unpredictable and
% distance between 2 switches should vary between 1 and 4
n_questions = size(numberVector,2);
[taskType,...
    sideAnswer] = deal( NaN(n_questions));
% define first trial task type
taskType(1)    = random_binary;


%% display color associated to each question type: pair/impair vs low/high than 5
oddEven_col = mental_n_col.oddEven;
lowHigh_col = mental_n_col.lowHigh;

%% display general instruction on the screen
% adapt text depending on if instructions included or not
switch instructions_disp
    case 0
        DrawFormattedText(window, xScreenCenter, yScreenCenter,...
            ['Vous allez pouvoir vous entraîner à effectuer la tâche mais vous ',...
            ' devrez vous rappeler de la correspondance entre les couleurs et la ',...
            ' question à laquelle il faut répondre.']);
        Screen(window,'Flip');
        WaitSecs(t_instructions);
    case 1
        DrawFormattedText(window, xScreenCenter, yScreenCenter,...
            ['Vous allez pouvoir vous entraîner à effectuer la tâche.']);
        DrawFormattedText(window, xScreenCenter, yScreenCenter*5/6,...
            'Appuyez quand vous êtes prêts à commencer la tâche.');
        Screen(window,'Flip');
        KbWait;
end

%% perform the learning
i_max_correct = 0;
i_question = 1;
while i_max_correct < n_max_to_reach
    % trial informations
    numberValue_tmp = numberVector(i_question);
    taskType_tmp = taskType(i_question);
    
    % display percentage of correct answers already
    % provided
    Screen('FillArc', window,...
        stim.difficulty.currLevelColor,...
        stim.difficulty.middle_center,...
        trial_startAngle,...
        endAngle - trial_startAngle);
    
    % display number to solve
    switch taskType_tmp
        case 0 % odd/even
            textColor = oddEven_col;
        case 1 % lower/higher than 5
            textColor = lowHigh_col;
    end
    DrawFormattedText(window, num2str(numberValue_tmp), xScreenCenter, yScreenCenter*(1/6), textColor);
    
    % display instructions
    if instructions_disp
        switch taskType_tmp
            case 0
                DrawFormattedText(window, 'Chiffre pair ou impair?',xScreenCenter, yScreenCenter/3, textColor);
            case 1
                DrawFormattedText(window, 'Chiffre < ou > 5?',xScreenCenter, yScreenCenter/3, textColor);
        end
        LGCM_mental_effort_task_question_display(scr, task_trialType, sideQuestion, textColor);
    end
    
    [~,onsetTrial] = Screen(window,'Flip');
    
    %% check key presses
    [keyisdown, timeAnswer, keyCode] = KbCheck();
    
    if (keyisdown == 1) &&...
            ((keyCode(key.left) == 1 && keyCode(key.right) == 0) ||...
            (keyCode(key.left) == 0 && keyCode(key.right) == 1))
        
        if keyCode(key.left) == 1 && keyCode(key.right) == 0 % left answer
            sideAnswer_tmp = -1;
        elseif keyCode(key.left) == 0 && keyCode(key.right) == 1 % right answer
            sideAnswer_tmp = 1;
        end % left or right answer? (ignore the rest = if another key has
        % been pressed or if both pressed in the same time for ex.)
        
        %% determine whether the response was correct or not
        [answerCorrect_tmp] = LGCM_mental_effort_answer_correct(taskType_tmp,...
            numberValue_tmp,...
            sideAnswer_tmp, sideQuestion);
        
        switch answerCorrect_tmp
            case 0 % error made
                % if wrong, set back to zero: needs to restart
                trial_startAngle = 0;
                i_max_correct = 0;
            case 1 % if correct, update the count and the display
                trial_startAngle = trial_startAngle + totalAngleDistance/n_max_to_reach;
                i_max_correct = i_max_correct + 1;
        end
        
        % record current task type
        taskType(i_question) = taskType_tmp;
        
        %% define the task type for the next question
        if i_max_correct >= 1
            % once a correct answer has been provided you need to assign the
            % moments when a switch will happen
            % = define a STAY/SWITCH sequence
            % (as long as errors are being made, keep the task easy)
            
            if i_max_correct == 1
                [task_seq_tmp] = LGCM_mental_effort_task_switches(taskType_tmp, n_max_to_reach, n_switch);
            end % first correct answer of a sequence
            
            % extract 
            taskType(i_question) = taskType_tmp;
            taskType(i_question + 1) = task_seq_tmp(i_max_correct + 1);
        end % at least one correct answer filter
        
        %% update variables of interest
        % record time to answer and side of the answer
        rt.trialsctions(i_question) = timeAnswer - onsetTrial;
        sideAnswer(i_question) = sideAnswer_tmp;
        
        %% update total count of answers in any case
        i_question = i_question + 1;
        
    end % has a relevant key been pressed
    
end % keep performing learning until number of subsequent answers reaches threshold predefined


end % function