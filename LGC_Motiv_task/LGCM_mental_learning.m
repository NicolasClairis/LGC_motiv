function[learning_summary] = LGCM_mental_learning(scr, stim, key,...
    n_max_to_reach,...
    mental_n_col, sideQuestion, switchPerc, instructions_disp)
% [learning_summary] = LGCM_mental_learning(scr, stim, key,...
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
% key: structure with key code for Psychtoolbox to identify which key
% corresponds to left and right cues
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
xScreenCenter = scr.xCenter;
yScreenCenter = scr.yCenter;
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
    rt,...
    sideAnswer,...
    goodOrBadAnswer] = deal( NaN(1, n_questions));
% define first trial task type
taskType(1)    = random_binary;


%% display color associated to each question type: pair/impair vs low/high than 5
oddEven_col = mental_n_col.oddEven;
lowHigh_col = mental_n_col.lowHigh;

blackCol = [0 0 0];
wrapat = 50; % go to line for DrawFormattedText when more than this number of characters

%% display general instruction on the screen
% adapt text depending on if instructions included or not

for iTimeLoop = 1:2
    switch instructions_disp
        case 1 % WITH INSTRUCTIONS
            DrawFormattedText(window,...
                ['Vous allez pouvoir vous entra�ner � effectuer la t�che. ',...
                'La couleur du chiffre repr�sente la question pos�e.',...
                'Vous pouvez r�pondre avec les boutons. ',...
                'La correspondance entre la position des boutons et la ',...
                'r�ponse que vous souhaitez donner est la suivante:'],...
                'center', yScreenCenter/3, blackCol, wrapat);
        case 0
            DrawFormattedText(window,...
                ['Vous allez pouvoir vous entra�ner � effectuer la t�che mais vous ',...
                ' devrez vous rappeler de la correspondance entre les couleurs et la ',...
                ' question � laquelle il faut r�pondre. /n Pour rappel:'],...
                'center', yScreenCenter/3, blackCol, wrapat);
    end
    
    % odd/even info
    DrawFormattedText(window,'Pair ou Impair?',...
        'center',yScreenCenter, oddEven_col);
    DrawFormattedText(window, 'Pair',...
        xScreenCenter/2, yScreenCenter*5/4, oddEven_col);
    DrawFormattedText(window, 'Impair',...
        xScreenCenter*3/2, yScreenCenter*5/4, oddEven_col);
    
    % lower/higher than 5 info
    DrawFormattedText(window,'< ou > 5?',...
        'center',yScreenCenter*6/4, lowHigh_col);
    DrawFormattedText(window, '< 5',...
        xScreenCenter/2, yScreenCenter*7/4, lowHigh_col);
    DrawFormattedText(window, '> 5',...
        xScreenCenter*3/2, yScreenCenter*7/4, lowHigh_col);
    
    if iTimeLoop == 1 % force them to read at first
        Screen(window,'Flip');
        WaitSecs(t_instructions);
    elseif iTimeLoop == 2 % after t_instructions seconds, they can manually start
        DrawFormattedText(window, 'Appuyez quand vous �tes pr�ts � commencer la t�che.',...
            'center', yScreenCenter*15/8, blackCol);
        Screen(window,'Flip');
        KbWait;
    end
end % loop over forced reading/manual pass loop

%% perform the learning
i_max_correct = 0;
i_question = 1;

% initial display to get timing of the start
[onsetTrial] = LGCM_mental_display_stim(scr, stim,...
    startAngle, endAngle,...
    sideQuestion, taskType(1), numberVector(1), mental_n_col, instructions_disp);
onset_question_tmp = onsetTrial; % for the first question

% loop until relevant number of subsequent correct answers has been reached
while i_max_correct < n_max_to_reach
    % trial informations
    numberValue_tmp = numberVector(i_question);
    taskType_tmp = taskType(i_question);
    
    LGCM_mental_display_stim(scr, stim,...
        startAngle, endAngle,...
        sideQuestion, taskType_tmp, numberValue_tmp, mental_n_col, instructions_disp);
    
    %% check key presses
    [keyisdown, timeAnswer, keyCode] = KbCheck();
    
    if (keyisdown == 1) &&...
            ((keyCode(key.left) == 1 && keyCode(key.right) == 0) ||...
            (keyCode(key.left) == 0 && keyCode(key.right) == 1)) % focus only when 1 single button
        % which belongs to the 2 buttons of interest has been pressed
        
        if keyCode(key.left) == 1 && keyCode(key.right) == 0 % left answer
            sideAnswer_tmp = -1;
        elseif keyCode(key.left) == 0 && keyCode(key.right) == 1 % right answer
            sideAnswer_tmp = 1;
        end % left or right answer? (ignore the rest = if another key has
        % been pressed or if both pressed in the same time for ex.)
        
        %% wait for participant to release the keys before updating the variables of interest
        KbReleaseWait(); % without this,
        % it can create weird bugs as for the next question it could reuse the previous answer
        
        %% determine whether the response was correct or not
        [answerCorrect_tmp] = LGCM_mental_effort_answer_correct(taskType_tmp,...
            numberValue_tmp,...
            sideAnswer_tmp, sideQuestion);
        
        switch answerCorrect_tmp
            case 0 % error made
                % if wrong, set back indicators to zero: needs to restart
                startAngle = 0;
                i_max_correct = 0;
            case 1 % if correct, update the count and the display
                startAngle = startAngle + totalAngleDistance/n_max_to_reach;
                i_max_correct = i_max_correct + 1;
        end
        
        %% define the task type for the next question
        if i_max_correct >= 1 && i_max_correct < n_max_to_reach
            % once a correct answer has been provided you need to assign the
            % moments when a switch will happen
            % = define a STAY/SWITCH sequence
            % (as long as errors are being made, keep the task easy)
            %
            % once the expected amount of correct answers has been reached,
            % there is no need to make more switches
            
            if i_max_correct == 1
                [task_seq_tmp] = LGCM_mental_effort_task_switches(taskType_tmp, n_max_to_reach, n_switch);
            end % first correct answer of a sequence
            
            % extract
            taskType(i_question + 1) = task_seq_tmp(i_max_correct + 1);
            
        elseif i_max_correct == 0 % keep doing the same task as long as the participant makes an error
            % so that the task is clearly understood before adding the
            % switches
            taskType(i_question + 1) = taskType(i_question);
        end % at least one correct answer filter
        
        %% update variables of interest
        % record time to answer and re-initialize the counter to get next
        % RT
        rt(i_question) = timeAnswer - onset_question_tmp;
        onset_question_tmp = timeAnswer; % for the next question, the onset corresponds to the answer of the previous question
        % record side and nature of the answer provided (correct/wrong?)
        sideAnswer(i_question) = sideAnswer_tmp;
        goodOrBadAnswer(i_question) = answerCorrect_tmp;
        
        %% update total count of answers in any case
        i_question = i_question + 1;
        
    end % has a relevant key been pressed
    
end % keep performing learning until number of subsequent answers reaches threshold predefined

%% record all in output
% keep only questions performed
questions_done = ~isnan(sideAnswer);
% record question parameters
learning_summary.numberVector   = numberVector(questions_done);
learning_summary.taskType       = taskType(questions_done);
learning_summary.sideAnswer     = sideAnswer(questions_done);
learning_summary.isGoodAnswer   = goodOrBadAnswer(questions_done);
learning_summary.rt             = rt(questions_done);
learning_summary.n_switch       = n_switch;
learning_summary.n_max_to_reach = n_max_to_reach;
% record output
learning_summary.n_trials_required_forLearning = i_question - 1;


end % function