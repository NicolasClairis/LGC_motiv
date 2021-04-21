function[mentalE_perf, trial_success, onsets] = LGCM_mental_effort_perf(scr, stim, key,...
    numberVector, mentalE_prm, n_max_to_reach,...
    learning_col, learning_instructions, time_limit, t_max)
%[mentalE_perf, trial_success, onsets] = LGCM_mental_effort_perf(scr, stim, key,...
%     numberVector, mentalE_prm, n_max_to_reach,...
%     curr_learning_col, curr_learning_instructions, time_limit, t_max)
%
% LGCM_mental_effort_perf corresponds to the actual performance. Can be
% used both for learning period (with or without instructions) and for the
% actual task. It corresponds to one mental effort trial (either learning,
% calibration or task trial).
%
% INPUTS
% scr: structure with screen parameters
%
% stim: stucture with stimuli parameters
%
% key: structure with key code for Psychtoolbox to identify which key
% corresponds to left and right cues
%
% numberVector: vector with big list of numbers for the current trial
%
% mentalE_prm: structure with main parameters for mental effort task:
%   .startAngle: start angle for displaying the effort arc
%
%   .mental_n_col: structure with the colour to use for the font
%       .oddEven: colour to use for odd or even question
%       .lowHigh: colour to use for lower/higher than 5 question
%
%   .sideQuestion: structure which tells where each answer is expected to be
%       .oE.pair: -1 means left button corresponds to the pair answer and
%       .oE.impair = 1 means right button corresponds to the impair answer
%       same logic applies to .hL.low and .hL.high fields
%
%   .switchPerc: percentage of switches required (based on total number of
%   subsequent correct answers you want)
%
% learning_col:
% 'col1': learning with colour 1 only
% 'col2': learning with colour 2 only
% 'all': learning with both colours
%
% learning_instructions
% 'fullInstructions': display instructions: ask if odd/even (or lower/higher than 5) and
% display also on the screen the relevant answer to each question
% 'partialInstructions': display only the two possible answers but not the
% question anymore
% 'noInstructions': no reminder of what the question is nor of where you should answer
%
% time_limit
% (0)/false = no time limit to perform the requested number of questions
% (1)/true = when time limit is reached, stop the trial, even if requested
% number not reached => consider the trial a failure
%
% t_max: maximum time allowed to reach the requested number of correct
% answers
%
% OUTPUTS
% mentalE_perf: structure with summary of mental effort performance
%   .nTrials: number of trials it took to reach a correct
% performance
%   .rt: reaction time (rt) for each question
%   .taskType: task type (0: odd/even; 1: lower/higher than 5)
%   .totalTime_success: time passed between trial onset and last correct
%   answer (if the trial was a success), NaN if trial was a failure
%
% trial_success:
% false: failure
% true: effort was correctly performed in the requested time
%
% onsets: structure with onset of each number on the screen
%   .nb_i: onset of (i) question during the current trial
%
% See also LGCM_mental_learning.m

%% extract relevant variables
% angle values
startAngle_currentTrial = mentalE_prm.startAngle;
endAngle = 360;
totalAngleDistance = endAngle - startAngle_currentTrial;

% extract main mental effort parameters
sideQuestion = mentalE_prm.sideQuestion;
mental_n_col = mentalE_prm.mental_n_col;
switchPerc = mentalE_prm.switchPerc;

% determine number of switches to implement in a given sequence
% with instructions
if strcmp(learning_col,'all')
    n_switch = switchPerc*n_max_to_reach;
    if n_switch ~= round(n_switch)
        error('something is wrong with your script, number of switch has to be integer');
    end
else % no switch if learning session focusing on one single colour
    n_switch = 0;
end

%% define moments when the task switches: should be unpredictable and
% distance between 2 switches should vary between 1 and 4
n_questions = size(numberVector,2);
[taskType,...
    rt,...
    sideAnswer,...
    goodOrBadAnswer] = deal( NaN(1, n_questions));

% define first trial task type
switch learning_col
    case 'all'
        taskType(1) = random_binary; % define randomly the nature of the first trial
    case 'col1' % start with colour 1
        switch mentalE_prm.mental_n_col.col1
            case 'oddEven'
                taskType(1) = 0;
            case 'lowHigh'
                taskType(1) = 1;
        end
    case 'col2' % start with colour 2
        switch mentalE_prm.mental_n_col.col2
            case 'oddEven'
                taskType(1) = 0;
            case 'lowHigh'
                taskType(1) = 1;
        end
end

%% initialize the counters
% number of subsequent correct answers
i_max_correct = 0;
% number of questions answered
i_question = 1;

%% wait all keys are released before starting
KbReleaseWait;
% you could even add a timer here in order to ask the participant to
% release the keys if he/she is currently pressing one key

%% initial display to get timing of the start
[onsetTrial] = LGCM_mental_display_stim(scr, stim,...
    startAngle_currentTrial, endAngle,...
    sideQuestion, taskType(1), numberVector(1), mental_n_col,...
    learning_instructions);
onset_question_tmp = onsetTrial; % for the first question
onsets.nb_1 = onsetTrial;
timeNow = onsetTrial;

startAngle = startAngle_currentTrial;

% keep track of number of errors made during the trial
jErrorsMade = 0;

% loop until relevant number of subsequent correct answers has been reached
% or that max time limit has been reached (if one time limit has been
% defined)
while (i_max_correct < n_max_to_reach) &&...
        ( ( (time_limit == true) && (timeNow < onsetTrial + t_max) ) ||...
        (time_limit == false) )
    timeNow = GetSecs;
    
    % trial informations
    numberValue_tmp = numberVector(i_question);
    taskType_tmp = taskType(i_question);
    
    % display instructions after 2 errors have been made (in case where no
    % instructions on screen)
    if strcmp(learning_instructions,'noInstructions') &&...
            i_question > 1 &&...
            mod(jErrorsMade, 2) == 0 &&...
            goodOrBadAnswer(i_question  - 1) == 0
        % will display instructions for the next question everytime after
        % 2 errors (can be adapted easily to remain for th whole trial if
        % necessary)
        learning_instructions_bis = 'fullInstructions';
    else
        learning_instructions_bis = learning_instructions;
    end
    
    onset_stim = LGCM_mental_display_stim(scr, stim,...
        startAngle, endAngle,...
        sideQuestion, taskType_tmp, numberValue_tmp, mental_n_col,...
        learning_instructions_bis);
    
    %% record onset
    if i_question > 1
        onsets.(['nb_',num2str(i_question)]) = onset_stim;
    end
    
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
        % without this, it can create weird bugs as for the next question
        % it could reuse the previous answer
        KbReleaseWait();
        
        %% determine whether the response was correct or not
        [answerCorrect_tmp] = LGCM_mental_effort_answer_correct(taskType_tmp,...
            numberValue_tmp,...
            sideAnswer_tmp, sideQuestion);
        
        switch answerCorrect_tmp
            case 0 % error made
                % version where you reset the timer whenever they make an
                % error
                % startAngle = startAngle_currentTrial; % re-initialize
                % i_max_correct = 0; % if wrong, set back indicators to zero: needs to restart
                
                % just (-1) decrement after an error (otherwise too hard)
                startAngle = startAngle - totalAngleDistance/n_max_to_reach;
                i_max_correct = i_max_correct - 1; % if wrong, decrement the total number of correct answers
                jErrorsMade = jErrorsMade + 1;
            case 1 % if correct, update the count and the display
                startAngle = startAngle + totalAngleDistance/n_max_to_reach;
                i_max_correct = i_max_correct + 1;
        end
        
        %% define the task type for the next question
        if i_max_correct == 1
            % once a correct answer has been provided you need to assign the
            % moments when a switch will happen
            % = define a STAY/SWITCH sequence
            [task_seq_tmp] = LGCM_mental_effort_task_switches(taskType_tmp, n_max_to_reach, n_switch);
        end % first correct answer of a sequence
        
        % define the task type for the next question
        if i_max_correct < n_max_to_reach
            % once the expected amount of correct answers has been reached,
            % there is no need to make more switches
            
            if answerCorrect_tmp == 0 % no task switch after an error to keep the task easy after an error has been made
                taskType(i_question + 1) = taskType(i_question);
            else
                taskType(i_question + 1) = task_seq_tmp(i_max_correct + 1);
            end
        end % no need to update task type for the next question if the end has been reached
        
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
    
end % keep performing until number of subsequent answers reaches threshold predefined or if timer has been reached

%% record all in output
% keep only questions performed
questions_done = ~isnan(sideAnswer);
% record question parameters
mentalE_perf.numberVector   = numberVector(questions_done);
mentalE_perf.taskType       = taskType(questions_done);
mentalE_perf.sideAnswer     = sideAnswer(questions_done);
mentalE_perf.isGoodAnswer   = goodOrBadAnswer(questions_done);
mentalE_perf.rt             = rt(questions_done);
mentalE_perf.n_switch       = n_switch;
mentalE_perf.n_max_to_reach = n_max_to_reach;
mentalE_perf.n_errorsMade   = jErrorsMade;
% record number of questions answered and how many were correct
mentalE_perf.n_questions_performed = i_question - 1;
mentalE_perf.n_questions_correct = i_max_correct;
% record if trial was achieved or interrompted due to time limit (=failure)
if i_max_correct == n_max_to_reach % reached the top
    trial_success = true;
    totalTime_success = timeAnswer - onsetTrial; % total time between start of the trial and last correct answer
elseif i_max_correct < n_max_to_reach % not enough good answers
    trial_success = false;
    totalTime_success = NaN;
end
mentalE_perf.success = trial_success;
mentalE_perf.totalTime_success = totalTime_success;

end % function