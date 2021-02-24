function[mentalE_perf, trial_success] = LGCM_mental_effort_perf(scr, stim, key,...
    numberVector, mentalE_prm, n_max_to_reach, instructions_disp, time_limit, t_max)
%[mentalE_perf, trial_success] = LGCM_mental_effort_perf(scr, stim, key,...
%     numberVector, mentalE_prm, n_max_to_reach, instructions_disp, time_limit, t_max)
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
% instructions_disp:
% (0) no reminder of what the question is nor of where you should answer
% (1) display instructions: ask if odd/even (or lower/higher than 5) and
% display also on the screen the relevant answer to each question
%
% OUTPUTS
% mentalE_perf: structure with summary of mental effort performance
%   .nTrials: number of trials it took to reach a correct
% performance
%   .rt: reaction time (rt) for each question
%   .taskType: task type (0: odd/even; 1: lower/higher than 5)
%
% See also LGCM_mental_learning.m

%% extract relevant variables
% angle values
startAngle = mentalE_prm.startAngle;
endAngle = 360;
totalAngleDistance = endAngle - startAngle;

% extract main mental effort parameters
sideQuestion = mentalE_prm.sideQuestion;
mental_n_col = mentalE_prm.mental_n_col;
switchPerc = mentalE_prm.switchPerc;

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
    startAngle, endAngle,...
    sideQuestion, taskType(1), numberVector(1), mental_n_col, instructions_disp);
onset_question_tmp = onsetTrial; % for the first question
timeNow = onsetTrial;

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
        % without this, it can create weird bugs as for the next question
        % it could reuse the previous answer
        KbReleaseWait();
        
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
% record number of questions answered and how many were correct
mentalE_perf.n_questions_performed = i_question - 1;
mentalE_perf.n_questions_correct = i_max_correct;
% record if trial was achieved or interrompted due to time limit (=failure)
if i_max_correct == n_max_to_reach
    trial_success = true;
elseif i_max_correct < n_max_to_reach
    trial_success = false;
end
mentalE_perf.success = trial_success;

end % function