function[trial_failed, perf, onsetEffortPeriod, onsetAnswers, onsetNumbers] = LGCM_mental_effort(scr, stim,...
    t_effort_max, n_correctAnswers,...
    R_chosen, R_or_P, E_chosen, n_E_levels,...
    numberVector, switchPerc, task_trialStart,...
    sideQuestion,...
    key)
% 
% LGCM_mental_effort sub-function of LGCM_choice_task_main.m to perform the
% mental effort task
%
% INPUTS
% scr: structure with screen parameters
%
% stim: stucture with stimuli parameters
% 
% t_effort_max: maximal time allowed to perform the effort
%
% n_correctAnswers: number of correct subsequent answers required for each
% difficulty level
%
% R_chosen: reward level chosen for this trial
%
% R_or_P: reward or punishment trial?
%
% E_chosen: effort level chosen for this trial => tells you how many
% correct answers have to be given subsequently to fill the trial
%
% n_E_levels: maximum number of effort levels
%
% numberVector: big vector with list of numbers
%
% switchPerc: percentage of switch for this task
%
% task_trialStart: string which indicates which is the first task to do
% (0) start with odd/even task
% (1) higher/lower than 5 task
%
% sideQuestion: structure with side for each answer (especially if you decide to vary it)
% sideQuestion.oE.pair, sideQuestion.oE.impair, sideQuestion.hL.low,
% sideQuestion.hL.high
% (-1) left
% (+1) right
%
% key: structure with relevant key codes for left and right key press
% 
% OUTPUTS
% trial_failed:
% (0): trial finished successfully
% (1): trial failed (participant did not solve enough questions in the time
% available)
%
% perf: vector with all the relevant information for each question
%
% onsetEffortPeriod: onset when participant can start performing the effort
%
% onsetAnswers: timing of every answer provided
%
% onsetNumbers: timing of every time a new number appears on screen

%% extract main variables of interest
window = scr.window;
% angle values
startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen)]);
endAngle = stim.difficulty.arcEndAngle;
totalAngleDistance = endAngle - startAngle;

% extract correspondancy between level of effort and number of subsequent
% correct answers to give
n_correctAnswersToGive = n_correctAnswers.(['effort_',num2str(E_chosen)]);
n_switch = switchPerc*n_correctAnswersToGive;
if n_switch ~= round(n_switch)
    error('something is wrong with your script, number of switch has to be integer');
end

[onsetAnswers,...
    onsetNumbers] = deal(NaN(1,30)); % initialize onset numbers

%% display all relevant information on screen
jNumber = 1; % index for the number to ask

% display initial information on the screen and record timing
LGCM_mental_effort_task_displayQuestion(scr, stim, startAngle, endAngle, task_trialStart, R_chosen, R_or_P, numberVector(jNumber));
[~,onsetEffortPeriod] = Screen(window,'Flip');

%% start the effort
effortPhaseOver = 0; % condition for when the trial ends (either time or
% number of subsequent correct answers)
iCorrect = 0; % index for number of subsequent correct answers provided
task_type_tmp = task_trialStart; % variable determining whether odd/even task or higher/lower than 5 task

while effortPhaseOver == 0
    %% check timing => task stops when maximal time achieved, 
    % even if the correct amount of answers has not been attained
    timeNow = GetSecs;
    if timeNow > onsetEffortPeriod + t_effort_max
        trial_failed = 1;
        effortPhaseOver = 1;
    end
    
    %% check whether total amount of subsequent correct answers has been reached or not
    if iCorrect == n_correctAnswersToGive
        trial_failed = 0; % trial is a success
        effortPhaseOver = 1; % effort phase is over
    end % end mental effort period if total amount of correct answers has been reached
    
    
    %% check key presses
    [keyisdown, timePress, keyCode] = KbCheck();
    
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
        [answerCorrect_tmp] = LGCM_mental_effort_answer_correct(task_type_tmp,...
            numberVector(jNumber),...
            sideAnswer_tmp, sideQuestion);
        switch answerCorrect_tmp
            case 0 % error made
                % re-initialize counter for number of subsequent correct
                % answers
                iCorrect = 0;
                % re-initialize the angle as well
                currentAngle = startAngle;
                % no change in 
                
            case 1 % correct
                % increase counter for number of correct answers
                iCorrect = iCorrect + 1;
                % adapt the angle
                currentAngle = startAngle + (iCorrect/n_correctAnswersToGive)*totalAngleDistance;
        end
        
        %% store performance
        %(when they answered and what they answered for each question)
        perf.taskType = [perf.taskType, task_type_tmp]; % which task
        perf.answerTime = [perf.answerTime, timePress]; % when did the subject answer
        perf.answer = [perf.answer, answerCorrect_tmp]; % was the answer correct?
        perf.sideAnswer = [perf.sideAnswer, sideAnswer_tmp]; % which side did he answer left or right?
        onsetAnswers(jNumber) = timePress;
        
        %% define task type for the next question
        if iCorrect >= 1
            % once a correct answer has been provided you need to assign the
            % moments when a switch will happen
            % = define a STAY/SWITCH sequence
            % (as long as errors are being made, keep the task easy)
            if iCorrect == 1
                [task_seq_tmp] = LGCM_mental_effort_task_switches(task_type_tmp, n_correctAnswersToGive, n_switch);
            end % first correct answer of a sequence
            
            % extract task type for the next number
            task_type_tmp = task_seq_tmp(iCorrect + 1);
        end % was the last answer correct?
        
        %% pass on to the next number of the list
        jNumber = jNumber + 1;
        
        %% display update of effort scale and of number value
        LGCM_mental_effort_task_displayQuestion(scr, stim, currentAngle, endAngle, task_type_tmp, R_chosen, R_or_P, numberVector(jNumber));
        [~,onsetNumbers(jNumber)] = Screen(window,'Flip');
    end % left or right key has been pressed
    
end % while effort period not finished

end % function