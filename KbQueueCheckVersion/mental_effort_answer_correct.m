function[answerCorrect] = mental_effort_answer_correct(task_type, numberVal, sideAnswer, sideQuestion)
%[answerCorrect] = mental_effort_answer_correct(task_type, numberVal, sideAnswer, sideQuestions)
% mental_effort_answer_correct will tell you if the current question
% was solved correctly or not
%
% INPUTS
%
% task_type
% (1) lower/higher than 5 task
% (2) whatever answer is ok
%
% numberVal: number value for the current question
%
% sideAnswer
% (-1): answered left
% (+1): answered right
%
% sideQuestion: structure with side for each answer (especially if you decide to vary it)
% sideQuestion.hL.low, sideQuestion.hL.high
% (-1) left
% (+1) right
%
% OUTPUTS
%
% answerCorrect
% (0) answer provided was wrong
% (1) answer provided was correct

switch task_type
    case 1 % lower/higher than 5 task
        
        if ( (sideAnswer == sideQuestion.hL.low) && (numberVal < 5) ) ||...
               (sideAnswer == sideQuestion.hL.high) && (numberVal > 5) % select lower than 5 and was correct or selected higher than 5 and was correct
            answerCorrect = 1;
        else
            answerCorrect = 0;
        end % answer correct or not?
        
    case 2 % for N-back version first questions (press any button is fine)
        if ismember(sideAnswer,[-1,1])
            answerCorrect = 1;
        end
end % task type switch

end % function