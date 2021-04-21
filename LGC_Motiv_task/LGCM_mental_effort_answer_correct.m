function[answerCorrect] = LGCM_mental_effort_answer_correct(task_type, numberVal, sideAnswer, sideQuestion)
%[answerCorrect] = LGCM_mental_effort_answer_correct(task_type, numberVal, sideAnswer, sideQuestions)
% LGCM_mental_effort_answer_correct will tell you if the current question
% was solved correctly or not
%
% INPUTS
%
% task_type
% (0) odd/even task
% (1) lower/higher than 5 task
%
% numberVal: number value for the current question
%
% sideAnswer
% (-1): answered left
% (+1): answered right
%
% sideQuestion: structure with side for each answer (especially if you decide to vary it)
% sideQuestion.oE.pair, sideQuestion.oE.impair, sideQuestion.hL.low,
% sideQuestion.hL.high
% (-1) left
% (+1) right
%
% OUTPUTS
%
% answerCorrect
% (0) answer provided was wrong
% (1) answer provided was correct

switch task_type
    case 0 % odd/even task
        
        if ( (sideAnswer == sideQuestion.oE.pair) && (mod(numberVal,2) == 0) ) ||...
               (sideAnswer == sideQuestion.oE.impair) && (mod(numberVal,2) ~= 0 ) % select pair or impair and was correct
            answerCorrect = 1;
        else
            answerCorrect = 0;
        end % answer correct or not?
        
    case 1 % lower/higher than 5 task
        
        if ( (sideAnswer == sideQuestion.hL.low) && (numberVal < 5) ) ||...
               (sideAnswer == sideQuestion.hL.high) && (numberVal > 5) % select lower than 5 and was correct or selected higher than 5 and was correct
            answerCorrect = 1;
        else
            answerCorrect = 0;
        end % answer correct or not?
        
end % task type switch

end % function