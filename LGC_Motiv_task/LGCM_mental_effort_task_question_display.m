function[] = LGCM_mental_effort_task_question_display(scr, task_trialType, sideQuestion, textCol)
% [] = LGCM_mental_effort_task_question_display(scr, task_trialStart, sideQuestion, textCol)
% LGCM_mental_effort_task_question_display will display question for the
% mental effort being made according to the current task.
%
% INPUTS
% scr: structure with screen information (window and x,y coordinates of the
% middle of the screen
%
% task_trialType
% (0) odd/even task
% (1) lower/higher than 5 task
%
% sideQuestion: structure with side for each answer (especially if you decide to vary it)
% sideQuestion.oE.pair, sideQuestion.oE.impair, sideQuestion.hL.low,
% sideQuestion.hL.high
% (-1) left
% (+1) right
%
% textCol: rgb code to know with which colour should the text be displayed
%
% See also LGCM_mental_effort.m

%% extract relevant info
window = scr.window;
xScreenCenter = scr.xCenter;
yScreenCenter = scr.yCenter;
blackCol = [0 0 0];

%% determine x,y coordinates
y_coord = yScreenCenter*(5/3);
x_left = xScreenCenter/2;
x_right = xScreenCenter*(3/2);

%% display on the screen
switch task_trialType
        case 0 % odd/even
            DrawFormattedText(window, 'Chiffre pair ou impair?',...
                    'center', yScreenCenter/3, textCol);
                
            if sideQuestion.oE.pair == -1 && sideQuestion.oE.impair == +1
                    x_pair = x_left;
                    x_impair = x_right;
            elseif sideQuestion.oE.pair == +1 && sideQuestion.oE.impair == -1
                    x_pair = x_right;
                    x_impair = x_left;
            else
                error('error in sideQuestion definition');
            end
            DrawFormattedText(window,'pair', x_pair, y_coord, textCol );
            DrawFormattedText(window,'OU', 'center', y_coord, blackCol );
            DrawFormattedText(window,'impair', x_impair, y_coord, textCol );
            
        case 1 % higher/lower than 5?
            DrawFormattedText(window, 'Chiffre < ou > 5?',...
                    'center', yScreenCenter/3, textCol);
                
            if sideQuestion.hL.low == -1 && sideQuestion.hL.high == +1
                    x_low = x_left;
                    x_high = x_right;
            elseif sideQuestion.hL.low == +1 && sideQuestion.hL.high == -1
                    x_low = x_right;
                    x_high = x_left;
            else
                error('error in sideQuestion definition');
            end
            DrawFormattedText(window,'< 5', x_low, y_coord, textCol );
            DrawFormattedText(window,'OU', 'center', y_coord, blackCol );
            DrawFormattedText(window,'> 5', x_high, y_coord, textCol );
end % task type
    
end % function