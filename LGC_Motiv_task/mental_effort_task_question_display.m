function mental_effort_task_question_display(scr, task_trialType, sideQuestion, textCol, learning_instructions)
% [] = mental_effort_task_question_display(scr, task_trialStart, sideQuestion, textCol, learning_instructions)
% mental_effort_task_question_display will display question for the
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
% learning_instructions
% 'fullInstructions': display instructions: ask if odd/even (or lower/higher than 5) and
% display also on the screen the relevant answer to each question
% 'partialInstructions': display only the two possible answers but not the
% question anymore
% 'noInstructions': no reminder of what the question is nor of where you should answer
%
% See also mental_effort.m

%% check no error in script
if ~ismember(learning_instructions,{'fullInstructions','partialInstructions'})
    error(['learning instructions should be equal to fullInstructions or partialInstructions',...
        ' but currently equal to ',learning_instructions,'. Please fix it']);
end

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
        if strcmp(learning_instructions,'fullInstructions')
            DrawFormattedText(window, 'Chiffre pair ou impair?',...
                'center', yScreenCenter/3, textCol);

        end
        
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
        DrawFormattedText(window,'OU', 'center', y_coord, textCol );
        DrawFormattedText(window,'impair', x_impair, y_coord, textCol );

    case 1 % higher/lower than 5?
        if strcmp(learning_instructions,'fullInstructions')
            DrawFormattedText(window, 'Chiffre < ou > 5?',...
                'center', yScreenCenter/3, textCol);
        end
        
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
        DrawFormattedText(window,'OU', 'center', y_coord, textCol );
        DrawFormattedText(window,'> 5', x_high, y_coord, textCol );
        
    case 2 % first trial for 1-back version
        if strcmp(learning_instructions,'fullInstructions')
            DrawFormattedText(window, 'Appuyer sur n''importe quel bouton',...
                'center', yScreenCenter/3, textCol);
%             DrawFormattedText(window, 'Press any button',...
%                 'center', yScreenCenter/3, textCol);
        end

        DrawFormattedText(window,'Appuyer', x_left, y_coord, textCol );
        DrawFormattedText(window,'OU', 'center', y_coord, textCol );
        DrawFormattedText(window,'Appuyer', x_right, y_coord, textCol );
end % task type
    
end % function