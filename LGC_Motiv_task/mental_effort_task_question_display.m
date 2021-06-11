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

%% display on the screen
switch task_trialType
    case 0 % odd/even
        if strcmp(learning_instructions,'fullInstructions')
            DrawFormattedText(window, stim.Em.oddORevenQuestion.text,...
                stim.Em.oddORevenQuestion.x, stim.Em.oddORevenQuestion.y, textCol);
        end
        
        if sideQuestion.oE.pair == -1 && sideQuestion.oE.impair == +1
            x_pair      = stim.Em.even_left.x;
            x_impair    = stim.Em.odd_right.x;
        elseif sideQuestion.oE.pair == +1 && sideQuestion.oE.impair == -1
            x_pair      = stim.Em.even_right.x;
            x_impair    = stim.Em.odd_left.x;
        else
            error('error in sideQuestion definition');
        end
        DrawFormattedText(window, stim.Em.even.text, x_pair, stim.Em.even.y, textCol );     % pair
        DrawFormattedText(window, stim.Em.OR.text, stim.Em.OR.x, stim.Em.OR.y, textCol );   % OR
        DrawFormattedText(window, stim.Em.odd.text, x_impair, stim.Em.odd.y, textCol );     % impair

    case 1 % higher/lower than 5?
        if strcmp(learning_instructions,'fullInstructions')
            DrawFormattedText(window, stim.Em.lowerORhigherQuestion.text,...
                stim.Em.lowerORhigherQuestion.x, stim.Em.lowerORhigherQuestion.y, textCol);
        end
        
        if sideQuestion.hL.low == -1 && sideQuestion.hL.high == +1
            x_low = stim.Em.lower_left.x;
            x_high = stim.Em.higher_right.x;
        elseif sideQuestion.hL.low == +1 && sideQuestion.hL.high == -1
            x_low = stim.Em.lower_right.x;
            x_high = stim.Em.higher_left.x;
        else
            error('error in sideQuestion definition');
        end
        DrawFormattedText(window, stim.Em.lower.text, x_low, stim.Em.lower.y, textCol );    % < 5
        DrawFormattedText(window, stim.Em.OR.text, stim.Em.OR.x, stim.Em.OR.y, textCol );   % OR
        DrawFormattedText(window, stim.Em.higher.text, x_high, stim.Em.higher.y, textCol ); % > 5
        
    case 2 % first trial for 1-back version
        if strcmp(learning_instructions,'fullInstructions')
            DrawFormattedText(window, stim.Em.pressAnyButtonQuestion.text,...
                stim.Em.pressAnyButtonQuestion.x, stim.Em.pressAnyButtonQuestion.y, textCol);
%             DrawFormattedText(window, 'Press any button',...
%                 'center', yScreenCenter/3, textCol);
        end

        DrawFormattedText(window, stim.Em.pressAnyButton.text, stim.Em.pressAnyButton_left.x, stim.Em.pressAnyButtonQuestion.y, textCol );
        DrawFormattedText(window, stim.Em.OR.text, stim.Em.OR.x, stim.Em.OR.y, textCol );   % OR
        DrawFormattedText(window, stim.Em.pressAnyButton.text, stim.Em.pressAnyButton_right.x, stim.Em.pressAnyButtonQuestion.y, textCol );
end % task type
    
end % function