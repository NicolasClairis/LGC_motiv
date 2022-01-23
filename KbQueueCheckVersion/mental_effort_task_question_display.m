function mental_effort_task_question_display(scr, stim, task_trialType, sideQuestion, textCol, learning_instructions)
% [] = mental_effort_task_question_display(scr, stim, task_trialStart, sideQuestion, textCol, learning_instructions)
% mental_effort_task_question_display will display question for the
% mental effort being made according to the current task.
%
% INPUTS
% scr: structure with screen information (window and x,y coordinates of the
% middle of the screen
%
% stim: structure with stimuli informations
%
% task_trialType
% (1) lower/higher than 5 task
% (2) whatever button (for first trials of the N-back task)
%
% sideQuestion: structure with side for each answer (especially if you decide to vary it)
% sideQuestion.hL.low, sideQuestion.hL.high: should i press left or right
% to indicate if the number is higher or lower than 5?
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
        
    case 2 % first trials for N-back version (where no correct answer needed)
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