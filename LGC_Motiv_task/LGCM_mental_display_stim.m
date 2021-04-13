function[onset_stim] = LGCM_mental_display_stim(scr, stim,...
    startAngle, endAngle,...
    sideQuestion, taskType, numberValue, mental_n_col,...
    learning_instructions)
% [onset_stim] = LGCM_mental_display_stim(scr, stim,...
%     startAngle, endAngle,...
%     sideQuestion, taskType, numberValue, mental_n_col,...
%     learning_instructions)
% LGCM_mental_display_stim will display the arc, number to solve,
% instructions and reward level (all relevant info) according to the inputs
%
% INPUTS
% scr: screen structure with relevant info about screen
%
% stim: structure with relevant info about stimulus display
%
% startAngle, endAngle: values for the angle of the arc showing how far you
% are from reaching the end
%
% sideQuestion: structure explaining which side corresponds to which answer
% for each task (odd/even; lower/higher than 5 vs left/right)
%
% taskType: (0) odd/even; (1): lower/higher than 5
%
% numberValue: number for the current question
%
% mental_n_col: structure with information about colour corresponding to
% each task type
%
% learning_instructions
% 'fullInstructions': display instructions: ask if odd/even (or lower/higher than 5) and
% display also on the screen the relevant answer to each question
% 'partialInstructions': display only the two possible answers but not the
% question anymore
% 'noInstructions': no reminder of what the question is nor of where you should answer
%
% OUTPUTS
% onset_stim: time when everything appears on screen
%
% See also LGCM_mental_learning.m

%% extract relevant parameters
window = scr.window;
xScreenCenter = scr.xCenter;
yScreenCenter = scr.yCenter;
arcCurrLevelColor = stim.difficulty.currLevelColor;
arcPosition = stim.difficulty.middle_center;

%% percentage of correct answers already provided
Screen('FillArc', window,...
    arcCurrLevelColor,...
    arcPosition,...
    startAngle,...
    endAngle - startAngle);

%% number to solve
switch taskType
    case 0 % odd/even
        textColor = mental_n_col.oddEven;
    case 1 % lower/higher than 5
        textColor = mental_n_col.lowHigh;
end
DrawFormattedText(window, num2str(numberValue),...
    xScreenCenter, yScreenCenter*(9/6), textColor);

%% instructions
switch learning_instructions
    case {'fullInstructions','partialInstructions'}
        LGCM_mental_effort_task_question_display(scr, taskType, sideQuestion, textColor, learning_instructions);
end

%% reward level

%% display on screen
[~,onset_stim] = Screen(window,'Flip');

end % function