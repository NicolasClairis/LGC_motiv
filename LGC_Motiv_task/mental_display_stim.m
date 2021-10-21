function[onset_stim] = mental_display_stim(scr, stim,...
    startAngle, endAngle,...
    sideQuestion, taskTypeDisplay, taskTypePerf, numberValue, mental_n_col,...
    learning_instructions, maxPerfUntilNowAngle)
% [onset_stim] = mental_display_stim(scr, stim,...
%     startAngle, endAngle,...
%     sideQuestion, taskTypeDisplay, taskTypePerf, numberValue, mental_n_col,...
%     learning_instructions, maxPerfUntilNowAngle)
% mental_display_stim will display the arc, number to solve,
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
% taskTypeDisplay: (0) odd/even; (1): lower/higher than 5; (2) last question (for N-back version can be ignored = 0 in white)
% (type of the question displayed on screen which determines the number colour)
%
% taskTypePerf: (0) odd/even; (1): lower/higher than 5; (2) first question: any button press
% (type of the question to answer currently independent of the current display)
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
% maxPerfUntilNowAngle: for calibration, add an orange bar where the
% maximum perf has been reached until now
%
% OUTPUTS
% onset_stim: time when everything appears on screen
%
% See also mental_learning.m

%% extract relevant parameters
window = scr.window;
arcCurrLevelColor = stim.difficulty.currLevelColor;
arcPosition = stim.difficulty.middle_center;

%% percentage of correct answers already provided
Screen('FillArc', window,...
    arcCurrLevelColor,...
    arcPosition,...
    startAngle,...
    endAngle - startAngle);

% %%
% if taskTypePerf == 2 % for first question, add a short text to precise that any button press is ok
%     DrawFormattedText(window, 'Appuyez pour commencer.',...
%             'center', yScreenCenter*15/8, selectedCol);
% end

%% number to solve
switch taskTypeDisplay
    case 0 % odd/even
        textColor = mental_n_col.oddEven;
    case 1 % lower/higher than 5
        textColor = mental_n_col.lowHigh;
    case 2 % last question
        textColor = mental_n_col.lastQuestion;
end
% increase text size for number
Screen('TextSize', window, scr.textSize.mentalNumber);
% display number on screen
DrawFormattedText(window, num2str(numberValue),...
    stim.Em.(['numberPerf_',num2str(numberValue)]).x,...
    stim.Em.(['numberPerf_',num2str(numberValue)]).y,...
    textColor);
% text size back to baseline
Screen('TextSize', window, scr.textSize.baseline);

%% display orange bar where is the best performance until now
if ~isempty(maxPerfUntilNowAngle)
    Screen('FillArc', window,...
        arcCurrLevelColor,...
        arcPosition,...
        maxPerfUntilNowAngle,...
        maxPerfUntilNowAngle);
    error('find how to draw a line correct place')
end

%% instructions
switch learning_instructions
    case {'fullInstructions','partialInstructions'}
        mental_effort_task_question_display(scr, stim, taskTypeDisplay, sideQuestion, textColor, learning_instructions);
end

%% display on screen
[~,onset_stim] = Screen(window,'Flip');

end % function