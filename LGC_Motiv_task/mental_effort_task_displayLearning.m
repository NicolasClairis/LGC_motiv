function mental_effort_task_displayLearning(scr, stim, startAngle, endAngle, task_type, numberValue)
%mental_effort_task_displayQuestion(scr, stim, startAngle, endAngle, task_type, R_chosen, R_or_P, numberValue)
% mental_effort_task_displayQuestion will display all relevant information for mental
% effort task
%
% INPUTS
% scr: structure with screen informations
%
% stim: structure with stimuli display relevant information
%
% startAngle: initial angle of the effort arc
%
% endAngle: end angle of the effort arc
%
% task_type
% (0) odd/even task
% (1) lower/higher than 5 task
%
% numberValue: number for the current question

%% extract initial information
% screen informations
window = scr.window;
xScreenCenter = scr.xScreenCenter;
yScreenCenter = scr.yScreenCenter;

white = scr.colours.white;
%% display amount of effort to be done on top of the incentive as an arc
Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
    stim.difficulty.middle_center,...
    startAngle,...
    endAngle - startAngle);

% display number to solve
DrawFormattedText(window,num2str(numberValue), xScreenCenter, yScreenCenter*(1/6),white);

% display question according to type of task on which to start
mental_effort_task_question_display(scr, task_type, sideQuestion);

end % function