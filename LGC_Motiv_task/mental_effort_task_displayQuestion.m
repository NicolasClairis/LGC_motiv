function[] = LGCM_mental_effort_task_displayQuestion(scr, stim, startAngle, endAngle, task_type, R_chosen, R_or_P, numberValue)
%LGCM_mental_effort_task_displayQuestion(scr, stim, startAngle, endAngle, task_type, R_chosen, R_or_P, numberValue)
% LGCM_mental_effort_task_displayQuestion will display all relevant information for mental
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
% R_chosen: Reward level for this trial
%
% R_or_P:
% 'R': rewarding trial
% 'P': punishing trial
%
% numberValue: number for the current question

%% extract initial information
% screen informations
window = scr.window;
xScreenCenter = scr.xScreenCenter;
yScreenCenter = scr.yScreenCenter;
% reward/punishment display informations
R_img = stim.reward.texture.(['reward_',num2str(R_chosen)]); % reward image
R_coords = stim.reward.middle_center.(['reward_',num2str(R_chosen)]); % reward coordinates

%% display reward (or punishment) display on the background
switch R_or_P
    case 'R'
        Screen('DrawTexture', window,...
            R_img,...
            [],...
            R_coords);
    case 'P'
        error('punishment background to be defined yet');
end

%% display amount of effort to be done on top of the incentive as an arc
Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
    stim.difficulty.middle_center,...
    startAngle,...
    endAngle - startAngle);

% display number to solve
DrawFormattedText(window,num2str(numberValue), 'center', yScreenCenter*(1/6));

% display question according to type of task on which to start
LGCM_mental_effort_task_question_display(scr, task_type);

end % function