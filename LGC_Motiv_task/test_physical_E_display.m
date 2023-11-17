% script to show different effort levels in the physical effort task

[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(0, 1);

[trainingTimes, calibTimes, learningTimes, taskTimes, mainTimes] = timings_definition({'RP'}, 54, 4, 'mental');
langage = 'fr';

n_E_levels = 4;
[stim] = stim_initialize(scr, n_E_levels, langage);
key = relevant_key_definition('mental', 0, 4);

F_threshold = 55; % force should be maintained above this threshold (expressed in % of MVC)
F_tolerance = 2.5; % tolerance allowed around the threshold (expressed in % of MVC)

%% main variables
E_chosen = 2;
R_or_P = 'R';
R_chosen = 0.65;
savePath = 'P:\boulot\postdoc_CarmenSandi\experiment_figures\task\';
%% angle values
startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen)]);
currentAngle = startAngle;
endAngle = stim.difficulty.arcEndAngle;
totalAngleDistance = endAngle - startAngle;

%% display effort scale on the left for the live display of the effort level
% add force threshold to maintain
force_initialDisplay = 55;
disp_realtime_force(scr, F_threshold, F_tolerance, force_initialDisplay, 'task');

%% display difficulty level on top of the reward as an arc
Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
    stim.difficulty.middle_center,...
    startAngle,...
    endAngle - startAngle);

% add levels of money
drawMoneyProportional(scr, stim, R_chosen, R_or_P);

Screen('Flip',window);

imageArray = Screen('GetImage', window);
imwrite(imageArray,[savePath,'test_physical_E_display_E',num2str(E_chosen),'.png']);

WaitSecs(1);
sca;