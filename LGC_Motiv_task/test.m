%% script to find where to put the numbers for the display of performance

[scr, xScreenCenter, yScreenCenter, window, baselineTextSize] = ScreenConfiguration(0, 1);
ShowCursor;
langage = 'fr';
n_E_levels = 4;
[key] = relevant_key_definition('mental', 0, 2);
[mentalE_prm] = mental_effort_parameters();
mentalE_prm.startAngle = 0;
sideQuestion = mentalE_prm.sideQuestion;
taskTypeDisplay = 1;
numberVector_calib = 1;
mental_n_col = mentalE_prm.mental_n_col;
calib_errorLimits_Em.useOfErrorMapping = false;
calib_errorLimits_Em.useOfErrorThreshold = true;
calib_errorLimits_Em.errorThreshold = 3;
white = scr.colours.white;

startAngle = 0;
endAngle = 360;
maxPerfUntilNowAngle = 270;

Rtop = 0.80;
Rmin = 0.00;

Ptop = -1.60;
Pmin = -0.80;

% extract coordinates
[stim] = stim_initialize(scr, n_E_levels, langage);

% display numbers 
DrawFormattedText(window, '+0.56',stim.leftMoneyWinEperf.x, stim.leftMoneyWinEperf.y, white);
DrawFormattedText(window, '+0.00',stim.rightMoneyWinEperf.x, stim.rightMoneyWinEperf.y, white);
% display the rest
[onset_stim] = mental_display_stim(scr, stim,...
    startAngle, endAngle,...
    sideQuestion, taskTypeDisplay, numberVector_calib, mental_n_col,...
    'noInstructions', maxPerfUntilNowAngle);