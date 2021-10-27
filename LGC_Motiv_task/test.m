[scr, xScreenCenter, yScreenCenter, window, baselineTextSize] = ScreenConfiguration(0, 1);
ShowCursor;
langage = 'fr';
n_E_levels = 3;
[stim] = stim_initialize(scr, n_E_levels, langage);
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

startAngle = 0;
endAngle = 360;
maxPerfUntilNowAngle = 270;


[onset_stim] = mental_display_stim(scr, stim,...
    startAngle, endAngle,...
    sideQuestion, taskTypeDisplay, numberVector_calib, mental_n_col,...
    'noInstructions', maxPerfUntilNowAngle);
