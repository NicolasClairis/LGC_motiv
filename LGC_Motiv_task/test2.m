[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(0, 1);

[trainingTimes, calibTimes, learningTimes, taskTimes, mainTimes] = timings_definition({'RP'}, 3, 3, 54, 'mental');
langage = 'fr';sca;

[stim] = stim_initialize(scr, 3, langage);
key = relevant_key_definition('mental', 0, 4);
mentalE_prm_calib = mental_effort_parameters();
mentalE_prm_calib.startAngle = 0; % for learning always start at zero
% no error threshold nor mapping of answers when errors are
% made
calib_errorLimits_Em.useOfErrorMapping = false;
calib_errorLimits_Em.useOfErrorThreshold = true;
calib_errorLimits_Em.errorThreshold = 5;

n_calibTrials = 2;
n_calibMax = 3;
n_errorsThreshold = 0;

[numberVector_calib] = mental_calibNumberVector(n_calibTrials, n_calibMax);
[n_mental_max_perTrial, calib_summary] = mental_calibNumbers(scr, stim, key,...
    numberVector_calib, mentalE_prm_calib, n_calibTrials, calibTimes, langage, n_errorsThreshold);
WaitSecs(3);
sca;