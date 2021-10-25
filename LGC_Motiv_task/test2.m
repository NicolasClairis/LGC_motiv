[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(0, 1);

[trainingTimes, calibTimes, learningTimes, taskTimes, mainTimes] = timings_definition({'RP'}, 3, 3, 54, 'mental');
langage = 'fr';sca;

[R_money] = R_amounts(3, 'yes');
[stim] = stim_initialize(scr, 3, langage, R_money);
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

[numberVector_calib] = mental_calibNumberVector(n_calibTrials, n_calibMax);
[n_mental_max_perTrial, calib_summary] = mental_calibNumbers(scr, stim, key,...
    numberVector_calib, mentalE_prm_calib, n_calibTrials, calibTimes);
WaitSecs(3);
sca;