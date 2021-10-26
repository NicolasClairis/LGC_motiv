% test choice period
IRM = 0;
testing_script = 1;
[scr, xScreenCenter, yScreenCenter, window, baselineTextSize] = ScreenConfiguration(IRM, testing_script);
ShowCursor;

n_R_levels = 4;
punishment_yn = 'yes';
IPdata.baselineR = 0.5;
IPdata.baselineP = 0.5;
IPdata.mentalDeltaIP = 0.3;
effort_type = 'mental';
[R_money] = R_amounts_IP(n_R_levels, punishment_yn, IPdata, effort_type);
langage = 'fr';
n_E_levels = 4;
[stim] = stim_initialize(scr, n_E_levels, langage, R_money);
R_left = R_money.R_2;
R_right = R_money.R_0;
E_left = 2;
E_right = 0;
R_or_P = 'R';
n_buttonsChoice = 4;
IRMbuttons = 1;
timeParameter.timeLimit = true;
timeParameter.t_choice = 5;
[key] = relevant_key_definition(effort_type, IRMbuttons, n_buttonsChoice);
confidenceDisp = false;

[choice_trial, onsetDispChoiceOptions, onsetChoice, stoptask] = choice_period(scr, stim,...
    R_left, R_right, E_left, E_right, R_or_P,...
    timeParameter, key, confidenceDisp);