%% script to display chosen option for describin experiment
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
[stim] = stim_initialize(scr, n_E_levels, langage);

n_buttonsChoice = 4;
IRMbuttons = 1;
[key] = relevant_key_definition(effort_type, IRMbuttons, n_buttonsChoice);
confidence.display = true;
confidence.lowOrHigh = 1;

R_chosen = 0.50; 
E_chosen = 0;
R_or_P = 'P';
choice_task_dispChosen(scr, stim,  -1, R_chosen, E_chosen, R_or_P, confidence);
img_array = Screen('GetImage',window);
imwrite(img_array, 'dispChosen.jpg');
imwrite(img_array, 'dispChosen.png');
WaitSecs(1);
KbWait;
sca;