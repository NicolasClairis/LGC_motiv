%% script to check that visual display is ok with the CIBM screen

IRM = 1;
%% load screen
[scr, xScreenCenter, yScreenCenter, window, baselineTextSize] = ScreenConfiguration(IRM, 0);
ShowCursor;


%% prepare the parameters
n_R_levels=4;
n_E_levels = 4;
langage = 'fr';
n_buttonsChoice=4;
confidenceDisp = 1;
[stim] = stim_initialize(scr, n_E_levels, langage);
key = relevant_key_definition('mental', IRM, n_buttonsChoice);
timeParameter.timeLimit = false;


%% start checking
keyboard_check_start(key, IRM);

%% display a choice trial on the screen (script will wait one of the 4 buttons to be pressed)
disp('please press one of the four buttons to move on');
choice_period(scr, stim,...
        1.6, 1.8, 1, 2, 'R',...
        timeParameter, key, confidenceDisp);

%% display mental effort on the screen to check if big enough
mentalE_prm = mental_effort_parameters();
[onsetTrial] = mental_display_stim(scr, stim,...
    0, 360,...
    mentalE_prm.sideQuestion, 1, 1, mentalE_prm.mental_n_col,...
    'noInstructions');
disp('please press space to move on'); 
[~, ~, keyCode] = KbCheck();
while(keyCode(key.space) ~= 1)
    % wait until the key has been pressed
    [~, ~, keyCode] = KbCheck();
end

%% shut down PTB
sca;
KbQueueStop;
KbQueueRelease;