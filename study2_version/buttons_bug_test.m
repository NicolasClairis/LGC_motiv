
%% Clear the workspace and the screen, instrreset resets the udp channels
ShowCursor;
sca; % close all PTB screens
close all; % close all windows
clearvars; % clear variables from memory
instrreset; % Disconnect and delete all instrument objects
clc;

%% initialize screen
IRM = 1;
testing_script = 1;
[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(IRM, testing_script);
ShowCursor;

%% initialize buttons
% how many possible answers
n_buttonsChoice = 4;
key = relevant_key_definition('mental', IRM, n_buttonsChoice);

%% press space before starting
if IRM == 1
    DrawFormattedText(window,'Press space to start.',...
        'center','center',scr.colours.white, scr.wrapat);
    Screen(window,'Flip');
    [~, ~, keyCode] = KbCheck();
    while(keyCode(key.space) ~= 1)
        % wait until the key has been pressed
        [~, ~, keyCode] = KbCheck();
    end
    DrawFormattedText(window,'Ok space was pressed.',...
        'center','center',scr.colours.white, scr.wrapat);
    Screen(window,'Flip');
    WaitSecs(1);
end

%% start checking keyboard presses
keyboard_check_start(key, IRM);

%% use the buttons
buttonTest(key, window);
% THIS FIRST CALL IS WHEN THE BUG GENERALLY HAPPENS

%%
DrawFormattedText(window,'Press space again.',...
    'center','center',scr.colours.white, scr.wrapat);
Screen(window, 'Flip');
[~, ~, keyCode] = KbCheck();
while(keyCode(key.space) ~= 1)
    % wait until the key has been pressed
    [~, ~, keyCode] = KbCheck();
end
DrawFormattedText(window,'Ok space was pressed.',...
    'center','center',scr.colours.white, scr.wrapat);
Screen(window,'Flip');
WaitSecs(1);

%% use the buttons again
buttonTest(key, window);

%% release key buffer
KbQueueStop;
KbQueueRelease;

%% Close the PTB screen
sca;