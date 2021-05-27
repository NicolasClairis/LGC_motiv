function[stim] = stim_initialize(scr, n_R_levels, n_E_levels)
%[stim] = stim_initialize(scr, n_R_levels, n_E_levels)
%stim_initialize will initialize the v
%
% INPUTS
% scr: structure with main screen informations (size, center, window, etc.
%
% n_R_levels: number of reward levels
%
% n_E_levels: number of difficulty levels
%
% OUTPUTS
% stim: structure with stimulus informations
%
% See also main_experiment.m

%% extract screen main informations
window          = scr.window;
xScreenCenter   = scr.xCenter;
yScreenCenter   = scr.yCenter;
% screenXpixels   = xScreenCenter*2;
screenYpixels   = yScreenCenter*2;

% colours
black = scr.colours.black;
white = scr.colours.white;
grey = scr.colours.grey;
% difficultyArcColor = [178 24 43];
difficultyArcColor = [255 210 0];
% white = [255 255 255];
% screen_background_colour = scr.background_colour;

%% Enable alpha blending for anti-aliasing of all our textures
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Money variables

% extract reward amount text size (for choice)
[~,~,textSizeR] = DrawFormattedText(window,'+0.00 CHF', xScreenCenter, yScreenCenter, white);
xSizeText = textSizeR(3) - textSizeR(1);
ySizeText = textSizeR(4) - textSizeR(2);
stim.reward.xSizeText = xSizeText;
stim.reward.ySizeText = ySizeText;
% define where the text will be displayed
stim.reward.text.top_center_start = [xScreenCenter - xSizeText/2, yScreenCenter*(2/3) - ySizeText/2]; % for display of chosen option
stim.reward.text.top_left_start = [xScreenCenter/2 - xSizeText/2, yScreenCenter/2 - ySizeText/2]; % for display of left option
stim.reward.text.top_right_start = [xScreenCenter*(3/2) - xSizeText/2, yScreenCenter/2 - ySizeText/2]; % for display of right option
% display on middle of the screen for performance feedback
stim.reward.text.middle_center_start = [xScreenCenter - xSizeText/2, yScreenCenter - ySizeText/2];

% add grey screen on top to be sure that this does not actually appear on
% the screen
Screen('FillRect',window, grey, [0 0 xScreenCenter*2 yScreenCenter*2]);
Screen(window,'Flip');

%% difficulty rings
% 
% position each ring on the screen (for choice task)
stim.difficulty.below_center = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter, yScreenCenter*(3/2));
stim.difficulty.below_left = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter/2, yScreenCenter*(3/2));
stim.difficulty.below_right = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter*(3/2), yScreenCenter*(3/2));
% position each ring on the screen (for performance task)
stim.difficulty.middle_center = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter, yScreenCenter);
stim.difficulty.arcEndAngle = 360;

% define the circle size for each difficulty level depending on the
% difficulty
% note level 1 = easiest level (ascending order)
for iDiff = 1:n_E_levels
    % extract name for subfield of the current difficulty level
    diff_level_nm = ['level_',num2str(iDiff)];
    
    % extract angle for the arc which will correspond to the difficulty
    % level: max circle = max difficulty level
    startAngle_tmp = stim.difficulty.arcEndAngle*((n_E_levels - iDiff)./n_E_levels);
    if startAngle_tmp < 360
        stim.difficulty.startAngle.(diff_level_nm) = startAngle_tmp;
    elseif startAngle_tmp == 360
        stim.difficulty.startAngle.(diff_level_nm) = 0;
    end
end % difficulty

%% extract text size
[~,~,textSizeWin] = DrawFormattedText(window,'Gagner',xScreenCenter,yScreenCenter,white);
stim.textRectSize.xSizeWin = textSizeWin(3) - textSizeWin(1);
stim.textRectSize.ySizeWin = textSizeWin(4) - textSizeWin(2);
[~,~,textSizeLose] = DrawFormattedText(window,'perdre',xScreenCenter,yScreenCenter,white);
stim.textRectSize.xSizeLose = textSizeLose(3) - textSizeLose(1);
stim.textRectSize.ySizeLose = textSizeLose(4) - textSizeLose(2);
[~,~,textSizeForEffort] = DrawFormattedText(window,'pour',xScreenCenter,yScreenCenter,white);
stim.textRectSize.xSizeForEffort = textSizeForEffort(3) - textSizeForEffort(1);
stim.textRectSize.ySizeForEffort = textSizeForEffort(4) - textSizeForEffort(2);
% add grey screen on top to be sure that this does not actually appear on
% the screen
Screen('FillRect',window, grey, [0 0 xScreenCenter*2 yScreenCenter*2]);
Screen(window,'Flip');

%% color used to represent the signal
% no use of monetary images anymore
alpha_punishment = 115;
stim.punishment.colourOverlay   = [255 50 0 alpha_punishment];
stim.difficulty.maxColor        = black;
stim.difficulty.currLevelColor  = difficultyArcColor;
stim.difficulty.ovalWidth       = 3;

%% square to display for the chosen option
stim.chosenOption.reward = stim.reward.top_center;
stim.chosenOption.difficulty = stim.difficulty.below_center;
stim.chosenOption.squareColour = black;
stim.chosenOption.squareRect = [xScreenCenter*(2/3),...
    yScreenCenter*(1/6),...
    xScreenCenter*(4/3),...
    yScreenCenter*(11/6)];
stim.chosenOption.squareWidth = 10;

%% force display

%% for the end of the performance period
% circle around the money to signify end of the trial (win or loss)
stim.endTrialcircle  = [0 0 moneySize+(moneySize/5) moneySize+(moneySize/5)];
stim.end_trial.middle_center = CenterRectOnPointd(stim.endTrialcircle, xScreenCenter, yScreenCenter);

%% define bar size for the waiting time
stim.barTimeWaitRect = [xScreenCenter*(1/2),...
    yScreenCenter*(3/4),...
    xScreenCenter*(3/2),...
    yScreenCenter];

%% fixation cross coordinates on the screen (code relative to screen Y size)
cross_length    = screenYpixels/6;
cross_thickness = 0.2*cross_length;
stim.cross.verticalLine = [xScreenCenter - (cross_thickness/2),...
    yScreenCenter - (cross_length/2),...
    xScreenCenter + (cross_thickness/2),...
    yScreenCenter + (cross_length/2)];
stim.cross.horizontalLine = [xScreenCenter - (cross_length/2),...
    yScreenCenter - (cross_thickness/2),...
    xScreenCenter + (cross_length/2),...
    yScreenCenter + (cross_thickness/2)];

end % function