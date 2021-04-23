function[stim] = LGCM_stim_initialize(scr, n_R_levels, n_E_levels, pics_folder)
%[stim] = LGCM_stim_initialize(scr, n_R_levels, n_E_levels, pics_folder)
%LGCM_stim_initialize will initialize the v
%
% INPUTS
% scr: structure with main screen informations (size, center, window, etc.
%
% n_R_levels: number of reward levels
%
% n_E_levels: number of difficulty levels
%
% pics_folder: path where pictures of reward levels are stored
%
% OUTPUTS
% stim: structure with stimulus informations
%
% See also LGCM_main_experiment.m

%% extract screen main informations
window          = scr.window;
xScreenCenter   = scr.xCenter;
yScreenCenter   = scr.yCenter;
% screenXpixels   = xScreenCenter*2;
screenYpixels   = yScreenCenter*2;

% colours
black = [0 0 0];
% difficultyArcColor = [178 24 43];
difficultyArcColor = [255 210 0];
% white = [255 255 255];
% screen_background_colour = scr.background_colour;

%% Enable alpha blending for anti-aliasing of all our textures
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Money variables
moneySize = yScreenCenter/2; % 214 in initial Arthur version
stim.reward.moneyRect  = [0 0 moneySize moneySize];
stim.reward.moneySize  = moneySize;

% Import the reward images and store them into memory
punishmentRecalibrateCoord = 8;
for iR = 1:n_R_levels
    % extract name for subfield of the current difficulty level
    R_level_nm = ['reward_',num2str(iR)];
    
    % import the reward image
    [image_tmp, ~, alpha_tmp] = imread([pics_folder 'SMT_Money_',num2str(iR),'.png']);
    
    % transform into texture
    image_tmp(:,:,4) = alpha_tmp;
    
    % store the image
    stim_tmp = Screen('MakeTexture', window, image_tmp);
    stim.reward.texture.(R_level_nm) = stim_tmp;
    % display on top of the screen (for choice)
    stim.reward.top_center.(R_level_nm) = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter, yScreenCenter/2);
    stim.reward.top_left.(R_level_nm) = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter/2, yScreenCenter/2);
    stim.reward.top_right.(R_level_nm) = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter*(3/2), yScreenCenter/2);
    % display on middle of the screen for performance
    stim.reward.middle_center.(R_level_nm) = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter, yScreenCenter);
    
    % corresponding coordinates for red circle overlay for punishments
    stim.punishment.circleOverlay.top_center.(R_level_nm)       = [stim.reward.top_center.(R_level_nm)(1)+punishmentRecalibrateCoord,...
        stim.reward.top_center.(R_level_nm)(2)+punishmentRecalibrateCoord,...
        stim.reward.top_center.(R_level_nm)(3)-punishmentRecalibrateCoord,...
        stim.reward.top_center.(R_level_nm)(4)-punishmentRecalibrateCoord];
    stim.punishment.circleOverlay.top_left.(R_level_nm)         = [stim.reward.top_left.(R_level_nm)(1)+punishmentRecalibrateCoord,...
        stim.reward.top_left.(R_level_nm)(2)+punishmentRecalibrateCoord,...
        stim.reward.top_left.(R_level_nm)(3)-punishmentRecalibrateCoord,...
        stim.reward.top_left.(R_level_nm)(4)-punishmentRecalibrateCoord];
    stim.punishment.circleOverlay.top_right.(R_level_nm)        = [stim.reward.top_right.(R_level_nm)(1)+punishmentRecalibrateCoord,...
        stim.reward.top_right.(R_level_nm)(2)+punishmentRecalibrateCoord,...
        stim.reward.top_right.(R_level_nm)(3)-punishmentRecalibrateCoord,...
        stim.reward.top_right.(R_level_nm)(4)-punishmentRecalibrateCoord];
    stim.punishment.circleOverlay.middle_center.(R_level_nm)    = [stim.reward.middle_center.(R_level_nm)(1)+punishmentRecalibrateCoord,...
        stim.reward.middle_center.(R_level_nm)(2)+punishmentRecalibrateCoord,...
        stim.reward.middle_center.(R_level_nm)(3)-punishmentRecalibrateCoord,...
        stim.reward.middle_center.(R_level_nm)(4)-punishmentRecalibrateCoord];
end

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
stim.textRectSize.xSizeWin = textSizeWin(4) - textSizeWin(1);
[~,~,textSizeLose] = DrawFormattedText(window,'perdre',xScreenCenter,yScreenCenter,white);
stim.textRectSize.xSizeLose = textSizeLose(4) - textSizeLose(1);
[~,~,textSizeForEffort] = DrawFormattedText(window,'pour',xScreenCenter,yScreenCenter,white);
stim.textRectSize.xSizeForEffort = textSizeForEffort(4) - textSizeForEffort(1);

%% color used to represent the signal
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