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
blue = [0 255 0];
% white = [255 255 255];
% screen_background_colour = scr.background_colour;

%% Enable alpha blending for anti-aliasing of all our textures
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Money variables
moneySize = yScreenCenter/3; % 214 in initial Arthur version
stim.reward.moneyRect  = [0 0 moneySize moneySize];
stim.reward.moneySize  = moneySize;

% Import the reward images and store them into memory
for iR = 0:n_R_levels
    % extract name for subfield of the current difficulty level
    R_level_nm = ['reward_',num2str(iR)];
    
    % import the reward image
    [image_tmp, ~, alpha_tmp] = imread([pics_folder 'SMT_Money_',num2str(iR),'.png']);
    
    % transform into texture
    image_tmp(:,:,4) = alpha_tmp;
    
    % store the image
    stim_tmp = Screen('MakeTexture', window, image_tmp);
    stim.reward.texture.(R_level_nm) = stim_tmp;
    stim.reward.middle.(R_level_nm) = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter, yScreenCenter);
    stim.reward.left.(R_level_nm) = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter/2, yScreenCenter);
    stim.reward.right.(R_level_nm) = CenterRectOnPointd(stim.reward.moneyRect, xScreenCenter*(3/2), yScreenCenter);
end

%% size of the different circles
% threshold of the maximal duration circle
max_dur_circle = [0 0 yScreenCenter yScreenCenter];
% define the circle size for each difficulty level depending on the
% distance between the money on the screen and the maximal duration
distBtwMoneyAndMaxDur = max_dur_circle(4) - moneySize;
circleForEachDifficultyLevel = [0 0 distBtwMoneyAndMaxDur distBtwMoneyAndMaxDur];
% define the circle size for each difficulty level depending on the di
for iDiff = 1:n_E_levels
    % extract name for subfield of the current difficulty level
    diff_level_nm = ['level_',num2str(iDiff)];
    
    if iDiff == n_E_levels
        diff_circle.(diff_level_nm) = max_dur_circle;
    elseif iDiff < n_E_levels
        diff_circle.(diff_level_nm) = (iDiff/n_E_levels).*circleForEachDifficultyLevel;
    end
    
    % position each ring on the screen
    stim.difficulty.center.(diff_level_nm) = CenterRectOnPointd(diff_circle.(diff_level_nm), xScreenCenter, yScreenCenter);
    stim.difficulty.left.(diff_level_nm) = CenterRectOnPointd(diff_circle.(diff_level_nm), xScreenCenter/2, yScreenCenter);
    stim.difficulty.right.(diff_level_nm) = CenterRectOnPointd(diff_circle.(diff_level_nm), xScreenCenter*(3/2), yScreenCenter);
end % difficulty

%% color used to represent the signal
stim.difficulty.maxColor        = black;
stim.difficulty.currLevelColor  = blue;
stim.difficulty.ovalWidth       = 3;

%% force display


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