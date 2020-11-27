function[stim] = LGCM_stim_initialize(scr)
%[stim] = LGCM_stim_initialize(scr)
%LGCM_stim_initialize will initialize the v
%
% INPUTS
% scr: structure with main screen informations (size, center, window, etc.
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

%% Enable alpha blending for anti-aliasing of all our textures
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% color used to represent the signal
stim.signalColor = [0 0.75 0.3];

%% Money variables
moneySize = warning('should be encoded relatively to screen size'); %214;
stim.money = [0 0 moneySize moneySize];
stim.moneySize      = moneySize;

%% thresholds for reaching target
% 30% of MVC defined to this size so that 70 can be at 500 pixels
threshold_1    = [0 0 moneySize moneySize];
% 70% of MVC
threshold_2 = warning('should be encoded relatively to screen size');% [0 0 500 500];

%% missed position
%
missTargetRect     = [0 0 moneySize moneySize];
% 60% of MVC
missThresholdRect = warning('should be encoded relatively to screen size');%[0 0 429 429];

%% Position of each ring/coin/signal on the screen
% force thresholds
stim.centeredthreshold_1 = CenterRectOnPointd(threshold_1, xScreenCenter, yScreenCenter);
stim.centeredthreshold_2 = CenterRectOnPointd(threshold_2, xScreenCenter, yScreenCenter);
% miss position on the screen
stim.missThreshold = CenterRectOnPointd(missThresholdRect, xScreenCenter, yScreenCenter);
stim.missTarget = CenterRectOnPointd(missTargetRect, xScreenCenter*2 - moneySize/2, yScreenCenter);

%% Import the images
[image_20Cent, ~, alpha_20Cent]   = imread([pics_folder 'SMT_Coin_20RP.png']);
[image_50Cent, ~, alpha_50Cent]   = imread([pics_folder 'SMT_Coin_50RP.png']);
[image_1FR, ~, alpha_1FR]         = imread([pics_folder 'SMT_Coin_1FR.png']);
[image_0FR, ~, alpha_0FR]         = imread([pics_folder 'SMT_Coin_0FR.png']);

%% Make the image into a texture
image_20Cent(:,:,4) = alpha_20Cent;
image_50Cent(:,:,4) = alpha_50Cent;
image_1FR(:,:,4)    = alpha_1FR;
image_0FR(:,:,4)    = alpha_0FR;
stim.imageTexture_20Cent    = Screen('MakeTexture', window, image_20Cent);
stim.imageTexture_50Cent    = Screen('MakeTexture', window, image_50Cent);
stim.imageTexture_1FR       = Screen('MakeTexture', window, image_1FR);
stim.imageTexture_0FR       = Screen('MakeTexture', window, image_0FR);

%% fixation cross coordinates on the screen (code relative to screen Y size)
cross_length    = screenYpixels/5;
cross_thickness = 0.2*cross_length;
stim.cross_verticalLine = [xScreenCenter - (cross_thickness/2),...
    yScreenCenter - (cross_length/2),...
    xScreenCenter + (cross_thickness/2),...
    yScreenCenter + (cross_length/2)];
stim.cross_horizontalLine = [xScreenCenter - (cross_length/2),...
    yScreenCenter - (cross_thickness/2),...
    xScreenCenter + (cross_length/2),...
    yScreenCenter + (cross_thickness/2)];

end % function