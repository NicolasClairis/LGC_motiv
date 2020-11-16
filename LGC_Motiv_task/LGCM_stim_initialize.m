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

%% Enable alpha blending for anti-aliasing of all our textures
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%%
stim.VSsignal = 0;
% The color used to represent the signal
stim.signalColor = [0 0.75 0.3];

%% Money variables
moneySize = 214;
stim.money = [0 0 moneySize moneySize];

%% thresholds for reaching target
% 30% of MVC defined to this size so that 70 can be at 500 pixels
threshold_1    = [0 0 moneySize moneySize];
% 70% of MVC
threshold_2 = [0 0 500 500];

%% missed position
%
missTarget     = [0 0 moneySize moneySize];
% 60% of MVC
missThreshold = [0 0 429 429];

% First appearance of the stimulus to wait and let the participant see the value of the coin
stim.first_appearence = false;

% Is it going to be reward or punishment first
isPositiveStimuli = [randperm(numel([1 0]))-1 2];
stim.isPositiveStimuli = isPositiveStimuli(1);

%% Position of each ring/coin/signal on the screen
% force thresholds
stim.centeredthreshold_1 = CenterRectOnPointd(threshold_1, xScreenCenter, yScreenCenter);
stim.centeredthreshold_2 = CenterRectOnPointd(threshold_2, xScreenCenter, yScreenCenter);
% miss position on the screen
missThreshold = CenterRectOnPointd(missThreshold, xScreenCenter, yScreenCenter);
missTarget = CenterRectOnPointd(missTarget, xScreenCenter*2 - moneySize/2, yScreenCenter);

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

%% Prepare incentive texture
stim.incentive = [1,...
    0.5,...
    0.2,...
    0];
% be careful to respect the same order for both variables
stim.imageTextures = [stim.imageTexture_1FR,...
    stim.imageTexture_50Cent,...
    stim.imageTexture_20Cent,...
    stim.imageTexture_0FR];

%% store all info in output structure
stim.moneySize      = moneySize;
stim.missTarget     = missTarget;
stim.missThreshold  = missThreshold;
stim.threshold_1    = threshold_1;
stim.threshold_2    = threshold_2;

end % function