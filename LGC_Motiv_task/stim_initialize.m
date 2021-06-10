function[stim] = stim_initialize(scr, n_E_levels)
%[stim] = stim_initialize(scr, n_E_levels)
%stim_initialize will initialize the v
%
% INPUTS
% scr: structure with main screen informations (size, center, window, etc.
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
leftBorder      = scr.leftBorder;
upperBorder     = scr.upperBorder;
visibleYsize = scr.visibleYsize;
visibleXsize = scr.visibleXsize;
wrapat = scr.wrapat;

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
stim.reward.textSizeForPTB = scr.textSize.reward;
Screen('TextSize', window, stim.reward.textSizeForPTB);
% extract reward amount text size (for choice)
[~,~,textSizeR] = DrawFormattedText(window,'+0.00 CHF', 'center', 'center', white);
xSizeText = textSizeR(3) - textSizeR(1);
ySizeText = textSizeR(4) - textSizeR(2);
stim.reward.xSizeText = xSizeText;
stim.reward.ySizeText = ySizeText;
% define where the text will be displayed
stim.reward.text.top_left_start     = [leftBorder + visibleXsize*(1/4) - xSizeText/2,   y_coordinates(upperBorder, visibleYsize, 2/5, textSizeR)]; % left option choice period
stim.reward.text.top_right_start    = [leftBorder + visibleXsize*(3/4) - xSizeText/2,   y_coordinates(upperBorder, visibleYsize, 2/5, textSizeR)]; % right option choice period
stim.reward.text.top_center_start   = [x_centerCoordinates(xScreenCenter, textSizeR),   y_coordinates(upperBorder, visibleYsize, 2/5, textSizeR)]; % chosen option display
% display on middle of the screen for performance feedback
stim.reward.text.middle_center_start = [x_centerCoordinates(xScreenCenter, textSizeR),  y_coordinates(upperBorder, visibleYsize, 1/2, textSizeR)]; % feedback

% define the colour to use for the text according to the condition
% (reward/punishment)
stim.reward.text.colour = white;
stim.punishment.text.colour = [239 138 98];

% set text back to baseline size
Screen('TextSize', window, scr.textSize.baseline);

%% difficulty rings
difficultyRectlinearSize = visibleYsize/4;
difficultyRectXYsize  = [0 0 difficultyRectlinearSize difficultyRectlinearSize];
stim.difficulty.rectSize  = difficultyRectXYsize;
% position each ring on the screen (for choice task)
stim.difficulty.below_center	= CenterRectOnPointd(difficultyRectXYsize, xScreenCenter,                   upperBorder + visibleYsize*(3/4));
stim.difficulty.below_left      = CenterRectOnPointd(difficultyRectXYsize, leftBorder + visibleXsize/4,     upperBorder + visibleYsize*(3/4));
stim.difficulty.below_right     = CenterRectOnPointd(difficultyRectXYsize, leftBorder + visibleXsize*(3/4), upperBorder + visibleYsize*(3/4));
% position each ring on the screen (for performance task)
stim.difficulty.middle_center   = CenterRectOnPointd(difficultyRectXYsize, xScreenCenter, yScreenCenter);
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
% win option
[~,~,textSizeWin] = DrawFormattedText(window,'Gagner','center','center',white);
xSizeWin = textSizeWin(3) - textSizeWin(1);
ySizeWin = textSizeWin(4) - textSizeWin(2);
stim.textRectSize.xSizeWin = xSizeWin;
stim.textRectSize.ySizeWin = ySizeWin;
% lose option
[~,~,textSizeLose] = DrawFormattedText(window,'Perdre','center','center',white);
xSizeLose = textSizeLose(3) - textSizeLose(1);
ySizeLose = textSizeLose(4) - textSizeLose(2);
stim.textRectSize.xSizeLose = xSizeLose;
stim.textRectSize.ySizeLose = ySizeLose;
% effort
[~,~,textSizeForEffort] = DrawFormattedText(window,'pour','center','center',white);
xSizeForEffort = textSizeForEffort(3) - textSizeForEffort(1);
ySizeForEffort = textSizeForEffort(4) - textSizeForEffort(2);
stim.textRectSize.xSizeForEffort = xSizeForEffort;
stim.textRectSize.ySizeForEffort = ySizeForEffort;
% extract x/y coordinates for the display of the corresponding text
stim.winRewardText.top_left         = [leftBorder + visibleXsize/4 - xSizeWin/2,            stim.reward.text.top_left_start(2) - ySizeWin*2.5];
stim.winRewardText.top_right        = [leftBorder + visibleXsize*(3/4) - xSizeWin/2,        stim.reward.text.top_right_start(2) - ySizeWin*2.5];
stim.winRewardText.top_center       = [xScreenCenter - xSizeWin/2,                          stim.reward.text.top_center_start(2) - ySizeWin*2.5];
stim.loseRewardText.top_left        = [leftBorder + visibleXsize/4 - xSizeLose/2,           stim.reward.text.top_left_start(2) - ySizeLose*2.5];
stim.loseRewardText.top_right       = [leftBorder + visibleXsize*(3/4) - xSizeLose/2,       stim.reward.text.top_right_start(2) - ySizeLose*2.5];
stim.loseRewardText.top_center      = [xScreenCenter - xSizeLose/2,                         stim.reward.text.top_center_start(2) - ySizeLose*2.5];
stim.effort_introText.bottom_left   = [leftBorder + visibleXsize/4 - xSizeForEffort/2,      stim.difficulty.below_left(2) - ySizeForEffort];
stim.effort_introText.bottom_right  = [leftBorder + visibleXsize*(3/4) - xSizeForEffort/2,  stim.difficulty.below_right(2)  - ySizeForEffort];
stim.effort_introText.bottom_center = [xScreenCenter - xSizeForEffort/2,                    stim.difficulty.below_center(2)  - ySizeForEffort];


%% color used to represent the signal
% no use of monetary images anymore
stim.difficulty.maxColor        = black;
stim.difficulty.currLevelColor  = difficultyArcColor;
stim.difficulty.ovalWidth       = 3;

%% square to display for the chosen option
[~,~,textSizeChosenMsg] = DrawFormattedText(window,'Vous avez choisi','center','center',white);
ySizeChosenMsg = textSizeChosenMsg(4) - textSizeChosenMsg(2);
stim.chosenOption.message_yCoord = y_coordinates(upperBorder, visibleYsize, 3/16, textSizeChosenMsg);
stim.chosenOption.reward = stim.reward.text.top_center_start;
stim.chosenOption.difficulty = stim.difficulty.below_center;
stim.chosenOption.squareColour = black;
stim.chosenOption.squareRect = [leftBorder + visibleXsize*(1/3),...
    upperBorder + visibleYsize*(3/16) + ySizeChosenMsg,...
    leftBorder + visibleXsize*(2/3),...
    upperBorder + visibleYsize*(11/12)];
stim.chosenOption.squareWidth = 10;

%% define bar size for the waiting time
stim.barTimeWaitRect = [leftBorder + visibleXsize*(1/4),...
    upperBorder + visibleYsize*(3/8),...
    leftBorder + visibleXsize*(3/4),...
    upperBorder + visibleYsize*(1/2)];
stim.barTimeWait.colour = white;

% accompanying text
remainingTimeText = 'Temps restant';
[~,~,remainingTimeTextSize] = DrawFormattedText(window,remainingTimeText,'center','center',white);
stim.remainingTime.text = remainingTimeText;
stim.remainingTime.x    = x_centerCoordinates(xScreenCenter, remainingTimeTextSize);
stim.remainingTime.y    = y_coordinates(upperBorder, visibleYsize, 1/4, remainingTimeTextSize);
stim.remainingTime.colour = white;


%% fixation cross coordinates on the screen (code relative to screen Y size)
cross_length    = visibleYsize/6;
cross_thickness = 0.2*cross_length;
stim.cross.verticalLine = [xScreenCenter - (cross_thickness/2),...
    yScreenCenter - (cross_length/2),...
    xScreenCenter + (cross_thickness/2),...
    yScreenCenter + (cross_length/2)];
stim.cross.horizontalLine = [xScreenCenter - (cross_length/2),...
    yScreenCenter - (cross_thickness/2),...
    xScreenCenter + (cross_length/2),...
    yScreenCenter + (cross_thickness/2)];

%% prepare all instructions
% task start
stim.expWillStart.text = 'L''experimentateur va bientot demarrer la tache.';
[~,~,textSizeExpWillStart] = DrawFormattedText(window,...
    stim.expWillStart.text,...
    'center', 'center', white, wrapat);
stim.expWillStart.x = x_centerCoordinates(xScreenCenter, textSizeExpWillStart);
stim.expWillStart.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeExpWillStart);

% reward training instructions
stim.training.R.text = ['Vous allez a present choisir entre deux options associees a differents niveaux de recompense et d''effort '...
                'l''option qui vous parait la plus interessante.'];
[~,~,textSizeRewardTraining] = DrawFormattedText(window,...
                stim.training.R.text,...
                'center', 'center', white, wrapat);
stim.training.R.x = x_centerCoordinates(xScreenCenter, textSizeRewardTraining);
stim.training.R.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeRewardTraining);
stim.training.R.colour = white;

% punishment training instructions
stim.training.P.text = ['Vous allez a present choisir entre deux options associees a differents niveaux de pertes et d''effort '...
                'l''option qui vous parait la moins penible.'];
[~,~,textSizePunishmentTraining] = DrawFormattedText(window,...
                stim.training.P.text,...
                'center', 'center', white, wrapat);
stim.training.P.x = x_centerCoordinates(xScreenCenter, textSizePunishmentTraining);
stim.training.P.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizePunishmentTraining);
stim.training.P.colour = white;

% reward + punishment training
stim.training.RP.text = ['Vous allez a present choisir entre deux options associees a differents niveaux de recompenses ou de pertes et d''effort '...
                'l''option qui vous parait preferable.'];
[~,~,textSizeRewardAndPunishmentTraining] = DrawFormattedText(window,...
                stim.training.RP.text,...
                'center', 'center', white, wrapat);
stim.training.RP.x = x_centerCoordinates(xScreenCenter, textSizeRewardAndPunishmentTraining);
stim.training.RP.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeRewardAndPunishmentTraining);
stim.training.RP.colour = white;

% message to press when ready
stim.pressWhenReady.text = 'Appuyez quand vous etes pret(e) a commencer la tache.';
[~,~,textSizePressWhenReady] = DrawFormattedText(window, stim.pressWhenReady.text, 'center', 'center', white);
stim.pressWhenReady.x = x_centerCoordinates(xScreenCenter, textSizePressWhenReady);
stim.pressWhenReady.y = y_coordinates(upperBorder, visibleYsize, 15/6, textSizePressWhenReady);
stim.pressWhenReady.colour = white;

% total gains end of session
[~,~,textSizeEndMsg] = DrawFormattedText(window,['Felicitations! Cette session est maintenant terminee.',...
            'Vous avez obtenu: 0.00 chf au cours de cette session.'],'center','center',white, wrapat);
stim.endSessionMessage.x = x_centerCoordinates(xScreenCenter, textSizeEndMsg);
stim.endSessionMessage.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeEndMsg);

%% MVC calibration for physical effort
% MVC instructions
stim.Ep.MVC.instructions.text = ['Avant de commencer l''experience, ',...
    'nous allons vous demander ',...
    'de serrer la poignee de force au maximum de vos capacites plusieurs ',...
    'fois d''affilee.'];
[~,~,textSizeMVCInstructions] = DrawFormattedText(window, stim.Ep.MVC.instructions.text, 'center','center', white, wrapat);
stim.Ep.MVC.instructions.x = x_centerCoordinates(xScreenCenter, textSizeMVCInstructions);
stim.Ep.MVC.instructions.y = y_coordinates(upperBorder, visibleYsize, 7/10, textSizeMVCInstructions);
stim.Ep.MVC.instructions.colour = white;
stim.Ep.MVC.instructions_bis.text = 'Tenez-vous pret a serrer la poignee.';
[~,~,textSizeMVCInstructions_bis] = DrawFormattedText(window, stim.Ep.MVC.instructions_bis.text, 'center', 'center', white);
stim.Ep.MVC.instructions_bis.x = x_centerCoordinates(xScreenCenter, textSizeMVCInstructions_bis);
stim.Ep.MVC.instructions_bis.y = y_coordinates(upperBorder, visibleYsize, 3/10, textSizeMVCInstructions_bis);
stim.Ep.MVC.instructions_bis.colour = white;

% GO instruction
stim.Ep.MVC.GO.text = 'GO !';
[~,~,textSizeGO] = DrawFormattedText(window, stim.Ep.MVC.GO.text, 'center', 'center', white);
stim.Ep.MVC.GO.x = x_centerCoordinates(xScreenCenter, textSizeGO);
stim.MVC_rest.y = y_coordinates(upperBorder, visibleYsize, 9/10, textSizeGO);
stim.Ep.MVC.GO.colour = white;

% effort scale


% post-effort rest
stim.MVC_rest.text = 'Reposez-vous quelques secondes.';
[~,~,textSizeRest] = DrawFormattedText(window, stim.MVC_rest.text, 'center', 'center', white);
stim.MVC_rest.x = x_centerCoordinates(xScreenCenter, textSizeRest);
stim.MVC_rest.y = y_coordinates(upperBorder, visibleYsize, 4/5, textSizeRest);
stim.MVC_rest.colour = white;

%% prepare feedback messages
% reward feedback
stim.feedback.reward.text = 'Vous avez obtenu';
[~,~,textSizeRewardFbkMsg] = DrawFormattedText(window, stim.feedback.reward.text,...
    'center', 'center',...
    white);
stim.feedback.reward.x = x_centerCoordinates(xScreenCenter, textSizeRewardFbkMsg);
stim.feedback.reward.y = y_coordinates(upperBorder, visibleYsize, 3/8, textSizeRewardFbkMsg);
stim.feedback.colour = white;

% punishment feedback
stim.feedback.punishment.text = 'Vous avez perdu';
[~,~,textSizePunishmentFbkMsg] = DrawFormattedText(window, stim.feedback.punishment.text,...
    'center', 'center',...
    white);
stim.feedback.punishment.x = x_centerCoordinates(xScreenCenter, textSizePunishmentFbkMsg);
stim.feedback.punishment.y = y_coordinates(upperBorder, visibleYsize, 3/8, textSizePunishmentFbkMsg);

% error: too slow feedback
stim.feedback.error_tooSlow.text = 'Trop lent!';
[~,~,textSizeErrorTooSlowFbkMsg] = DrawFormattedText(window, stim.feedback.error_tooSlow.text,...
    'center', 'center',...
    white);
stim.feedback.error_tooSlow.x = x_centerCoordinates(xScreenCenter, textSizeErrorTooSlowFbkMsg);
stim.feedback.error_tooSlow.y = y_coordinates(upperBorder, visibleYsize, 3/8, textSizeErrorTooSlowFbkMsg); % used to be 1/5*yScreenCenter

% error too many errors feedback
stim.feedback.error_tooManyErrors.text = 'Trop d''erreurs!';
[~,~,textSizeErrorTooManyErrorsFbkMsg] = DrawFormattedText(window, stim.feedback.error_tooManyErrors.text,...
    'center', 'center',...
    white);
stim.feedback.error_tooManyErrors.x = x_centerCoordinates(xScreenCenter, textSizeErrorTooManyErrorsFbkMsg);
stim.feedback.error_tooManyErrors.y = y_coordinates(upperBorder, visibleYsize, 3/8, textSizeErrorTooManyErrorsFbkMsg); % used to be 1/5*yScreenCenter


% for the end of the performance period circle to signify end of the trial (win or loss)
stim.endTrialcircle  = [0, 0, (difficultyRectlinearSize + (difficultyRectlinearSize/5)), (difficultyRectlinearSize + (difficultyRectlinearSize/5) )];
stim.end_trial.middle_center = CenterRectOnPointd(stim.endTrialcircle, xScreenCenter, yScreenCenter);


%% add grey screen on top to be sure that this does not actually appear on
% the screen
Screen('FillRect',window, grey, [0 0 xScreenCenter*2 yScreenCenter*2]);
Screen(window,'Flip');

end % function

function[x] = x_centerCoordinates(xScreenCenter, textSize)
% to center the coordinates on the screen
x = xScreenCenter - (textSize(3) - textSize(1))/2;
end

function[y] = y_coordinates(upperBorder, visibleYsize, YpercentageLocation, textSize)
% define where the stimulus will be displayed on the Y axis: start after
% the non-visible part of the upper border, then define the location on the
% visible part of the screen and center the text on this location

y = upperBorder + visibleYsize*YpercentageLocation - (textSize(4) - textSize(2)/2);

end