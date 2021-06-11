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
stim.difficulty.middle_center   = CenterRectOnPointd(difficultyRectXYsize, xScreenCenter, y_coordinates(upperBorder, visibleYsize, 1/2, difficultyRectXYsize));
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
stim.cross.colour = white;

%% prepare all instructions
% task start
stim.expWillStart.text = 'L''experimentateur va bientot demarrer la tache.';
[~,~,textSizeExpWillStart] = DrawFormattedText(window,...
    stim.expWillStart.text,...
    'center', 'center', white, wrapat);
stim.expWillStart.x = x_centerCoordinates(xScreenCenter, textSizeExpWillStart);
stim.expWillStart.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeExpWillStart);

% titles for each period
titleTextSize = scr.textSize.taskPeriodsTitles;
baselineTextSize = scr.textSize.baseline;
Screen('TextSize', window, titleTextSize);
% learning physical
stim.Ep.learning.title.text = 'Apprentissage tache physique';
[~,~,textSizeEpLearningTitle] = DrawFormattedText(window, stim.Ep.learning.title.text,...
            'center','center',white);
stim.Ep.learning.title.x = x_centerCoordinates(xScreenCenter, textSizeEpLearningTitle);
stim.Ep.learning.title.y = y_coordinates(upperBorder, visibleYsize, 1/4, textSizeEpLearningTitle);
stim.Ep.learning.title.colour = white;
% learning mental
stim.Em.learning.title.text = 'Apprentissage tache mentale';
[~,~,textSizeEmLearningTitle] = DrawFormattedText(window, stim.Em.learning.title.text,...
            'center','center',white);
stim.Em.learning.title.x = x_centerCoordinates(xScreenCenter, textSizeEmLearningTitle);
stim.Em.learning.title.y = y_coordinates(upperBorder, visibleYsize, 1/4, textSizeEmLearningTitle);
stim.Em.learning.title.colour = white;
% training physical
stim.Ep.training.title.text = 'Entrainement tache physique';
[~,~,textSizeEpTrainingTitle] = DrawFormattedText(window, stim.Ep.training.title.text,...
            'center','center',white);
stim.Ep.training.title.x = x_centerCoordinates(xScreenCenter, textSizeEpTrainingTitle);
stim.Ep.training.title.y = y_coordinates(upperBorder, visibleYsize, 1/4, textSizeEpTrainingTitle);
stim.Ep.training.title.colour = white;
% training mental
stim.Em.training.title.text = 'Entrainement tache mentale';
[~,~,textSizeEmTrainingTitle] = DrawFormattedText(window, stim.Em.training.title.text,...
            'center','center',white);
stim.Em.training.title.x = x_centerCoordinates(xScreenCenter, textSizeEmTrainingTitle);
stim.Em.training.title.y = y_coordinates(upperBorder, visibleYsize, 1/4, textSizeEmTrainingTitle);
stim.Em.training.title.colour = white;
% task physical
stim.Ep.task.title.text = 'Tache physique';
[~,~,textSizeEpTaskTitle] = DrawFormattedText(window, stim.Ep.task.title.text,...
            'center','center',white);
stim.Ep.task.title.x = x_centerCoordinates(xScreenCenter, textSizeEpTaskTitle);
stim.Ep.task.title.y = y_coordinates(upperBorder, visibleYsize, 1/4, textSizeEpTaskTitle);
stim.Ep.task.title.colour = white;
% task mental
stim.Em.task.title.text = 'Tache mentale';
[~,~,textSizeEmTaskTitle] = DrawFormattedText(window, stim.Em.task.title.text,...
            'center','center',white);
stim.Em.task.title.x = x_centerCoordinates(xScreenCenter, textSizeEmTaskTitle);
stim.Em.task.title.y = y_coordinates(upperBorder, visibleYsize, 1/4, textSizeEmTaskTitle);
stim.Em.task.title.colour = white;

% set back baseline text size
Screen('TextSize', window, baselineTextSize);

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

% end of physical training
stim.training.Ep.endMsg.text = 'Bravo! Votre entrainement physique est termine.';
[~,~,textSizeEpEndTraining] = DrawFormattedText(window,stim.training.Ep.endMsg.text,...
        'center','center',white, wrapat);
stim.training.Ep.endMsg.x = x_centerCoordinates(xScreenCenter, textSizeEpEndTraining);
stim.training.Ep.endMsg.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEpEndTraining);
stim.training.Ep.endMsg.colour = white;

% mental learning instructions
% full help instructions
stim.Em.learning.fullInstructions.text = ['Vous allez pouvoir vous entrainer a effectuer la tache. ',...
    'La couleur du chiffre represente la question posee. ',...
    'Vous pouvez repondre avec les boutons. ',...
    'La correspondance entre la position des boutons et la ',...
    'reponse que vous souhaitez donner est la suivante:'];
[~,~,textSizeEmLearningFullInstructions] = DrawFormattedText(window,...
    stim.Em.learning.fullInstructions.text,...
    'center', 'center', white, wrapat);
stim.Em.learning.fullInstructions.x = x_centerCoordinates(xScreenCenter, textSizeEmLearningFullInstructions);
stim.Em.learning.fullInstructions.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearningFullInstructions);
stim.Em.learning.fullInstructions.colour = white;
% partial help instructions
stim.Em.learning.partialInstructions.text = ['Vous allez pouvoir vous entrainer a effectuer la tache mais vous ',...
    'devrez vous rappeler de la correspondance entre les couleurs et la ',...
    'question a laquelle il faut repondre. Pour rappel:'];
[~,~,textSizeEmLearningPartialInstructions] = DrawFormattedText(window,...
    stim.Em.learning.partialInstructions.text,...
    'center', 'center', white, wrapat);
stim.Em.learning.partialInstructions.x = x_centerCoordinates(xScreenCenter, textSizeEmLearningPartialInstructions);
stim.Em.learning.partialInstructions.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearningPartialInstructions);
stim.Em.learning.partialInstructions.colour = white;
% no help instructions
stim.Em.learning.noInstructions.text = ['Vous allez pouvoir vous entrainer a effectuer la tache mais vous ',...
    'devrez vous rappeler de la correspondance entre les couleurs et la ',...
    'question a laquelle il faut repondre. Pour rappel:'];
[~,~,textSizeEmLearningNoInstructions] = DrawFormattedText(window,...
    stim.Em.learning.noInstructions.text,...
    'center', 'center', white, wrapat);
stim.Em.learning.noInstructions.x = x_centerCoordinates(xScreenCenter, textSizeEmLearningNoInstructions);
stim.Em.learning.noInstructions.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearningNoInstructions);
stim.Em.learning.noInstructions.colour = white;
% extended learning instructions
stim.Em.learning.extendedLearning.text = ['Attention, desormais, vous allez devoir repondre au chiffre ',...
    'qui vient d''etre affiche, pas au chiffre qui est a l''ecran. ',...
    ' Appuyez sur n''importe quel bouton pour le premier chiffre. ',...
    'Vous ne pourrez jamais repondre au dernier chiffre affiche. ',...
    'Pour rappel:'];
[~,~,textSizeEmLearningextendedLearning] = DrawFormattedText(window,...
    stim.Em.learning.extendedLearning.text,...
    'center', 'center', white, wrapat);
stim.Em.learning.extendedLearning.x = x_centerCoordinates(xScreenCenter, textSizeEmLearningextendedLearning);
stim.Em.learning.extendedLearning.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearningextendedLearning);
stim.Em.learning.extendedLearning.colour = white;
            
% end of mental learning
stim.training.Em.endMsg.text = 'Bravo! Votre entrainement a la tache d''effort mental est termine.';
[~,~,textSizeEmEndTraining] = DrawFormattedText(window,stim.training.Em.endMsg.text,...
        'center','center',white, wrapat);
stim.training.Em.endMsg.x = x_centerCoordinates(xScreenCenter, textSizeEmEndTraining);
stim.training.Em.endMsg.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEmEndTraining);
stim.training.Em.endMsg.colour = white;

% mental learning end of learning trial
stim.training.Em.endTrialMsg.text = 'Bravo!';
[~,~,textSizeEmTrialEndTraining] = DrawFormattedText(window,stim.training.Em.endTrialMsg.text,...
        'center','center',white, wrapat);
stim.training.Em.endTrialMsg.x = x_centerCoordinates(xScreenCenter, textSizeEmTrialEndTraining);
stim.training.Em.endTrialMsg.y = y_coordinates(upperBorder, visibleYsize, 1/4, textSizeEmTrialEndTraining);
stim.training.Em.endTrialMsg.colour = white;
stim.training.Em.endTrialMsg_bis.text = 'Au suivant!';
[~,~,textSizeEmTrialEndTraining_bis] = DrawFormattedText(window,stim.training.Em.endTrialMsg_bis.text,...
        'center','center',white, wrapat);
stim.training.Em.endTrialMsg_bis.x = x_centerCoordinates(xScreenCenter, textSizeEmTrialEndTraining_bis);
stim.training.Em.endTrialMsg_bis.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEmTrialEndTraining_bis);
stim.training.Em.endTrialMsg_bis.colour = white;

% message to press when ready
stim.pressWhenReady.text = 'Appuyez quand vous etes pret(e) a commencer.';
[~,~,textSizePressWhenReady] = DrawFormattedText(window, stim.pressWhenReady.text, 'center', 'center', white);
stim.pressWhenReady.x = x_centerCoordinates(xScreenCenter, textSizePressWhenReady);
stim.pressWhenReady.y = y_coordinates(upperBorder, visibleYsize, 15/16, textSizePressWhenReady);
stim.pressWhenReady.colour = white;

% total gains end of session
[~,~,textSizeEndMsg] = DrawFormattedText(window,['Felicitations! Cette session est maintenant terminee.',...
            'Vous avez obtenu: 0.00 chf au cours de cette session.'],'center','center',white, wrapat);
stim.endSessionMessage.x = x_centerCoordinates(xScreenCenter, textSizeEndMsg);
stim.endSessionMessage.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEndMsg);

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

% post-effort rest
stim.MVC_rest.text = 'Reposez-vous quelques secondes.';
[~,~,textSizeRest] = DrawFormattedText(window, stim.MVC_rest.text, 'center', 'center', white);
stim.MVC_rest.x = x_centerCoordinates(xScreenCenter, textSizeRest);
stim.MVC_rest.y = y_coordinates(upperBorder, visibleYsize, 4/5, textSizeRest);
stim.MVC_rest.colour = white;

% post-main task MVC calibration instructions
stim.postTaskMVCmeasurement.text = ['Pour finir cette session, nous allons vous demander ',...
        'd''essayer a nouveau de battre votre record.'];
[~,~,textSizePostTaskMVC] = DrawFormattedText(window, stim.postTaskMVCmeasurement.text,...
        'center', 'center', white, wrapat);
stim.postTaskMVCmeasurement.x = x_centerCoordinates(xScreenCenter, textSizePostTaskMVC);
stim.postTaskMVCmeasurement.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizePostTaskMVC);
stim.postTaskMVCmeasurement.colour = white;

%% mental calibration
stim.mentalCalibInstructions.text = ['Desormais vous devrez repondre dans un temps limite. Essayez de completer',...
    ' le cercle en repondant aussi vite que possible et correctement aux questions posees.'];
[~,~,textSizeMentalCalibInstructions] = DrawFormattedText(window, stim.mentalCalibInstructions.text,...
    'center', 'center', white, wrapat);
stim.mentalCalibInstructions.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibInstructions);
stim.mentalCalibInstructions.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeMentalCalibInstructions);
stim.mentalCalibInstructions.colour = white;

% calibration feedback
% success
[~,~,textSizeMentalCalibSuccess] = DrawFormattedText(window, ['Bravo vous avez tout resolu dans le temps imparti!',...
    ' Votre meilleur temps est de 0.0000 s.'],'center', 'center', white, wrapat);
stim.mentalCalibSuccessFbk.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibSuccess);
stim.mentalCalibSuccessFbk.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibSuccess);
stim.mentalCalibSuccessFbk.colour = white;
% success bis
stim.mentalCalibSuccessFbk_bis.text = 'Bravo vous avez tout resolu dans le temps imparti!';
[~,~,textSizeMentalCalibSuccess_bis] = DrawFormattedText(window, stim.mentalCalibSuccessFbk_bis.text,...
    'center', 'center', white, wrapat);
stim.mentalCalibSuccessFbk_bis.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibSuccess_bis);
stim.mentalCalibSuccessFbk_bis.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibSuccess_bis);
stim.mentalCalibSuccessFbk_bis.colour = white;
% failure
stim.mentalCalibFailureFbk.text = 'Nous allons refaire cette etape, essayez de faire mieux!';
[~,~,textSizeMentalCalibFail] = DrawFormattedText(window, stim.mentalCalibFailureFbk.text,...
    'center', 'center', white, wrapat);
stim.mentalCalibFailureFbk.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibFail);
stim.mentalCalibFailureFbk.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibFail);
stim.mentalCalibFailureFbk.colour = white;

% end of calibration
[~,~,textSizeMentalCalibEnd] = DrawFormattedText(window, ['Bravo!',...
    ' Votre meilleur temps est de 0.0000 s.'],'center', 'center', white, wrapat);
stim.mentalCalibEnd.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibEnd);
stim.mentalCalibEnd.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibEnd);
stim.mentalCalibEnd.colour = white;

%% color used to represent the effort signal
% no use of monetary images anymore
stim.difficulty.maxColor        = black;
stim.difficulty.currLevelColor  = difficultyArcColor;
stim.difficulty.ovalWidth       = 3;

%% choice period
stim.choice.choiceQuestion.text = 'QUE PREFEREZ-VOUS?';
stim.choice.choiceOR.text = 'OU';
[~,~,textSizeChoiceQuestion] = DrawFormattedText(window, stim.choice.choiceQuestion.text, 'center','center',white);
stim.choice.choiceQuestion.x = x_centerCoordinates(xScreenCenter, textSizeChoiceQuestion);
stim.choice.choiceQuestion.y = y_coordinates(upperBorder, visibleYsize, 1/8, textSizeChoiceQuestion);
stim.choice.choiceQuestion.colour = white;
[~,~,textSizeOR] = DrawFormattedText(window, stim.choice.choiceOR.text, 'center','center',white);
stim.choice.choiceOR.x = x_centerCoordinates(xScreenCenter, textSizeOR);
stim.choice.choiceOR.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeOR);
stim.choice.choiceOR.colour = white;

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
stim.loseRewardText.top_left        = [leftBorder + visibleXsize/4 - xSizeLose/2,           stim.reward.text.top_left_start(2) - ySizeLose*2.5];
stim.loseRewardText.top_right       = [leftBorder + visibleXsize*(3/4) - xSizeLose/2,       stim.reward.text.top_right_start(2) - ySizeLose*2.5];
stim.effort_introText.bottom_left   = [leftBorder + visibleXsize/4 - xSizeForEffort/2,      stim.difficulty.below_left(2) - ySizeForEffort];
stim.effort_introText.bottom_right  = [leftBorder + visibleXsize*(3/4) - xSizeForEffort/2,  stim.difficulty.below_right(2)  - ySizeForEffort];

%% release buttons message
stim.releaseButtonsMsg.text = 'Relachez les boutons svp';
[~,~,textSizeReleaseButtons] = DrawFormattedText(scr.window, stim.releaseButtonsMsg.text,'center','center',white);
stim.releaseButtonsMsg.x = x_centerCoordinates(xScreenCenter, textSizeReleaseButtons);
stim.releaseButtonsMsg.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeReleaseButtons);
stim.releaseButtonsMsg.colour = white;

%% display of the chosen option
stim.chosenOptionMsg.text = 'Vous avez choisi';
[~,~,textSizeChosenMsg] = DrawFormattedText(window, stim.chosenOptionMsg.text,'center','center',white);
stim.chosenOptionMsg.x = x_centerCoordinates(xScreenCenter, textSizeChosenMsg);
stim.chosenOptionMsg.y = y_coordinates(upperBorder, visibleYsize, 3/16, textSizeChosenMsg);
% place reward amount and difficulty level accordingly
ySizeChosenMsg = textSizeChosenMsg(4) - textSizeChosenMsg(2);
stim.chosenOption.reward = stim.reward.text.top_center_start;
stim.chosenOption.difficulty = stim.difficulty.below_center;
stim.chosenOption.squareColour = black;
stim.chosenOption.squareRect = [leftBorder + visibleXsize*(1/3),...
    upperBorder + visibleYsize*(3/16) + ySizeChosenMsg,...
    leftBorder + visibleXsize*(2/3),...
    upperBorder + visibleYsize*(11/12)];
stim.chosenOption.squareWidth = 10;
stim.winRewardText.top_center       = [xScreenCenter - xSizeWin/2,                          stim.reward.text.top_center_start(2) - ySizeWin*2.5];
stim.loseRewardText.top_center      = [xScreenCenter - xSizeLose/2,                         stim.reward.text.top_center_start(2) - ySizeLose*2.5];
stim.effort_introText.bottom_center = [xScreenCenter - xSizeForEffort/2,                    stim.difficulty.below_center(2)  - ySizeForEffort];

%% mental effort performance
% display of the relevant instructions
stim.Em.oddORevenQuestion.text = 'Chiffre pair ou impair?';
[~,~,textSizeOddORevenQuestion] = DrawFormattedText(window, stim.Em.oddORevenQuestion.text,...
    'center', 'center', white);
stim.Em.oddORevenQuestion.x = x_centerCoordinates(xScreenCenter, textSizeOddORevenQuestion);
stim.Em.oddORevenQuestion.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeOddORevenQuestion);
stim.Em.oddORevenQuestionInstructions.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeOddORevenQuestion);
% EVEN
stim.Em.even.text = 'pair';
[~,~,textSizeEven] = DrawFormattedText(window, stim.Em.even.text, 'center', 'center', white );
stim.Em.even_left.x = leftBorder + visibleXsize*(1/4) - (textSizeEven(3) - textSizeEven(1))/2;
stim.Em.even_right.x = leftBorder + visibleXsize*(3/4) - (textSizeEven(3) - textSizeEven(1))/2;
stim.Em.even.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeEven);
stim.Em.evenInstructions.y = y_coordinates(upperBorder, visibleYsize, 5/8, textSizeEven);
% OR
stim.Em.OR.text = 'OU';
[~,~,textSizeOR] = DrawFormattedText(window, stim.Em.OR.text, 'center', 'center', white );
stim.Em.OR.x = x_centerCoordinates(xScreenCenter, textSizeOR);
stim.Em.OR.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeOR);
% ODD
stim.Em.odd.text = 'impair';
[~,~,textSizeOdd] = DrawFormattedText(window, stim.Em.odd.text, 'center', 'center', white );
stim.Em.odd_left.x = leftBorder + visibleXsize*(1/4) - (textSizeOdd(3) - textSizeOdd(1))/2;
stim.Em.odd_right.x = leftBorder + visibleXsize*(3/4) - (textSizeOdd(3) - textSizeOdd(1))/2;
stim.Em.odd.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeOdd);
stim.Em.oddInstructions.y = y_coordinates(upperBorder, visibleYsize, 5/8, textSizeOdd);
% question < 5 or > 5
stim.Em.lowerORhigherQuestion.text = 'Chiffre < ou > 5?';
[~,~,textSizeLowerHigherQuestion] = DrawFormattedText(window, stim.Em.lowerORhigherQuestion.text,...
    'center', 'center', white);
stim.Em.lowerORhigherQuestion.x = x_centerCoordinates(xScreenCenter, textSizeLowerHigherQuestion);
stim.Em.lowerORhigherQuestion.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeLowerHigherQuestion);
stim.Em.lowerORhigherQuestionInstructions.y = y_coordinates(upperBorder, visibleYsize, 3/4, textSizeLowerHigherQuestion);
% < 5
stim.Em.lower.text = '< 5';
[~,~,textSizeLower] = DrawFormattedText(window, stim.Em.lower.text, 'center', 'center', white );
stim.Em.lower_left.x = leftBorder + visibleXsize*(1/4) - (textSizeLower(3) - textSizeLower(1))/2;
stim.Em.lower_right.x = leftBorder + visibleXsize*(3/4) - (textSizeLower(3) - textSizeLower(1))/2;
stim.Em.lower.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeLower);
stim.Em.lowerInstructions.y = y_coordinates(upperBorder, visibleYsize, 7/8, textSizeLower);
% > 5
stim.Em.higher.text = '< 5';
[~,~,textSizeHigher] = DrawFormattedText(window,'> 5', 'center', 'center', white );
stim.Em.higher_left.x = leftBorder + visibleXsize*(1/4) - (textSizeHigher(3) - textSizeHigher(1))/2;
stim.Em.higher_right.x = leftBorder + visibleXsize*(3/4) - (textSizeHigher(3) - textSizeHigher(1))/2;
stim.Em.higher.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeHigher);
stim.Em.higherInstructions.y = y_coordinates(upperBorder, visibleYsize, 7/8, textSizeHigher);
% press any button
stim.Em.pressAnyButtonQuestion.text = 'Appuyer sur n''importe quel bouton';
[~,~,textSizePressAnyButtonQuestion] = DrawFormattedText(window, stim.Em.pressAnyButtonQuestion.text,...
    'center', 'center', white);
stim.Em.pressAnyButtonQuestion.x = x_centerCoordinates(xScreenCenter, textSizePressAnyButtonQuestion);
stim.Em.pressAnyButtonQuestion.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizePressAnyButtonQuestion);
% press
stim.Em.pressAnyButton.text = 'Appuyer';
[~,~,textSizePressAnyButton] = DrawFormattedText(window, stim.Em.pressAnyButton.text,...
    'center', 'center', white);
stim.Em.pressAnyButton_left.x = leftBorder + visibleXsize*(1/4) - (textSizePressAnyButton(3) - textSizePressAnyButton(1))/2;
stim.Em.pressAnyButton_right.x = leftBorder + visibleXsize*(3/4) - (textSizePressAnyButton(3) - textSizePressAnyButton(1))/2;
stim.Em.pressAnyButtonQuestion.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizePressAnyButton);

% display of the number to solve
% in case you want to adjust center for each number individually
Screen('TextSize', window, scr.textSize.mentalNumber);
for iNber = [1:4, 6:9]
    n_str = num2str(iNber);
    [~,~,textSizeEmNumber] = DrawFormattedText(window, n_str,...
        'center', 'center', white);
    stim.Em.(['numberPerf_',n_str]).x = x_centerCoordinates(xScreenCenter, textSizeEmNumber);
    stim.Em.(['numberPerf_',n_str]).y = y_coordinates(upperBorder, visibleYsize, 3/4, textSizeEmNumber);
    stim.Em.(['numberPerfLearning_',n_str]).y = y_coordinates(upperBorder, visibleYsize, 1/12, textSizeEmNumber);
end
Screen('TextSize', window, scr.textSize.baseline);

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

y = upperBorder + visibleYsize*YpercentageLocation - (textSize(4) - textSize(2))/2;

% check if y is off-screen
if (y < upperBorder) || (y > (upperBorder + visibleYsize))
    error('wtf with these coordinates?');
end

end