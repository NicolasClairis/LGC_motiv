function[stim] = stim_initialize(scr, n_E_levels, langage)
%[stim] = stim_initialize(scr, n_E_levels, langage)
%stim_initialize will initialize most of the visual stimuli used in the
%task.
%
% INPUTS
% scr: structure with main screen informations (size, center, window, etc.
%
% n_E_levels: number of difficulty levels
%
% langage:
% 'fr': display instructions in french
% 'engl': display instructions in english
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
orange = scr.colours.orange;
red = scr.colours.red;
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
xSizeTextR = textSizeR(3) - textSizeR(1);
ySizeTextR = textSizeR(4) - textSizeR(2);
stim.reward.xSizeText = xSizeTextR;
stim.reward.ySizeText = ySizeTextR;
% define where the text will be displayed during choice
stim.reward.text.top_left_start     = [leftBorder + visibleXsize*(1/4) - xSizeTextR/2,   y_coordinates(upperBorder, visibleYsize, 2/5, textSizeR)]; % left option choice period
stim.reward.text.top_right_start    = [leftBorder + visibleXsize*(3/4) - xSizeTextR/2,   y_coordinates(upperBorder, visibleYsize, 2/5, textSizeR)]; % right option choice period
% display on middle of the screen for performance feedback
stim.reward.text.middle_center_start = [x_centerCoordinates(xScreenCenter, textSizeR),  y_coordinates(upperBorder, visibleYsize, 1/2, textSizeR)]; % feedback
% colour for the text
stim.reward.text.colour = white;

% same for punishments
[~,~,textSizeP] = DrawFormattedText(window,'-0.00 CHF', 'center', 'center', white);
xSizeTextP = textSizeP(3) - textSizeP(1);
ySizeTextP = textSizeP(4) - textSizeP(2);
stim.punishment.xSizeText = xSizeTextP;
stim.punishment.ySizeText = ySizeTextP;
% define where the text will be displayed during choice
stim.punishment.text.top_left_start     = [leftBorder + visibleXsize*(1/4) - xSizeTextR/2,   y_coordinates(upperBorder, visibleYsize, 2/5, textSizeP)]; % left option choice period
stim.punishment.text.top_right_start    = [leftBorder + visibleXsize*(3/4) - xSizeTextR/2,   y_coordinates(upperBorder, visibleYsize, 2/5, textSizeP)]; % right option choice period
% display on middle of the screen for performance feedback
stim.punishment.text.middle_center_start = [x_centerCoordinates(xScreenCenter, textSizeP),  y_coordinates(upperBorder, visibleYsize, 1/2, textSizeP)]; % feedback
% stim.punishment.text.colour = [239 138 98];
stim.punishment.text.colour = white;

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
% note level 0 = easiest level (ascending order)
% first: load timings for physical effort so as to calibrate everything
% based on this
[Ep_time_levels] = physical_effortLevels(n_E_levels);
% normalize all values by maximal effort
E_maxDuration = Ep_time_levels.(['level_',num2str(n_E_levels-1)]);
for iDiff = 0:(n_E_levels-1)
    % extract name for subfield of the current difficulty level
    diff_level_nm = ['level_',num2str(iDiff)];
    E_durPerc_tmp = Ep_time_levels.(['level_',num2str(iDiff)]);
    
    % extract angle for the arc which will correspond to the difficulty
    % level: max circle = max difficulty level
    startAngle_tmp = stim.difficulty.arcEndAngle*( (E_maxDuration - E_durPerc_tmp)/E_maxDuration);
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
switch langage
    case 'fr'
        stim.expWillStart.text = 'L''experimentateur va bientot demarrer la tache.';
    case 'engl'
        stim.expWillStart.text = 'The experimenter will soon start the task.';
end
[~,~,textSizeExpWillStart] = DrawFormattedText(window,...
    stim.expWillStart.text,...
    'center', 'center', white, wrapat);
stim.expWillStart.x = x_centerCoordinates(xScreenCenter, textSizeExpWillStart);
stim.expWillStart.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeExpWillStart);

% press space
switch langage
    case 'fr'
        stim.pressSpace.text = 'Appuyez sur espace quand vous etes pret(e) a demarrer.';
    case 'engl'
        stim.pressSpace.text = 'Press space key when you are ready to start.';
end
[~,~,textSizePressSpace] = DrawFormattedText(window,...
    stim.pressSpace.text,...
    'center', 'center', white, wrapat);
stim.pressSpace.x = x_centerCoordinates(xScreenCenter, textSizePressSpace);
stim.pressSpace.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizePressSpace);
        
%% titles for each period
titleTextSize = scr.textSize.taskPeriodsTitles;
baselineTextSize = scr.textSize.baseline;
Screen('TextSize', window, titleTextSize);
% learning physical
switch langage
    case 'fr'
        stim.Ep.learning.title.text = 'Apprentissage tache d''effort physique';
    case 'engl'
        stim.Ep.learning.title.text = 'Learning physical effort task';
end
[~,~,textSizeEpLearningTitle] = DrawFormattedText(window, stim.Ep.learning.title.text,...
    'center','center',white);
stim.Ep.learning.title.x = x_centerCoordinates(xScreenCenter, textSizeEpLearningTitle);
stim.Ep.learning.title.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEpLearningTitle);
stim.Ep.learning.title.colour = white;
% learning mental
switch langage
    case 'fr'
        stim.Em.learning.title.text = 'Apprentissage tache d''effort mental';
    case 'engl'
        stim.Em.learning.title.text = 'Learning mental effort task';
end
[~,~,textSizeEmLearningTitle] = DrawFormattedText(window, stim.Em.learning.title.text,...
    'center','center',white);
stim.Em.learning.title.x = x_centerCoordinates(xScreenCenter, textSizeEmLearningTitle);
stim.Em.learning.title.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEmLearningTitle);
stim.Em.learning.title.colour = white;
% training physical
switch langage
    case 'fr'
        stim.Ep.training.title.text = 'Entrainement tache d''effort physique';
    case 'engl'
        stim.Ep.training.title.text = 'Training physical effort task';
end

[~,~,textSizeEpTrainingTitle] = DrawFormattedText(window, stim.Ep.training.title.text,...
    'center','center',white);
stim.Ep.training.title.x = x_centerCoordinates(xScreenCenter, textSizeEpTrainingTitle);
stim.Ep.training.title.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEpTrainingTitle);
stim.Ep.training.title.colour = white;
% training mental
switch langage
    case 'fr'
        stim.Em.training.title.text = 'Entrainement tache mentale';
    case 'engl'
        stim.Em.training.title.text = 'Training mental effort task';
end
[~,~,textSizeEmTrainingTitle] = DrawFormattedText(window, stim.Em.training.title.text,...
    'center','center',white);
stim.Em.training.title.x = x_centerCoordinates(xScreenCenter, textSizeEmTrainingTitle);
stim.Em.training.title.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEmTrainingTitle);
stim.Em.training.title.colour = white;
% task physical
switch langage
    case 'fr'
        stim.Ep.task.title.text = 'Tache d''effort physique';
    case 'engl'
        stim.Ep.task.title.text = 'Physical effort task';
end

[~,~,textSizeEpTaskTitle] = DrawFormattedText(window, stim.Ep.task.title.text,...
    'center','center',white);
stim.Ep.task.title.x = x_centerCoordinates(xScreenCenter, textSizeEpTaskTitle);
stim.Ep.task.title.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEpTaskTitle);
stim.Ep.task.title.colour = white;
% task mental
switch langage
    case 'fr'
        stim.Em.task.title.text = 'Tache d''effort mental';
    case 'engl'
        stim.Em.task.title.text = 'Mental effort task';
end
[~,~,textSizeEmTaskTitle] = DrawFormattedText(window, stim.Em.task.title.text,...
    'center','center',white);
stim.Em.task.title.x = x_centerCoordinates(xScreenCenter, textSizeEmTaskTitle);
stim.Em.task.title.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEmTaskTitle);
stim.Em.task.title.colour = white;

% set back baseline text size
Screen('TextSize', window, baselineTextSize);

%% instructions for each period of the training

% learning instructions
switch langage
    case 'fr'
        stim.EpEm_learning.text = ['Vous allez maintenant faire une serie d''efforts sans pression temporelle. ',...
            'Cela va vous permettre de voir quel niveau d''effort correspond a chaque taille du cercle dans la tache principale.'];
    case 'engl'
        stim.EpEm_learning.text = ['You will now do a series of efforts without temporal pressure. ',...
            'This will allow you to understand what level of effort corresponds to each size of the circle in the main task.'];
end
[~,~,textSizeLearning] = DrawFormattedText(window,...
    stim.EpEm_learning.text,...
    'center', 'center', white, wrapat);
stim.EpEm_learning.x = x_centerCoordinates(xScreenCenter, textSizeLearning);
stim.EpEm_learning.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeLearning);
stim.EpEm_learning.colour = white;

% reward training instructions
switch langage
    case 'fr'
        stim.training.R.text = ['Vous allez a present choisir entre deux options associees a differents niveaux d''argent et d''effort '...
            'l''option qui vous parait la plus interessante.'];
    case 'engl'
        stim.training.R.text = ['You will now choose between two options associated with different levels of money and effort ',...
            'the option which seems the most interesting for you.'];
end
[~,~,textSizeRewardTraining] = DrawFormattedText(window,...
    stim.training.R.text,...
    'center', 'center', white, wrapat);
stim.training.R.x = x_centerCoordinates(xScreenCenter, textSizeRewardTraining);
stim.training.R.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeRewardTraining);
stim.training.R.colour = white;

% punishment training instructions
switch langage
    case 'fr'
        stim.training.P.text = ['Vous allez a present choisir entre deux options associees a differents niveaux d''argent et d''effort '...
            'l''option qui vous parait la moins penible.'];
    case 'engl'
        stim.training.P.text = ['You will now choose between two options associated with different levels of money and effort ',...
            'the option which seems the least aversive for you.'];
end
[~,~,textSizePunishmentTraining] = DrawFormattedText(window,...
    stim.training.P.text,...
    'center', 'center', white, wrapat);
stim.training.P.x = x_centerCoordinates(xScreenCenter, textSizePunishmentTraining);
stim.training.P.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizePunishmentTraining);
stim.training.P.colour = white;

% reward + punishment training (with 2 buttons, no confidence)
switch langage
    case 'fr'
        stim.training.RP.text = ['Vous allez a present choisir entre deux options associees a differents niveaux d''argent et d''effort ',...
            'l''option qui vous parait preferable.'];
    case 'engl'
        stim.training.RP.text = ['You will now choose between two options associated with different levels of money and effort ',...
            'the option which seems the best for you.'];
end
[~,~,textSizeRewardAndPunishmentTraining] = DrawFormattedText(window,...
    stim.training.RP.text,...
    'center', 'center', white, wrapat);
stim.training.RP.x = x_centerCoordinates(xScreenCenter, textSizeRewardAndPunishmentTraining);
stim.training.RP.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeRewardAndPunishmentTraining);
stim.training.RP.colour = white;

% reward + punishment training with confidence mapping
switch langage
    case 'fr'
        stim.training.RP_withConfMapping.text = ['Vous allez a present choisir entre deux options associees a differents niveaux d''argent et d''effort. ',...
            'Vous aurez un temps limite pour repondre et pour faire l''effort. ',...
            'Utilisez les quatre boutons pour exprimer votre certitude d''avoir choisi l''option la plus interessante pour vous.'];
    case 'engl'
        stim.training.RP_withConfMapping.text = ['You will now choose between two options associated with different levels of money and effort. ',...
            'You will have a limited time to answer and to do the effort.',...
            'Use the four buttons to express your degree of confidence on having selected the best option for you.'];
end
[~,~,textSizeRewardAndPunishmentConfMapTraining] = DrawFormattedText(window,...
    stim.training.RP_withConfMapping.text,...
    'center', 'center', white, wrapat);
stim.training.RP_withConfMapping.x = x_centerCoordinates(xScreenCenter, textSizeRewardAndPunishmentConfMapTraining);
stim.training.RP_withConfMapping.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeRewardAndPunishmentConfMapTraining);
stim.training.RP_withConfMapping.colour = white;

% reward + punishment training without confidence mapping but 4 buttons
switch langage
    case 'fr'
        stim.training.RP_withoutConfMapping.text = ['Vous allez maintenant refaire la meme chose, ',...
            'mais la correspondance avec les boutons ne sera plus affichee. ',...
            'Pensez bien a utiliser les quatre boutons pour repondre.'];
    case 'engl'
        stim.training.RP_withoutConfMapping.text = ['You will now do the same, ',...
            'but the mapping with the buttons will not be displayed anymore. ',...
            'Keep using the four buttons when you answer.'];
end
[~,~,textSizeRewardAndPunishmentWihoutConfMapTraining] = DrawFormattedText(window,...
    stim.training.RP.text,...
    'center', 'center', white, wrapat);
stim.training.RP_withoutConfMapping.x = x_centerCoordinates(xScreenCenter, textSizeRewardAndPunishmentWihoutConfMapTraining);
stim.training.RP_withoutConfMapping.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeRewardAndPunishmentWihoutConfMapTraining);
stim.training.RP_withoutConfMapping.colour = white;

% end of physical training
switch langage
    case 'fr'
        stim.training.Ep.endMsg.text = 'Bravo! Votre entrainement a la tache d''effort physique est termine.';
    case 'engl'
        stim.training.Ep.endMsg.text = 'Congratulations! Your training for the physical effort task is now completed.';
end
[~,~,textSizeEpEndTraining] = DrawFormattedText(window,stim.training.Ep.endMsg.text,...
    'center','center',white, wrapat);
stim.training.Ep.endMsg.x = x_centerCoordinates(xScreenCenter, textSizeEpEndTraining);
stim.training.Ep.endMsg.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEpEndTraining);
stim.training.Ep.endMsg.colour = white;

% mental learning instructions
% full help instructions
switch langage
    case 'fr'
        stim.Em.learning.fullInstructions.text = ['Indiquez si le chiffre a l''ecran est inferieur a 5 (a gauche) ',...
            'ou superieur (a droite)'];
        %'La couleur du chiffre represente la question posee. ',... only
        %for task switching version
    case 'engl'
        stim.Em.learning.fullInstructions.text = ['Indicate whether the number on screen is below 5 (to the left) ',...
            'or above 5 (to the right).'];
%         'The colour of the number represents the nature of the question.
%         ',... only for task switching version
end
[~,~,textSizeEmLearningFullInstructions] = DrawFormattedText(window,...
    stim.Em.learning.fullInstructions.text,...
    'center', 'center', white, wrapat);
stim.Em.learning.fullInstructions.x = x_centerCoordinates(xScreenCenter, textSizeEmLearningFullInstructions);
stim.Em.learning.fullInstructions.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearningFullInstructions);
stim.Em.learning.fullInstructions.colour = white;
% partial help instructions
switch langage
    case 'fr'
        stim.Em.learning.partialInstructions.text = ['Dorenavant, vous ',...
            'devrez vous rappeler de la correspondance entre chaque bouton et la ',...
            'reponse correspondante. Pour rappel:'];
    case 'engl'
stim.Em.learning.partialInstructions.text = ['Now you ',...
            'will need to remember the mapping between each button and the corresponding answer. ',...
            'Here is a quick reminder before starting:'];
end
[~,~,textSizeEmLearningPartialInstructions] = DrawFormattedText(window,...
    stim.Em.learning.partialInstructions.text,...
    'center', 'center', white, wrapat);
stim.Em.learning.partialInstructions.x = x_centerCoordinates(xScreenCenter, textSizeEmLearningPartialInstructions);
stim.Em.learning.partialInstructions.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearningPartialInstructions);
stim.Em.learning.partialInstructions.colour = white;
% no help instructions
switch langage
    case 'fr'
        stim.Em.learning.noInstructions.text = ['Dorenavant, vous ',...
            'devrez vous rappeler de la correspondance entre chaque bouton et la ',...
            'reponse correspondante. Pour rappel:'];
    case 'engl'
        stim.Em.learning.noInstructions.text = ['You can now train to perform the task but you ',...
            'will need to remember the mapping between each button and the corresponding answer. ',...
            'Here is a quick reminder before starting:'];
end
[~,~,textSizeEmLearningNoInstructions] = DrawFormattedText(window,...
    stim.Em.learning.noInstructions.text,...
    'center', 'center', white, wrapat);
stim.Em.learning.noInstructions.x = x_centerCoordinates(xScreenCenter, textSizeEmLearningNoInstructions);
stim.Em.learning.noInstructions.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearningNoInstructions);
stim.Em.learning.noInstructions.colour = white;
% extended learning instructions
switch langage
    case 'fr'
        stim.Em.learning.learning_Nback0.text = ['Vous allez a present vous entrainer ',...
            'a nouveau sur les differents niveaux de difficulte que vous rencontrerez dans la tache. ',...
            'Pour rappel:'];
        stim.Em.learning.learning_Nback1.text = ['Attention, vous allez devoir repondre au chiffre ',...
            'qui vient d''etre affiche, pas au chiffre qui est a l''ecran. ',...
            ' Appuyez sur n''importe quel bouton pour le premier chiffre. '];
        stim.Em.learning.learning_Nback2.text = ['Attention, ',...
            'vous allez devoir repondre au chiffre affiche 2 chiffres plus tot, ',...
            'pas au chiffre qui est a l''ecran. ',...
            ' Appuyez sur n''importe quel bouton pour les deux premiers chiffres. '];
        stim.Em.learning.learning_Nback2_bis.text = ['Il faut toujours repondre avec un decalage de 2. ',...
            'Il faudra aussi aller le plus vite possible. Le trait rouge indique le minimum a atteindre a chaque essai. ',...
            'Le trait orange indiquera votre meilleure performance.'];
    case 'engl'
        stim.Em.learning.learning_Nback0.text = ['You will train ',...
            'again on the different levels of difficulty that you will encounter in the task. ',...
            'Here is a quick reminder of the mapping:'];
        stim.Em.learning.learning_Nback1.text = ['Attention, ',...
            'you will have to respond to the number which has just been displayed, ',...
            'not to the number which is currently on the screen. ',...
            'Press any button for the first digit. '];
        stim.Em.learning.learning_Nback2.text = ['Attention, ',...
            'you will have to respond to the number which has been displayed 2 numbers before, ',...
            'not to the number which is currently on the screen. ',...
            'Press any button for the two first digits. '];
        stim.Em.learning.learning_Nback2_bis.text = ['You still have to answer with a delay of 2. ',...
            'You also have to be as fast as possible. The red line indicates the minimum to reach on every trial. ',...
            'The orange line indicates your best performance.'];
end
[~,~,textSizeEmLearning_Nback0] = DrawFormattedText(window,...
    stim.Em.learning.learning_Nback0.text,...
    'center', 'center', white, wrapat);
[~,~,textSizeEmLearning_Nback1] = DrawFormattedText(window,...
    stim.Em.learning.learning_Nback1.text,...
    'center', 'center', white, wrapat);
[~,~,textSizeEmLearning_Nback2] = DrawFormattedText(window,...
    stim.Em.learning.learning_Nback2.text,...
    'center', 'center', white, wrapat);
[~,~,textSizeEmLearning_Nback2_bis] = DrawFormattedText(window,...
    stim.Em.learning.learning_Nback2.text,...
    'center', 'center', white, wrapat);
stim.Em.learning.learning_Nback0.x = x_centerCoordinates(xScreenCenter, textSizeEmLearning_Nback0);
stim.Em.learning.learning_Nback0.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearning_Nback0);
stim.Em.learning.learning_Nback0.colour = white;
% Nback = 1
stim.Em.learning.learning_Nback1.x = x_centerCoordinates(xScreenCenter, textSizeEmLearning_Nback1);
stim.Em.learning.learning_Nback1.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearning_Nback1);
stim.Em.learning.learning_Nback1.colour = white;
% Nback = 2
stim.Em.learning.learning_Nback2.x = x_centerCoordinates(xScreenCenter, textSizeEmLearning_Nback2);
stim.Em.learning.learning_Nback2.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearning_Nback2);
stim.Em.learning.learning_Nback2.colour = white;
% Nback = 2 bis
stim.Em.learning.learning_Nback2_bis.x = x_centerCoordinates(xScreenCenter, textSizeEmLearning_Nback2_bis);
stim.Em.learning.learning_Nback2_bis.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeEmLearning_Nback2_bis);
stim.Em.learning.learning_Nback2_bis.colour = white;

% end of mental learning
switch langage
    case 'fr'
        stim.training.Em.endMsg.text = 'Bravo! Votre entrainement a la tache d''effort mental est termine.';
    case 'engl'
        stim.training.Em.endMsg.text = 'Congratulations! Your training for the mental effort task is now completed.';
end
[~,~,textSizeEmEndTraining] = DrawFormattedText(window,stim.training.Em.endMsg.text,...
    'center','center',white, wrapat);
stim.training.Em.endMsg.x = x_centerCoordinates(xScreenCenter, textSizeEmEndTraining);
stim.training.Em.endMsg.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEmEndTraining);
stim.training.Em.endMsg.colour = white;

% mental learning end of learning trial
switch langage
    case 'fr'
        stim.training.Em.endTrialMsg.text = 'Bravo!';
    case 'engl'
        stim.training.Em.endTrialMsg.text = 'Congratulations!';
end
[~,~,textSizeEmTrialEndTraining] = DrawFormattedText(window,stim.training.Em.endTrialMsg.text,...
    'center','center',white, wrapat);
stim.training.Em.endTrialMsg.x = x_centerCoordinates(xScreenCenter, textSizeEmTrialEndTraining);
stim.training.Em.endTrialMsg.y = y_coordinates(upperBorder, visibleYsize, 1/4, textSizeEmTrialEndTraining);
stim.training.Em.endTrialMsg.colour = white;
switch langage
    case 'fr'
        stim.training.Em.endTrialMsg_bis.text = 'Au suivant!';
    case 'engl'
        stim.training.Em.endTrialMsg_bis.text = 'Next!';
end
[~,~,textSizeEmTrialEndTraining_bis] = DrawFormattedText(window,stim.training.Em.endTrialMsg_bis.text,...
    'center','center',white, wrapat);
stim.training.Em.endTrialMsg_bis.x = x_centerCoordinates(xScreenCenter, textSizeEmTrialEndTraining_bis);
stim.training.Em.endTrialMsg_bis.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEmTrialEndTraining_bis);
stim.training.Em.endTrialMsg_bis.colour = white;

% message to press when ready
switch langage
    case 'fr'
        stim.pressWhenReady.text = 'Appuyez quand vous etes pret(e) a commencer.';
    case 'engl'
        stim.pressWhenReady.text = 'Press when you are ready to start.';
end
[~,~,textSizePressWhenReady] = DrawFormattedText(window, stim.pressWhenReady.text, 'center', 'center', white);
stim.pressWhenReady.x = x_centerCoordinates(xScreenCenter, textSizePressWhenReady);
stim.pressWhenReady.y = y_coordinates(upperBorder, visibleYsize, 15/16, textSizePressWhenReady);
stim.pressWhenReady.colour = white;

% end of session (before last calibration)
switch langage
    case 'fr'
        stim.endfMRIMessage.text = ['Nous allons maintenant vous demander ',...
            'de refaire votre maximum apres quelques secondes de pause.'];
    case 'engl'
        stim.endfMRIMessage.text = ['We will now ask you ',...
            'to perform your maximum after a few seconds of break.'];
end
[~,~,textSizeEndfMRIMsg] = DrawFormattedText(window,stim.endfMRIMessage.text,'center','center',white, wrapat);
stim.endfMRIMessage.x = x_centerCoordinates(xScreenCenter, textSizeEndfMRIMsg);
stim.endfMRIMessage.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEndfMRIMsg);

% total gains end of session
switch langage
    case 'fr'
        [~,~,textSizeEndMsg] = DrawFormattedText(window,['Felicitations! Cette session est maintenant terminee.',...
            'Vous avez obtenu: 0.00 chf au cours de cette session.'],'center','center',white, wrapat);
    case 'engl'
        [~,~,textSizeEndMsg] = DrawFormattedText(window,['Congratulations! This session is now completed.',...
            'You got: 0.00 chf during this session.'],'center','center',white, wrapat);
end
stim.endSessionMessage.x = x_centerCoordinates(xScreenCenter, textSizeEndMsg);
stim.endSessionMessage.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeEndMsg);

%% MVC calibration for physical effort
% MVC instructions
switch langage
    case 'fr'
        stim.Ep.MVC.instructions.text = ['Serrez la poignee de force au maximum de vos capacites. ',...
            'Chaque essai sera remunere.'];
    case 'engl'
        stim.Ep.MVC.instructions.text = ['Tighten the grip at your maximum. ',...
            'Each trial will be rewarded.'];
end
[~,~,textSizeMVCInstructions] = DrawFormattedText(window, stim.Ep.MVC.instructions.text, 'center','center', white, wrapat);
stim.Ep.MVC.instructions.x = x_centerCoordinates(xScreenCenter, textSizeMVCInstructions);
stim.Ep.MVC.instructions.y = y_coordinates(upperBorder, visibleYsize, 3/10, textSizeMVCInstructions);
stim.Ep.MVC.instructions.colour = white;
switch langage
    case 'fr'
        stim.Ep.MVC.instructions_bis.text = 'Tenez-vous pret a serrer la poignee.';
    case 'engl'
        stim.Ep.MVC.instructions_bis.text = 'Be ready to squeeze the grip.';
end
[~,~,textSizeMVCInstructions_bis] = DrawFormattedText(window, stim.Ep.MVC.instructions_bis.text, 'center', 'center', white);
stim.Ep.MVC.instructions_bis.x = x_centerCoordinates(xScreenCenter, textSizeMVCInstructions_bis);
stim.Ep.MVC.instructions_bis.y = y_coordinates(upperBorder, visibleYsize, 7/10, textSizeMVCInstructions_bis);
stim.Ep.MVC.instructions_bis.colour = white;

% GO instruction
stim.Ep.MVC.GO.text = 'GO !';
[~,~,textSizeGO] = DrawFormattedText(window, stim.Ep.MVC.GO.text, 'center', 'center', white);
stim.Ep.MVC.GO.x = x_centerCoordinates(xScreenCenter, textSizeGO);stim.Ep.MVC.GO.x = x_centerCoordinates(xScreenCenter, textSizeGO);
stim.Ep.MVC.GO.y = y_coordinates(upperBorder, visibleYsize, 9/10, textSizeGO);
stim.Ep.MVC.GO.colour = white;

% post-effort rest
switch langage
    case 'fr'
        stim.MVC_rest.text = 'Reposez-vous quelques secondes.';
    case 'engl'
        stim.MVC_rest.text = 'Rest for a few seconds.';
end
[~,~,textSizeRest] = DrawFormattedText(window, stim.MVC_rest.text, 'center', 'center', white);
stim.MVC_rest.x = x_centerCoordinates(xScreenCenter, textSizeRest);
stim.MVC_rest.y = y_coordinates(upperBorder, visibleYsize, 4/5, textSizeRest);
stim.MVC_rest.colour = white;

% post-main task MVC calibration instructions
switch langage
    case 'fr'
        stim.postTaskMVCmeasurement.text = ['Pour finir cette session, nous allons vous demander ',...
            'd''essayer a nouveau de battre votre record.'];
    case 'engl'
        stim.postTaskMVCmeasurement.text = ['To end this session, ',...
            'we are going to ask you to try to beat your record again.'];
end
[~,~,textSizePostTaskMVC] = DrawFormattedText(window, stim.postTaskMVCmeasurement.text,...
    'center', 'center', white, wrapat);
stim.postTaskMVCmeasurement.x = x_centerCoordinates(xScreenCenter, textSizePostTaskMVC);
stim.postTaskMVCmeasurement.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizePostTaskMVC);
stim.postTaskMVCmeasurement.colour = white;

%% mental calibration
switch langage
    case 'fr'
        stim.mentalCalibInstructions.text = ['Repondez ',...
            'aussi vite et ',...
            'aussi correctement que possible. Chaque essai sera remunere.'];
    case 'engl'
        stim.mentalCalibInstructions.text = ['Answer ',...
            'as quickly and correctly as possible. Each trial will be rewarded.'];
end
[~,~,textSizeMentalCalibInstructions] = DrawFormattedText(window, stim.mentalCalibInstructions.text,...
    'center', 'center', white, wrapat);
stim.mentalCalibInstructions.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibInstructions);
stim.mentalCalibInstructions.y = y_coordinates(upperBorder, visibleYsize, 1/3, textSizeMentalCalibInstructions);
stim.mentalCalibInstructions.colour = white;

% calibration feedback
% success
switch langage
    case 'fr'
        [~,~,textSizeMentalCalibSuccess] = DrawFormattedText(window, ['Bravo vous avez tout resolu dans le temps imparti! ',...
            'Votre meilleur temps est de 0.000 s.'],'center', 'center', white, wrapat);
    case 'engl'
        [~,~,textSizeMentalCalibSuccess] = DrawFormattedText(window, ['Well done, you solved everything in the allotted time! ',...
            'Your best timing is 0.000 s '],'center', 'center', white, wrapat);
end
stim.mentalCalibSuccessFbk.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibSuccess);
stim.mentalCalibSuccessFbk.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibSuccess);
stim.mentalCalibSuccessFbk.colour = white;
% success bis
switch langage
    case 'fr'
        stim.mentalCalibSuccessFbk_bis.text = 'Bravo vous avez tout resolu dans le temps imparti!';
    case 'engl'
        stim.mentalCalibSuccessFbk_bis.text = 'Well done, you solved everything in the allotted time!';
end
[~,~,textSizeMentalCalibSuccess_bis] = DrawFormattedText(window, stim.mentalCalibSuccessFbk_bis.text,...
    'center', 'center', white, wrapat);
stim.mentalCalibSuccessFbk_bis.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibSuccess_bis);
stim.mentalCalibSuccessFbk_bis.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibSuccess_bis);
stim.mentalCalibSuccessFbk_bis.colour = white;
% failure (not enough good answers)
switch langage
    case 'fr'
        stim.mentalCalibFailureFbk.text = 'Nous allons refaire cet essai. Essayez de faire mieux!';
    case 'engl'
        stim.mentalCalibFailureFbk.text = 'We will do this trial again. Try to do better!';
end
[~,~,textSizeMentalCalibFail] = DrawFormattedText(window, stim.mentalCalibFailureFbk.text,...
    'center', 'center', white, wrapat);
stim.mentalCalibFailureFbk.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibFail);
stim.mentalCalibFailureFbk.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibFail);
stim.mentalCalibFailureFbk.colour = white;
% failure (too many errors)
switch langage
    case 'fr'
        stim.mentalCalibFailureTooManyErrorsFbk.text = 'Trop d''erreurs! Nous allons refaire cet essai. Essayez de faire mieux!';
    case 'engl'
        stim.mentalCalibFailureTooManyErrorsFbk.text = 'Too many errors! We will do this trial again. Try to do better!';
end
[~,~,textSizeMentalCalibFailTooManyErrors] = DrawFormattedText(window, stim.mentalCalibFailureTooManyErrorsFbk.text,...
    'center', 'center', white, wrapat);
stim.mentalCalibFailureTooManyErrorsFbk.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibFailTooManyErrors);
stim.mentalCalibFailureTooManyErrorsFbk.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibFailTooManyErrors);
stim.mentalCalibFailureTooManyErrorsFbk.colour = white;


% number version
switch langage
    case 'fr'
        [~,~,textSizeMentalCalibFbk] = DrawFormattedText(window, 'Bravo! Votre meilleur score jusque-la est de X bonnes reponses.','center', 'center', white, wrapat);
    case 'engl'
        [~,~,textSizeMentalCalibFbk] = DrawFormattedText(window, 'Well done! Your best score until now is X correct answers.','center', 'center', white, wrapat);
end
stim.mentalCalibFbk.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibFbk);
stim.mentalCalibFbk.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibFbk);
stim.mentalCalibFbk.colour = white;

% end of calibration
switch langage
    case 'fr'
        [~,~,textSizeMentalCalibEnd] = DrawFormattedText(window, ['Bravo! ',...
            'Votre meilleur temps est de 0.000 s.'],'center', 'center', white, wrapat);
    case 'engl'
        [~,~,textSizeMentalCalibEnd] = DrawFormattedText(window, ['Well done! ',...
            'Your best timing is 0.000 s.'],'center', 'center', white, wrapat);
end
stim.mentalCalibEnd.x = x_centerCoordinates(xScreenCenter, textSizeMentalCalibEnd);
stim.mentalCalibEnd.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizeMentalCalibEnd);
stim.mentalCalibEnd.colour = white;
%% Staircase information

switch langage
    case 'fr'
        stim.staircase.text = 'Attention! Vous allez maintenant jouer pour de l''argent reel.';
    case 'engl'
        stim.staircase.text = 'Warning! You will now play for real money.';
end
[~,~,textSizeStaircaseInfo] = DrawFormattedText(window,stim.staircase.text,...
    'center','center',white, wrapat);
stim.staircase.x = x_centerCoordinates(xScreenCenter, textSizeStaircaseInfo);
stim.staircase.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeStaircaseInfo);
stim.staircase.colour = white;
%% color used to represent the effort signal
% no use of monetary images anymore
stim.difficulty.maxColor        = black;
stim.difficulty.currLevelColor  = difficultyArcColor;
stim.difficulty.ovalWidth       = 3;

%% parameters for trait indicating best performance until now for mental calibration
stim.calibBestUntilNow.color = orange;
arcPosition     = stim.difficulty.middle_center;
stim.calibBestUntilNow.circleRadius    = difficultyRectlinearSize/2;
stim.calibBestUntilNow.xCircleCenter = arcPosition(1) + (arcPosition(3) - arcPosition(1))/2;
stim.calibBestUntilNow.yCircleCenter = arcPosition(2) + (arcPosition(4) - arcPosition(2))/2;
stim.calibBestUntilNow.lineWidth = 3;

%% parameters for trait indicating minimal performance to reach for mental calibration
stim.calibMinToReach.color = red;
arcPosition     = stim.difficulty.middle_center;
stim.calibMinToReach.circleRadius    = difficultyRectlinearSize/2;
stim.calibMinToReach.xCircleCenter = arcPosition(1) + (arcPosition(3) - arcPosition(1))/2;
stim.calibMinToReach.yCircleCenter = arcPosition(2) + (arcPosition(4) - arcPosition(2))/2;
stim.calibMinToReach.lineWidth = 3;

%% choice period
switch langage
    case 'fr'
        stim.choice.choiceQuestion.text = 'QUE PREFEREZ-VOUS?';
        stim.choice.choiceOR.text = 'OU';
    case 'engl'
        stim.choice.choiceQuestion.text = 'WHAT DO YOU PREFER?';
        stim.choice.choiceOR.text = 'OR';
end
[~,~,textSizeChoiceQuestion] = DrawFormattedText(window, stim.choice.choiceQuestion.text, 'center','center',white);
stim.choice.choiceQuestion.x = x_centerCoordinates(xScreenCenter, textSizeChoiceQuestion);
stim.choice.choiceQuestion.y = y_coordinates(upperBorder, visibleYsize, 1/8, textSizeChoiceQuestion);
stim.choice.choiceQuestion.colour = white;
[~,~,textSizeOR] = DrawFormattedText(window, stim.choice.choiceOR.text, 'center','center',white);
stim.choice.choiceOR.x = x_centerCoordinates(xScreenCenter, textSizeOR);
stim.choice.choiceOR.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeOR);
stim.choice.choiceOR.colour = white;

% win option
switch langage
    case 'fr'
        stim.choice.win.text = 'Gagner';
    case 'engl'
        stim.choice.win.text = 'Win';
end
[~,~,textSizeWin] = DrawFormattedText(window,stim.choice.win.text,'center','center',white);
xSizeWin = textSizeWin(3) - textSizeWin(1);
ySizeWin = textSizeWin(4) - textSizeWin(2);
stim.textRectSize.xSizeWin = xSizeWin;
stim.textRectSize.ySizeWin = ySizeWin;
% lose option
switch langage
    case 'fr'
        stim.choice.lose.text = 'Perdre';
    case 'engl'
        stim.choice.lose.text = 'Lose';
end
[~,~,textSizeLose] = DrawFormattedText(window, stim.choice.lose.text, 'center','center',white);
xSizeLose = textSizeLose(3) - textSizeLose(1);
ySizeLose = textSizeLose(4) - textSizeLose(2);
stim.textRectSize.xSizeLose = xSizeLose;
stim.textRectSize.ySizeLose = ySizeLose;
% effort
switch langage
    case 'fr'
        stim.choice.for.text = 'pour';
    case 'engl'
        stim.choice.for.text = 'for';
end
[~,~,textSizeForEffort] = DrawFormattedText(window,stim.choice.for.text,'center','center',white);
xSizeForEffort = textSizeForEffort(3) - textSizeForEffort(1);
ySizeForEffort = textSizeForEffort(4) - textSizeForEffort(2);
stim.textRectSize.xSizeForEffort = xSizeForEffort;
stim.textRectSize.ySizeForEffort = ySizeForEffort;
% extract x/y coordinates for the display of the corresponding text
stim.winRewardText.top_left         = [leftBorder + visibleXsize/4 - xSizeWin/2,            stim.reward.text.top_left_start(2) - ySizeWin*2.5];
stim.winRewardText.top_right        = [leftBorder + visibleXsize*(3/4) - xSizeWin/2,        stim.reward.text.top_right_start(2) - ySizeWin*2.5];
stim.loseRewardText.top_left        = [leftBorder + visibleXsize/4 - xSizeLose/2,           stim.punishment.text.top_left_start(2) - ySizeLose*2.5];
stim.loseRewardText.top_right       = [leftBorder + visibleXsize*(3/4) - xSizeLose/2,       stim.punishment.text.top_right_start(2) - ySizeLose*2.5];
stim.effort_introText.bottom_left   = [leftBorder + visibleXsize/4 - xSizeForEffort/2,      stim.difficulty.below_left(2) - ySizeForEffort];
stim.effort_introText.bottom_right  = [leftBorder + visibleXsize*(3/4) - xSizeForEffort/2,  stim.difficulty.below_right(2)  - ySizeForEffort];
% display of confidence mapping
switch langage
    case 'fr'
        stim.leftSure.text      = 'SUR';
        stim.leftUnsure.text    = 'PEU SUR';
        stim.rightUnsure.text   = 'PEU SUR';
        stim.rightSure.text     = 'SUR';
    case 'engl'
        stim.leftSure.text      = 'SURE';
        stim.leftUnsure.text    = 'NOT SURE';
        stim.rightUnsure.text   = 'NOT SURE';
        stim.rightSure.text     = 'SURE';
end
% left sure
[~,~,textSizeLeftSure] = DrawFormattedText(window,stim.leftSure.text,'center','center',white);
xSizeLeftSure = textSizeLeftSure(3) - textSizeLeftSure(1);
ySizeLeftSure = textSizeLeftSure(4) - textSizeLeftSure(2);
stim.leftSure.x = leftBorder + visibleXsize*(1/4) - xSizeLeftSure*(3/2);
stim.leftSure.y = upperBorder + visibleYsize*(19/20) - ySizeLeftSure/2;
% stim.leftSure.colour = [0 255 0]; % colour corresponding to extreme left button: GREEN
stim.leftSure.colour = [0 0 255]; % colour corresponding to extreme left button: BLUE
% left unsure
[~,~,textSizeLeftUnsure] = DrawFormattedText(window,stim.leftUnsure.text,'center','center',white);
% xSizeLeftUnsure = textSizeLeftUnsure(3) - textSizeLeftUnsure(1);
ySizeLeftUnsure = textSizeLeftUnsure(4) - textSizeLeftUnsure(2);
stim.leftUnsure.x = leftBorder + visibleXsize*(1/4) + xSizeLeftSure/2;
stim.leftUnsure.y = upperBorder + visibleYsize*(19/20) - ySizeLeftUnsure/2;
% stim.leftUnsure.colour = [255 0 0]; % colour corresponding to middle left button: RED
stim.leftUnsure.colour = [255 255 0]; % colour corresponding to middle left button: YELLOW
% right unsure
[~,~,textSizeRightUnsure] = DrawFormattedText(window,stim.rightUnsure.text,'center','center',white);
xSizeRightUnsure = textSizeRightUnsure(3) - textSizeRightUnsure(1);
ySizeRightUnsure = textSizeRightUnsure(4) - textSizeRightUnsure(2);
stim.rightUnsure.x = leftBorder + visibleXsize*(3/4) - xSizeRightUnsure*(3/2);
stim.rightUnsure.y = upperBorder + visibleYsize*(19/20) - ySizeRightUnsure/2;
% stim.rightUnsure.colour = [0 0 255]; % colour corresponding to middle right button: BLUE
stim.rightUnsure.colour = [0 255 0]; % colour corresponding to middle right button: GREEN
% right sure
[~,~,textSizeRightSure] = DrawFormattedText(window,stim.rightSure.text,'center','center',white);
% xSizeRightSure = textSizeRightSure(3) - textSizeRightSure(1);
ySizeRightSure = textSizeRightSure(4) - textSizeRightSure(2);
stim.rightSure.x = leftBorder + visibleXsize*(3/4) + xSizeRightUnsure/2;
stim.rightSure.y = upperBorder + visibleYsize*(19/20) - ySizeRightSure/2;
% stim.rightSure.colour = [255 255 0]; % colour corresponding to extreme right button: YELLOW
stim.rightSure.colour = [255 0 0]; % colour corresponding to extreme right button: RED
%% release buttons message
switch langage
    case 'fr'
        stim.releaseButtonsMsg.text = 'Relachez les boutons svp';
    case 'engl'
        stim.releaseButtonsMsg.text = 'Release the buttons please';
end
[~,~,textSizeReleaseButtons] = DrawFormattedText(scr.window, stim.releaseButtonsMsg.text,'center','center',white);
stim.releaseButtonsMsg.x = x_centerCoordinates(xScreenCenter, textSizeReleaseButtons);
stim.releaseButtonsMsg.y = y_coordinates(upperBorder, visibleYsize, 1/2, textSizeReleaseButtons);
stim.releaseButtonsMsg.colour = white;

%% display of the chosen option
switch langage
    case 'fr'
        stim.chosenOptionMsg.text = 'Vous avez choisi';
    case 'engl'
        stim.chosenOptionMsg.text = 'You selected';
end
[~,~,textSizeChosenMsg] = DrawFormattedText(window, stim.chosenOptionMsg.text,'center','center',white);
stim.chosenOptionMsg.x = x_centerCoordinates(xScreenCenter, textSizeChosenMsg);
stim.chosenOptionMsg.y = y_coordinates(upperBorder, visibleYsize, 3/16, textSizeChosenMsg);

% text when they were too slow and then they are forced to perform the default option
switch langage
    case 'fr'
        stim.noChoiceMadeMsg.text = 'Trop lent!';
    case 'engl'
        stim.noChoiceMadeMsg.text = 'Too slow!';
end
[~,~,textSizeNoChoiceMsg] = DrawFormattedText(window, stim.noChoiceMadeMsg.text,'center','center',white);
stim.noChoiceMadeMsg.x = x_centerCoordinates(xScreenCenter, textSizeNoChoiceMsg);
stim.noChoiceMadeMsg.y = y_coordinates(upperBorder, visibleYsize, 3/16, textSizeNoChoiceMsg);

% place reward amount and difficulty level accordingly
ySizeChosenMsg = textSizeChosenMsg(4) - textSizeChosenMsg(2);

% square surrounding chosen option
stim.chosenOption.squareRect = [leftBorder + visibleXsize*(1/3),...
    upperBorder + visibleYsize*(3/16) + ySizeChosenMsg,...
    leftBorder + visibleXsize*(2/3),...
    upperBorder + visibleYsize*(11/12)];
stim.chosenOption.squareColour = black;
stim.chosenOption.squareWidth = 10;
% dotted lines square surrounding chosen option
lineLength = visibleXsize/30;
stim.chosenOption.dottedSquare.xyLines = [];
for iVerticalLines = (stim.chosenOption.squareRect(2)+lineLength/2):(2*lineLength):(stim.chosenOption.squareRect(4) - lineLength)
    xVerticalLeft = stim.chosenOption.squareRect(1); % same as for square
    yStartVertical = iVerticalLines;
    xVerticalRight = stim.chosenOption.squareRect(3);
    yEndVertical = yStartVertical + lineLength;
    stim.chosenOption.dottedSquare.xyLines = [stim.chosenOption.dottedSquare.xyLines,...
        [xVerticalLeft, xVerticalLeft, xVerticalRight, xVerticalRight;...
        yStartVertical, yEndVertical, yStartVertical, yEndVertical]];
end % vertical lines
for iHorizontalLines = (stim.chosenOption.squareRect(1)+lineLength/2):(2*lineLength):(stim.chosenOption.squareRect(3) - lineLength)
    xStartHorizontal = iHorizontalLines;
    yHorizontalTop = stim.chosenOption.squareRect(2); % same as for square
    xEndHorizontal = xStartHorizontal + lineLength;
    yHorizontalBottom = stim.chosenOption.squareRect(4); % same as for square
    stim.chosenOption.dottedSquare.xyLines = [stim.chosenOption.dottedSquare.xyLines,...
        [xStartHorizontal, xEndHorizontal, xStartHorizontal, xEndHorizontal;...
        yHorizontalTop, yHorizontalTop, yHorizontalBottom, yHorizontalBottom]];
end % horizontal lines

% Win/Lose text message
stim.winRewardText.top_center       = [xScreenCenter - xSizeWin/2,  stim.chosenOption.squareRect(2) + ySizeWin*1.5];
stim.loseRewardText.top_center      = [xScreenCenter - xSizeLose/2, stim.chosenOption.squareRect(2) + ySizeWin*1.5];
% amount of money to Win/Lose
stim.chosenOption.reward   = [x_centerCoordinates(xScreenCenter, textSizeR),  stim.winRewardText.top_center(2) + ySizeTextR*1.5]; % chosen option display
stim.reward.text.top_center_start = stim.chosenOption.reward;
stim.chosenOption.punishment   = [x_centerCoordinates(xScreenCenter, textSizeP),  stim.loseRewardText.top_center(2) + ySizeTextP*1.5]; % chosen option display
stim.punishment.text.top_center_start = stim.chosenOption.punishment;
% effort informations
stim.chosenOption.difficulty = stim.difficulty.below_center;
stim.effort_introText.bottom_center = [xScreenCenter - xSizeForEffort/2, stim.difficulty.below_center(2)  - ySizeForEffort];

%% mental effort performance
% display of the relevant instructions
% OR
switch langage
    case 'fr'
        stim.Em.OR.text = 'OU';
    case 'engl'
        stim.Em.OR.text = 'OR';
end
[~,~,textSizeOR] = DrawFormattedText(window, stim.Em.OR.text, 'center', 'center', white );
stim.Em.OR.x = x_centerCoordinates(xScreenCenter, textSizeOR);
stim.Em.OR.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeOR);
% question < 5 or > 5
switch langage
    case 'fr'
        stim.Em.lowerORhigherQuestion.text = 'Chiffre < ou > 5?';
    case 'engl'
        stim.Em.lowerORhigherQuestion.text = 'Is the number < or > than 5?';
end
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
stim.Em.higher.text = '> 5';
[~,~,textSizeHigher] = DrawFormattedText(window,'> 5', 'center', 'center', white );
stim.Em.higher_left.x = leftBorder + visibleXsize*(1/4) - (textSizeHigher(3) - textSizeHigher(1))/2;
stim.Em.higher_right.x = leftBorder + visibleXsize*(3/4) - (textSizeHigher(3) - textSizeHigher(1))/2;
stim.Em.higher.y = y_coordinates(upperBorder, visibleYsize, 5/6, textSizeHigher);
stim.Em.higherInstructions.y = y_coordinates(upperBorder, visibleYsize, 7/8, textSizeHigher);
% press any button
switch langage
    case 'fr'
        stim.Em.pressAnyButtonQuestion.text = 'Appuyer sur n''importe quel bouton';
    case 'engl'
        stim.Em.pressAnyButtonQuestion.text = 'Press any button';
end
[~,~,textSizePressAnyButtonQuestion] = DrawFormattedText(window, stim.Em.pressAnyButtonQuestion.text,...
    'center', 'center', white);
stim.Em.pressAnyButtonQuestion.x = x_centerCoordinates(xScreenCenter, textSizePressAnyButtonQuestion);
stim.Em.pressAnyButtonQuestion.y = y_coordinates(upperBorder, visibleYsize, 1/6, textSizePressAnyButtonQuestion);
% press
switch langage
    case 'fr'
        stim.Em.pressAnyButton.text = 'Appuyer';
    case 'engl'
        stim.Em.pressAnyButton.text = 'Press';
end
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
    stim.Em.(['numberPerfLearning_',n_str]).y = y_coordinates(upperBorder, visibleYsize, 1/10, textSizeEmNumber);
end
Screen('TextSize', window, scr.textSize.baseline);

%% display of the amount of money associated to the current trial
% left win (maximal performance amount)
[~,~,textSizeMaxPerfMoneyWin] = DrawFormattedText(window, '+0.00',...
    'center', 'center', white);
stim.leftMoneyWinEperf.x = x_centerCoordinates(xScreenCenter - difficultyRectlinearSize/2, textSizeMaxPerfMoneyWin);
yCoordMoneyPerf = upperBorder + visibleYsize*(1/2) - (textSizeMaxPerfMoneyWin(4) - textSizeMaxPerfMoneyWin(2))/2 - difficultyRectlinearSize;
if (yCoordMoneyPerf < upperBorder) || (yCoordMoneyPerf > (upperBorder + visibleYsize))
    error('wtf with these coordinates?');
end
stim.leftMoneyWinEperf.y = yCoordMoneyPerf;
% left loss (maximal performance amount)
[~,~,textSizeMaxPerfMoneyLoss] = DrawFormattedText(window, '-0.00',...
    'center', 'center', white);
stim.leftMoneyLoseEperf.x = x_centerCoordinates(xScreenCenter - difficultyRectlinearSize/2, textSizeMaxPerfMoneyLoss);
stim.leftMoneyLoseEperf.y = yCoordMoneyPerf;
% right win (minimal performance amount)
[~,~,textSizeMinPerfMoneyWin] = DrawFormattedText(window, '+0.00',...
    'center', 'center', white);
stim.rightMoneyWinEperf.x = x_centerCoordinates(xScreenCenter + difficultyRectlinearSize/2, textSizeMinPerfMoneyWin);
stim.rightMoneyWinEperf.y = yCoordMoneyPerf;
% right lose (minimal performance amount)
[~,~,textSizeMinPerfMoneyLose] = DrawFormattedText(window, '-0.00',...
    'center', 'center', white);
stim.rightMoneyLoseEperf.x = x_centerCoordinates(xScreenCenter + difficultyRectlinearSize/2, textSizeMinPerfMoneyLose);
stim.rightMoneyLoseEperf.y = yCoordMoneyPerf;

%% prepare feedback messages
% reward feedback
switch langage
    case 'fr'
        stim.feedback.reward.text = 'Vous avez obtenu';
    case 'engl'
        stim.feedback.reward.text = 'You got';
end
[~,~,textSizeRewardFbkMsg] = DrawFormattedText(window, stim.feedback.reward.text,...
    'center', 'center',...
    white);
stim.feedback.reward.x = x_centerCoordinates(xScreenCenter, textSizeRewardFbkMsg);
stim.feedback.reward.y = y_coordinates(upperBorder, visibleYsize, 3/8, textSizeRewardFbkMsg);
stim.feedback.colour = white;

% punishment feedback
switch langage
    case 'fr'
        stim.feedback.punishment.text = 'Vous avez perdu';
    case 'engl'
        stim.feedback.punishment.text = 'You lost';
end
[~,~,textSizePunishmentFbkMsg] = DrawFormattedText(window, stim.feedback.punishment.text,...
    'center', 'center',...
    white);
stim.feedback.punishment.x = x_centerCoordinates(xScreenCenter, textSizePunishmentFbkMsg);
stim.feedback.punishment.y = y_coordinates(upperBorder, visibleYsize, 3/8, textSizePunishmentFbkMsg);

% error: too slow feedback
switch langage
    case 'fr'
        stim.feedback.error_tooSlow.text = 'Trop lent!';
    case 'engl'
        stim.feedback.error_tooSlow.text = 'Too slow!';
end
[~,~,textSizeErrorTooSlowFbkMsg] = DrawFormattedText(window, stim.feedback.error_tooSlow.text,...
    'center', 'center',...
    white);
stim.feedback.error_tooSlow.x = x_centerCoordinates(xScreenCenter, textSizeErrorTooSlowFbkMsg);
stim.feedback.error_tooSlow.y = y_coordinates(upperBorder, visibleYsize, 3/8, textSizeErrorTooSlowFbkMsg);

% error too many errors feedback
switch langage
    case 'fr'
        stim.feedback.error_tooManyErrors.text = 'Trop d''erreurs!';
    case 'engl'
        stim.feedback.error_tooManyErrors.text = 'Too many errors!';
end
[~,~,textSizeErrorTooManyErrorsFbkMsg] = DrawFormattedText(window, stim.feedback.error_tooManyErrors.text,...
    'center', 'center',...
    white);
stim.feedback.error_tooManyErrors.x = x_centerCoordinates(xScreenCenter, textSizeErrorTooManyErrorsFbkMsg);
stim.feedback.error_tooManyErrors.y = y_coordinates(upperBorder, visibleYsize, 3/8, textSizeErrorTooManyErrorsFbkMsg);


% error try again feedback, can be displayed with too many errors and too slow
switch langage
    case 'fr'
        stim.feedback.error_tryAgain.text = 'Concentrez-vous et reessayez!';
    case 'engl'
        stim.feedback.error_tryAgain.text = 'Focus and try again!';
end
[~,~,textSizeErrorTryAgainFbkMsg] = DrawFormattedText(window, stim.feedback.error_tryAgain.text,...
    'center', 'center',...
    white);
stim.feedback.error_tryAgain.x = x_centerCoordinates(xScreenCenter, textSizeErrorTryAgainFbkMsg);
stim.feedback.error_tryAgain.y = y_coordinates(upperBorder, visibleYsize, 4/8, textSizeErrorTryAgainFbkMsg);

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
switch langage
    case 'fr'
        remainingTimeText = 'Temps restant';
    case 'engl'
        remainingTimeText = 'Remaining time';
end
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

function[x] = x_centerCoordinates(xCenter, textSize)
% to center the coordinates on the screen
x = xCenter - (textSize(3) - textSize(1))/2;

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