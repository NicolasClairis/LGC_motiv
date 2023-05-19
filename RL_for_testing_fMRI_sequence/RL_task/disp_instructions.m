function [timeDisplay] = disp_instructions(window, x,y, IRM, taskName)
%disp_instructions common function for all fMRI/MBB battery scripts to
%display instructions on screen
%
% INPUTS
% window: window where to display instructions
% x,y: coordinates of the center of the screen
% IRM: IRM (1) or outside IRM (may be training or other condition) (0)? => timings and instructions vary
% depending on which condition
% taskName: different instructions depending on task => taskName says which
% task we are at now
%
% made by Nicolas Clairis - march 2017

% instructions
if strcmp(taskName,'gripRP')
    if IRM == 0
        instruction = ['Vous allez devoir serrer une poignee. ',...
            'Plus vous serrez fort, plus vos gains seront importants ou vos pertes minimes. ',...
            'Une image a gauche de l''ecran vous indiquera le montant en jeu. ',...
            'Attendez que l''echelle apparaisse avant de demarrer votre effort. ',...
            'Un trait rouge indiquera le maximum de force atteint jusque la et le gain correspondant. ',...
            'La barre orange indique votre force actuelle. ',...
            'Seule le maximum de force exercee compte, et non pas la duree de l''effort.'];
    elseif IRM == 1
        instruction = 'Serrez la poignee pour maximiser vos gains!';
    end
elseif strcmp(taskName,'mentalRP')
    if IRM == 0
        instruction = ['Vous allez devoir comparer des paires de chiffres. ',...
            'Chaque chiffre d''une paire a une valeur et une taille differentes. ',...
            'Vous devez selectionner les chiffres qui ont la plus grande valeur numerique en ignorant les differences de tailles. ',...
            'Pour cela, vous devez utiliser les fleches du clavier. ',...
            'Les chiffres d''une meme paire sont toujours sur la meme ligne horizontale. Il faut toujours demarrer la comparaison en partant du bas. ',...
            'Une image a gauche de l''ecran vous indiquera le montant en jeu. ',...
            'Une barre orange surmontee d''un trait rouge indiquera combien de paires vous avez resolu jusque-la et le gain correspondant. ',...
            'Votre remuneration dependra du nombre de paires resolues dans le temps imparti. ',...
            'Il y a toujours 10 paires et c''est a vous de les resoudre le plus vite que possible. ',...
            'Si vous faites une erreur, vous aurez un court delai de penalite (0.5 sec), pendant lequel vous ne pourrez plus appuyer.'];
    elseif IRM == 1
        instruction = 'Resolvez le plus de paires que possible pour maximiser vos gains';
    end
elseif strcmp(taskName,'learning')
    if IRM == 0
        instruction = ['Vous allez devoir faire des choix entre differents symboles. ',...
            'Deux symboles vous seront presentes a chaque essai, de part et d''autre d''une croix de fixation. ',...
            'Pour selectionner un symbole, vous devez utiliser les fleches gauche et droite du clavier. ',...
            'Votre choix sera enregistre au bout d''un delai de 3 secondes. ',...
            'Il faut donc maintenir la touche enfoncee jusqu''a ce que le symbole choisi soit entoure d''un carre rouge. ',...
            'Si aucune touche n''est enfoncee a la fin du delai, l''essai sera considere comme un echec et vous perdrez de l''argent. ',...
            'Vous verrez ensuite le resultat de votre choix: ce peut etre soit un gain (+10euros), soit une perte (-10euros), soit rien (0euro). ',...
            'Certains symboles font gagner (+10euros), ou perdre (-10euros) plus souvent que d''autres. ',...
            'Vous devez donc trouver par essai / erreur quels sont les symboles les plus avantageux pour vous. ',...
            'Sachant que vous ne pourrez pas gagner a tous les coups. ',...
            'La signification de chaque symbole ne change pas en cours de route, et ne depend pas de la position a l''ecran (gauche ou droite).'];
    elseif IRM == 1
        instruction = 'Faites les bons choix pour maximiser vos gains!';
    end
end

% text parameters
if IRM == 0
Screen('TextSize',window,20);
elseif IRM == 1
    Screen('TextSize',window,30);
end
white = [255 255 255];
wrapat = 60; % number of characters before breaking the line
vspacing = 2; % space between lines

% position of display
DrawFormattedText(window, instruction, 'center', 'center', white, wrapat, [], [], vspacing);

% add "press any button to continue" text to screen
if IRM == 0
    textPleasePress = 'appuyer sur une touche pour continuer';
    Screen('TextSize', window, 20);
    [width,hight] = RectSize(Screen('TextBounds',window,textPleasePress));
    Screen('DrawText',window, textPleasePress, x-width/2, 9*y/5, [100 100 100]);
end

% display on screen + extract timing when displayed
[~,timeDisplay,~,~,~] = Screen(window,'Flip');

end

