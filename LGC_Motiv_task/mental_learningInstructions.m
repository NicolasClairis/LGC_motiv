function[onset_Press] = mental_learningInstructions(scr, learning_col, learning_instructions, mentalE_prm)
% [onset_Press] = mental_learningInstructions(scr, learning_col, learning_instructions, mentalE_prm)
% mental_learningInstructions will display instructions before learning starts.
%
% INPUTS
% scr: structure with screen parameters
%
% learning_col:
% 'col1': learning with colour 1 only
% 'col2': learning with colour 2 only
% 'all': learning with both colours
%
% learning_instructions
% 'fullInstructions': display instructions: ask if odd/even (or lower/higher than 5) and
% display also on the screen the relevant answer to each question
% 'partialInstructions': display only the two possible answers but not the
% question anymore
% 'noInstructions': no reminder of what the question is nor of where you should answer
%
% mentalE_prm: structure with mental effort parameters
%   .mental_n_col: field with colour for each type of question
%
% OUTPUTS
% onset_Press: time when participant finished reading instructions.


%% extract and define main variables of interest
window = scr.window;
xScreenCenter = scr.xCenter;
yScreenCenter = scr.yCenter;

t_instructions = 5; % minimal amount to see the instructions


%% display color associated to each question type: pair/impair vs low/high than 5
oddEven_col = mentalE_prm.mental_n_col.oddEven;
lowHigh_col = mentalE_prm.mental_n_col.lowHigh;

blackCol = scr.colours.black;
wrapat = scr.wrapat; % go to line for DrawFormattedText when more than this number of characters

%% display general instruction on the screen
% adapt text depending on if instructions included or not

for iTimeLoop = 1:2
    switch learning_instructions
        case 'fullInstructions' % WITH INSTRUCTIONS
            DrawFormattedText(window,...
                ['Vous allez pouvoir vous entraîner à effectuer la tâche. ',...
                'La couleur du chiffre représente la question posée. ',...
                'Vous pouvez répondre avec les boutons. ',...
                'La correspondance entre la position des boutons et la ',...
                'réponse que vous souhaitez donner est la suivante:'],...
                'center', yScreenCenter/3, blackCol, wrapat);
        case 'partialInstructions'
            DrawFormattedText(window,...
                ['Vous allez pouvoir vous entraîner à effectuer la tâche mais vous ',...
                'devrez vous rappeler de la correspondance entre les couleurs et la ',...
                'question à laquelle il faut répondre. Pour rappel:'],...
                'center', yScreenCenter/3, blackCol, wrapat);
        case 'noInstructions'
            DrawFormattedText(window,...
                ['Vous allez pouvoir vous entraîner à effectuer la tâche mais vous ',...
                'devrez vous rappeler de la correspondance entre les couleurs et la ',...
                'question à laquelle il faut répondre. Pour rappel:'],...
                'center', yScreenCenter/3, blackCol, wrapat);
        case 'extendedLearning'
            DrawFormattedText(window,...
                ['Vous allez à présent vous entraîner à nouveau sur les ',...
                'différents niveaux de difficulté que vous rencontrerez dans la tâche. ',...
                'Pour rappel:'],...
                'center',yScreenCenter/3,blackCol, wrapat);
    end
    
    % odd/even info
    if strcmp(learning_col,'all') ||...
            ( strcmp(learning_col, 'col1') && strcmp(mentalE_prm.mental_n_col.col1,'oddEven')) ||...
            ( strcmp(learning_col, 'col2') && strcmp(mentalE_prm.mental_n_col.col2,'oddEven'))
        DrawFormattedText(window,'Pair ou Impair?',...
            'center',yScreenCenter, oddEven_col);
        DrawFormattedText(window, 'Pair',...
            xScreenCenter/2, yScreenCenter*5/4, oddEven_col);
        DrawFormattedText(window, 'Impair',...
            xScreenCenter*3/2, yScreenCenter*5/4, oddEven_col);
    end
    
    % lower/higher than 5 info
    if strcmp(learning_col,'all') ||...
            ( strcmp(learning_col, 'col1') && strcmp(mentalE_prm.mental_n_col.col1,'lowHigh')) ||...
            ( strcmp(learning_col, 'col2') && strcmp(mentalE_prm.mental_n_col.col2,'lowHigh'))
        DrawFormattedText(window,'< ou > 5?',...
            'center',yScreenCenter*6/4, lowHigh_col);
        DrawFormattedText(window, '< 5',...
            xScreenCenter/2, yScreenCenter*7/4, lowHigh_col);
        DrawFormattedText(window, '> 5',...
            xScreenCenter*3/2, yScreenCenter*7/4, lowHigh_col);
    end
    
    if iTimeLoop == 1 % force them to read at first
        Screen(window,'Flip');
        WaitSecs(t_instructions);
    elseif iTimeLoop == 2 % after t_instructions seconds, they can manually start
        DrawFormattedText(window, 'Appuyez quand vous êtes prêt(e) à commencer la tâche.',...
            'center', yScreenCenter*15/8, blackCol);
        Screen(window,'Flip');
        onset_Press = KbWait;
    end
end % loop over forced reading/manual pass loop

end % function