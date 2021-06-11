function[onset_Press] = mental_learningInstructions(scr, stim, learning_col, learning_instructions, mentalE_prm)
% [onset_Press] = mental_learningInstructions(scr, stim, learning_col, learning_instructions, mentalE_prm)
% mental_learningInstructions will display instructions before learning starts.
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with stimuli informations
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
% 'extendedLearning': extended learning, to train on each level of
% difficulty
%
% mentalE_prm: structure with mental effort parameters
%   .mental_n_col: field with colour for each type of question
%
% OUTPUTS
% onset_Press: time when participant finished reading instructions.


%% extract and define main variables of interest
window = scr.window;

t_instructions = 5; % minimal amount to see the instructions


%% display color associated to each question type: pair/impair vs low/high than 5
oddEven_col = mentalE_prm.mental_n_col.oddEven;
lowHigh_col = mentalE_prm.mental_n_col.lowHigh;

% select side of each question
sideQuestion = mentalE_prm.sideQuestion;

wrapat = scr.wrapat; % go to line for DrawFormattedText when more than this number of characters

%% display general instruction on the screen
% adapt text depending on if instructions included or not

for iTimeLoop = 1:2
    % display instructions
    DrawFormattedText(window,...
        stim.Em.learning.(learning_instructions).text,...
        stim.Em.learning.(learning_instructions).x, stim.Em.learning.(learning_instructions).y,...
        stim.Em.learning.(learning_instructions).colour, wrapat);
    
    % odd/even info
    if strcmp(learning_col,'all') ||...
            ( strcmp(learning_col, 'col1') && strcmp(mentalE_prm.mental_n_col.col1,'oddEven')) ||...
            ( strcmp(learning_col, 'col2') && strcmp(mentalE_prm.mental_n_col.col2,'oddEven'))
        DrawFormattedText(window, stim.Em.oddORevenQuestion.text,...
            stim.Em.oddORevenQuestion.x, stim.Em.oddORevenQuestionInstructions.y, oddEven_col);
        if sideQuestion.oE.pair == -1 && sideQuestion.oE.impair == +1
            x_pair      = stim.Em.even_left.x;
            x_impair    = stim.Em.odd_right.x;
        elseif sideQuestion.oE.pair == +1 && sideQuestion.oE.impair == -1
            x_pair      = stim.Em.even_right.x;
            x_impair    = stim.Em.odd_left.x;
        else
            error('error in sideQuestion definition');
        end
        DrawFormattedText(window, stim.Em.even.text, x_pair, stim.Em.evenInstructions.y, oddEven_col );     % pair
        DrawFormattedText(window, stim.Em.odd.text, x_impair, stim.Em.oddInstructions.y, oddEven_col );     % impair
    end
    
    % lower/higher than 5 info
    if strcmp(learning_col,'all') ||...
            ( strcmp(learning_col, 'col1') && strcmp(mentalE_prm.mental_n_col.col1,'lowHigh')) ||...
            ( strcmp(learning_col, 'col2') && strcmp(mentalE_prm.mental_n_col.col2,'lowHigh'))
        DrawFormattedText(window, stim.Em.lowerORhigherQuestion.text,...
            stim.Em.lowerORhigherQuestion.x, stim.Em.lowerORhigherQuestionInstructions.y, lowHigh_col);
        if sideQuestion.hL.low == -1 && sideQuestion.hL.high == +1
            x_low = stim.Em.lower_left.x;
            x_high = stim.Em.higher_right.x;
        elseif sideQuestion.hL.low == +1 && sideQuestion.hL.high == -1
            x_low = stim.Em.lower_right.x;
            x_high = stim.Em.higher_left.x;
        else
            error('error in sideQuestion definition');
        end
        DrawFormattedText(window, stim.Em.lower.text, x_low, stim.Em.lowerInstructions.y, lowHigh_col );    % < 5
        DrawFormattedText(window, stim.Em.higher.text, x_high, stim.Em.higherInstructions.y, lowHigh_col ); % > 5
    end
    
    if iTimeLoop == 1 % force them to read at first
        Screen(window,'Flip');
        WaitSecs(t_instructions);
    elseif iTimeLoop == 2 % after t_instructions seconds, they can manually start
        DrawFormattedText(window, stim.pressWhenReady.text,...
            stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour);
        Screen(window,'Flip');
        onset_Press = KbWait;
    end
end % loop over forced reading/manual pass loop

end % function