function [key] = relevant_key_definition(effort_type, IRM, n_buttonsChoice)
%[key] = relevant_key_definition(effort_type, IRM, n_buttonsChoice)
% relevant_key_definition defines keyboard, handgrip and TTL relevant
% code.
%
% INPUTS
% effort_type:
% 'physical': physical effort task (will need handgrip on top of keyboards
% to answer to the choice)
% 'mental': mental effort task
%
% IRM:
% (0) use keyboard arrows
% (1) use button pad from the MRI
%
% n_buttonsChoice: define the number of buttons to use for the choice (2 or
% 4?) and define buttons accordingly
%
% OUTPUTS
% key: structure with mapping for the keys to consider for the experiment
%

%% keyboard keys configuration + waiting and recording first TTL for fMRI
if IRM == 0
    %% key configuration
    KbName('UnifyKeyNames');
    switch n_buttonsChoice
        case 2
            key.left = KbName('LeftArrow');
            key.right = KbName('RightArrow');
        case 4
            key.leftSure    = KbName('1');
            key.leftUnsure  = KbName('2');
            key.rightUnsure = KbName('3');
            key.rightSure   = KbName('4');
            % also need the left/right buttons for the mental effort task
            %             if strcmp(effort_type,'mental')
            key.left = KbName('LeftArrow');
            key.right = KbName('RightArrow');
            %             end
    end
    key.space = KbName('Space');
    key.escape= KbName('escape');
elseif IRM == 1
    %% fMRI key configuration
    KbName('UnifyKeyNames');
    switch n_buttonsChoice
        case 2
            key.left = 66; % LEFT BUTTON = blue, letter 'b'
            key.right = 90; % RIGHT BUTTON = yellow, letter 'z' with swiss or english keyboard
            %     key.right = 89; % RIGHT BUTTON = yellow, letter 'y' with french keyboard
        case 4
            key.leftSure    = 71; % green button (g)
            key.leftUnsure  = 82; % red button (r)
            key.rightUnsure = 66; % blue button (b)
            key.rightSure   = 90; % or 89 depending on the language of the keyboard
            % corresponds to 'y' or 'z' (for yellow)
            
            % also need the left/right buttons for the mental effort task
            key.left = 71;
            key.right = 82;
    end
    key.space = KbName('Space');
    key.escape = KbName('escape');
    key.trigger_id = 84; % trigger value corresponding to the TTL code (letter 't')
end

%% store how many buttons to answer there are
key.n_buttonsChoice = n_buttonsChoice;

end % function