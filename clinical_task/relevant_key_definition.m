function [key] = relevant_key_definition(IRM, n_buttonsChoice)
%[key] = relevant_key_definition(IRM, n_buttonsChoice)
% relevant_key_definition defines keyboard and TTL relevant code.
%
% INPUTS
% IRM:
% (0) use keyboard arrows
% (1) use button pad from the MRI
%
% n_buttonsChoice: define the number of buttons to use for the choice (2 or
% 4?) and define buttons accordingly
%
% OUTPUTS
% key: structure with mapping for the keys to consider for the experiment

% define button mode
button_mode = '4buttons_numeric'; % '2buttons_letters'/'4buttons_letters'/'4buttons_numeric'

%% keyboard keys configuration + waiting and recording first TTL for fMRI
KbName('UnifyKeyNames');
key.escape = KbName('escape');
key.space = KbName('space');
    %% fMRI key configuration
    switch n_buttonsChoice
        case 2
            switch button_mode
                case '2buttons_letters'
                    key.left = 71; % green button ('g')
                    key.right = 82; % red button ('r')
                case '4buttons_letters'
                    key.left = 66; % blue button (letter 'b')
                    key.right = 90; % yellow button (letter 'z' with swiss or english keyboard)
                    %     key.right = 89; % yellow button (letter 'y' with french keyboard)
                case '4buttons_numeric'
                    key.left = 49; % button (1)
                    key.right = 50; % button (2)
            end
        case 4
            switch button_mode
                case '2buttons_letters'
                    %% old buttons (green/red/blue/yellow) (letters mode)
                    key.leftSure    = 71; % green button (g)
                    key.leftUnsure  = 82; % red button (r)
                    key.rightUnsure = 66; % blue button (b)
                    key.rightSure   = 90; % or 89 depending on the language of the keyboard
                    % corresponds to 'y' or 'z' (for yellow)
                    
                    % also need the left/right buttons for the mental effort task
                    key.left = 71;
                    key.right = 82;
                    
                case '4buttons_letters'
                    %% new buttons (blue/yellow/green/red) (letters mode)
                    key.leftSure    = 66; % blue button (b)
                    key.leftUnsure  = 90; % yellow button (y/z depending on langage used
                    key.rightUnsure = 71; % green button (g)
                    key.rightSure   = 82; % red button (r)
                    
                    % also need the left/right buttons for the mental effort task
                    key.left = 66;
                    key.right = 90;
                    
                case '4buttons_numeric'
                    %% new buttons (blue/yellow/green/red) (numeric mode: 1/2/3/4)
                    key.leftSure    = 49; % blue button (1)
                    key.leftUnsure  = 50; % yellow button (2)
                    key.rightUnsure = 51; % green button (3)
                    key.rightSure   = 52; % red button (4)
                    
                    % also need the left/right buttons for the mental effort task
                    key.left = 49;
                    key.right = 50;
            end
    end
    
    %% fMRI TTL triggers
    if IRM == 1
    switch button_mode
        case {'2buttons_letters','4buttons_letters'}
            %% mode 0 (letters)
            key.trigger_id = 84; % trigger value corresponding to the TTL code (letter 't')
        case '4buttons_numeric'
            %% mode 4 (numeric)
            key.trigger_id = 54; % trigger value corresponding to the TTL code (number '6')
    end
end

%% store how many buttons to answer there are
key.n_buttonsChoice = n_buttonsChoice;

end % function