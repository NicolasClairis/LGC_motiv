function [key, dq] = relevant_key_definition(effort_type, IRM, n_buttonsChoice)
%[key, dq] = relevant_key_definition(effort_type, IRM, n_buttonsChoice)
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
% dq: (physical effort task only) structure with relevant info for
% dynamometer recording with a NI card

% define button mode
button_mode = '4buttons_numeric'; % '2buttons_letters'/'4buttons_letters'/'4buttons_numeric'

%% if physical effort task => requires opening the connection to the BioPac
switch effort_type
    case 'physical'
        % open communication with National Instruments card (which should
        % be connected to the Biopac module, itself connected to the
        % handgrip. Be careful to set everything properly in the channel 1
        % as input (grip) and output (to NI card)
        
        % you can use the following functions to check the name of the
        % device
        % d = daqlist;
        % d{1, "DeviceInfo"}
        
        % initialize NI input
        dq = daq("ni");
        dq.Rate = 20; % define the acquisition rate of the module
        NI_module_nm = "cDAQ1Mod1"; % define module to use on the NI card
        NI_channel_output = "ai0"; % define the channel output from the NI card
        addinput(dq, NI_module_nm, NI_channel_output,"Voltage");
    case 'mental'
        dq = []; % need to be initialized even though remains empty
end

%% keyboard keys configuration + waiting and recording first TTL for fMRI
KbName('UnifyKeyNames');
key.escape = KbName('escape');
key.space = KbName('space');
if IRM == 0
    %% key configuration
    switch n_buttonsChoice
        case 2
            key.left = KbName('LeftArrow');
            key.right = KbName('RightArrow');
        case 4
            key.leftSure    = KbName('c');
            key.leftUnsure  = KbName('v');
            key.rightUnsure = KbName('b');
            key.rightSure   = KbName('n');
            % also need the left/right buttons for the mental effort task
            %             if strcmp(effort_type,'mental')
            %             key.left = KbName('LeftArrow');
            %             key.right = KbName('RightArrow');
            key.left = KbName('c');
            key.right = KbName('v');
            %             end
    end
elseif IRM == 1
    %% fMRI key configuration
    switch n_buttonsChoice
        case 2
            key.left = KbName('LeftArrow');
            key.right = KbName('RightArrow');
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