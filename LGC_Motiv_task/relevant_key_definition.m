function [key, dq] = relevant_key_definition(effort_type, IRM)
%[key, dq] = relevant_key_definition(effort_type, IRM)
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
% OUTPUTS
% key: structure with mapping for the keys to consider for the experiment
%
% dq: (physical effort task only) structure with relevant info for
% dynamometer recording with a NI card

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
end

%% keyboard keys configuration + waiting and recording first TTL for fMRI
if IRM == 0
    %% key configuration
    KbName('UnifyKeyNames');
    key.left = KbName('LeftArrow');
    key.right = KbName('RightArrow');
    key.space = KbName('Space');
    key.escape= KbName('escape');
elseif IRM == 1
    %% fMRI key configuration
    KbName('UnifyKeyNames');
    key.left = 49;        % 49 % DROITE  bleu, left press
    key.right = 50;       %50   %% GAUCHE JAUNE
    key.space = KbName('Space');
    key.escape = KbName('escape');
    key.trigger_id = 84; % trigger value corresponding to the TTL code
end

end % function