function[dq] = grip_initialization()
%[dq] = grip_initialization()
% grip_initialization opens communication with National Instruments card
% (which should be connected to the Biopac module, itself connected to the
% handgrip.) Be careful to set everything properly in the channel 1
% as input (grip) and output (to NI card).
%
% OUTPUTS
% dq: structure with relevant informations for the handgrip

% you can use the following functions to check the name of the
% device
% d = daqlist;
% d{1, "DeviceInfo"}

%% initialize NI input
dq = daq("ni");
dq.Rate = 20; % define the acquisition rate of the module
NI_module_nm = "cDAQ1Mod1"; % define module to use on the NI card
NI_channel_output = "ai0"; % define the channel output from the NI card
addinput(dq, NI_module_nm, NI_channel_output,"Voltage");

end % function