function[force_newtons, force_kgf] = grip_biopac_volts_to_newtons_conversion(force_volts)
% [force_newtons, force_kgf] = grip_biopac_volts_to_newtons_conversion(force_volts)
% grip_biopac_volts_to_newtons_conversion will convert the force value in
% volts extracted from the biopac into a force in newtons and in kgf.
%
% The values for the conversion are based on this:
% https://www.biopac.com/wp-content/uploads/TSD121B-MRI.pdf
% and on the conversion from kgf to Newtons (1kgf = 9.80665 N)
%
% INPUTS
% force_volts: value of the force in volts (as 1 value or a 1*n vector)
%
% OUTPUTS
% force_newtons: value of the force in newtons (as 1 value or a 1*n vector)
%
% For info: normal range should be between 0 and 300 N (see (Cramer et al,
% 2002) for example)
%
% force_kgf: value of the force in kgf (as 1 value or a 1*n vector)
% function developed by Nicolas Clairis & Arthur Barakat - 23/01/2023

%% general constants
kgf_to_N_conversion  = 9.80665;
biopac_volt_to_kgf = 3.128;

%% convert volts into kgf
force_kgf = force_volts.*biopac_volt_to_kgf;
%% convert volts into newtons
force_newtons = force_kgf.*kgf_to_N_conversion;

end % function