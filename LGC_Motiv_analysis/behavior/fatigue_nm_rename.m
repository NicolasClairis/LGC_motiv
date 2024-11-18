function[F_nm2] = fatigue_nm_rename(F_nm1)
% [F_nm2] = fatigue_nm_rename(F_nm1)
% fatigue_nm_rename will rename the fatigue variable entered in input to
% make it compatible with a visual output in a Matlab figure.
%
% INPUTS
% F_nm1: string with the name for the fatigue variable to rename
%
% OUTPUTS
% F_nm2: string with the updated name for F_nm1.
%
% See also fatigue_measurements_crosscorrel.m and fatigue_f_metabolism.m
% where this function is used.

switch F_nm1
    case 'MPSTEFS_physical'
        F_nm2 = 'MPSTEFS_p';
    case 'MPSTEFS_mental'
        F_nm2 = 'MPSTEFS_m';
    case 'MPSTEFS_physical_energy'
        F_nm2 = 'MPSTEFS e_p';
    case 'MPSTEFS_mental_energy'
        F_nm2 = 'MPSTEFS e_m';
    case 'MPSTEFS_physical_fatigue'
        F_nm2 = 'MPSTEFS fat_p';
    case 'MPSTEFS_mental_fatigue'
        F_nm2 = 'MPSTEFS fat_m';
    case 'n_covid'
        F_nm2 = 'n.covid';
    case 'F4_F1'
        F_nm2 = 'F4-F1';
    case 'F4_F3'
        F_nm2 = 'F4-F3';
    case 'MVC_R1a'
        F_nm2 = 'MVC r1a';
    case 'MVC_R1b'
        F_nm2 = 'MVC r1b';
    case 'MVC_R2a'
        F_nm2 = 'MVC r2a';
    case 'MVC_R2b'
        F_nm2 = 'MVC r2b';
    case 'MVC_dR1'
        F_nm2 = 'MVC Δr1';
    case 'MVC_dR2'
        F_nm2 = 'MVC Δr2';
    case 'MVC_dTask'
        F_nm2 = 'MVC ΔTask';
    case 'perf_decrease_slope_E0'
        F_nm2 = 'perf dE0';
    case 'perf_decrease_slope_E1'
        F_nm2 = 'perf dE1';
    case 'perf_decrease_slope_E2'
        F_nm2 = 'perf dE2';
    case 'perf_decrease_slope_E3'
        F_nm2 = 'perf dE3';
    case 'choices_kFp'
        F_nm2 = 'kFp';
    otherwise
        F_nm2 = F_nm1;
end % fatigue variable

end % function