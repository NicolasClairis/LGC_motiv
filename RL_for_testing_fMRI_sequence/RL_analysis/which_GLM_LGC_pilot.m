function[GLMprm] = which_GLM_LGC_pilot(GLM)
% [GLMprm] = which_GLM_LGC_pilot()
% GLM parameters
%
% INPUTS
% GLM: GLM number
%
% OUTPUTS
% GLMprm: GLM parameters

%% initialize everything to 0
%% list of all possible onsets and regressors
[list_all_GLMprm, n_potReg] = LGC_list_potential_onsets_and_modulators;

%% set all parameters to zero by default

% general
for iGal = 1:n_potReg.gal
    curr_gal_prm = list_all_GLMprm.gal{iGal};
    GLMprm.gal.(curr_gal_prm) = 0;
end
GLMprm.gal.orth_vars = 1; % default = orthogonalize variables

% RL
% onsets & durations
for iORL = 1:n_potReg.RL.onsets
    curr_RL_onset_prm   = list_all_GLMprm.RL.onsets{iORL};
    curr_RL_dur_prm     = list_all_GLMprm.RL.durations{iORL};
    GLMprm.(curr_RL_onset_prm) = 0;
    GLMprm.(curr_RL_dur_prm) = 0;
end
% loop through modulators
for i_RL_mod = 1:n_potReg.RL.mods
    curr_mod = list_all_GLMprm.RL.modulators{i_RL_mod};
    for iRL_regs = 1:n_potReg.RL.(curr_mod) % loop through regressors
        curr_RL_prm = list_all_GLMprm.RL.(curr_mod){iRL_regs};
        GLMprm.(curr_mod).(curr_RL_prm) = 0;
    end
end

% special case for model type which is a string and not a numerical value
GLMprm.mod_stim.mdl_type    = '';
GLMprm.mod_chosen.mdl_type  = '';
GLMprm.mod_fbk.mdl_type     = '';

%% set up main attributes
switch GLM
    case 1
        % RL: Val = pAQA+pBQB/Conf = [p(left)-0.5]^2/DT
        GLMprm.gal.add_drv = 0;
        GLMprm.gal.grey_mask = 0;
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zscore = 1;
        % RL
        GLMprm.o_stim = 3;
        GLMprm.dur_stim = 1;
        GLMprm.mod_stim.mdl_type = 'Nico';
        GLMprm.mod_stim.mdl_n = 6;
        GLMprm.mod_stim.SV = 1;
        GLMprm.mod_stim.pBest = 2;
        GLMprm.mod_stim.RT = 1;
        GLMprm.o_chosen = 1;
        GLMprm.dur_chosen = 1;
        GLMprm.o_fbk = 5;
        GLMprm.dur_fbk = 0;
        GLMprm.mod_fbk.mdl_type = GLMprm.mod_stim.mdl_type;
        GLMprm.mod_fbk.mdl_n = GLMprm.mod_stim.mdl_n; % same as for stim
        GLMprm.mod_fbk.PE = 1;
    case 2 % same as GLM 1 but with grey matter mask
        % RL: Val = pAQA+pBQB/Conf = [p(left)-0.5]^2/DT
        GLMprm.gal.add_drv = 0;
        GLMprm.gal.grey_mask = 1;
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zscore = 1;
        % RL
        GLMprm.o_stim = 3;
        GLMprm.dur_stim = 1;
        GLMprm.mod_stim.mdl_type = 'Nico';
        GLMprm.mod_stim.mdl_n = 6;
        GLMprm.mod_stim.SV = 1;
        GLMprm.mod_stim.pBest = 2;
        GLMprm.mod_stim.RT = 1;
        GLMprm.o_chosen = 1;
        GLMprm.dur_chosen = 1;
        GLMprm.o_fbk = 5;
        GLMprm.dur_fbk = 0;
        GLMprm.mod_fbk.mdl_type = GLMprm.mod_stim.mdl_type;
        GLMprm.mod_fbk.mdl_n = GLMprm.mod_stim.mdl_n; % same as for stim
        GLMprm.mod_fbk.PE = 1;
end

end % function