function[GLMprm] = which_GLM_LGC_pilot(GLM)
% [GLMprm] = which_GLM_LGC_pilot()
% GLM parameters
%
% INPUTS
% GLM: GLM number
%
% OUTPUTS
% GLMprm: GLM parameters

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
end

end % function