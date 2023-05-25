function [con_names, con_vec] = LGC_RL_load_con(GLMprm, sub_nm)
%[con_names, con_vec] = LGC_RL_load_con(GLMprm, sub_nm)
% LGC_RL_load_con extracts the name (in con_names) and the corresponding
% vectors (in con_vec) of each regressor that will be applied in the first
% and the second level.
%
% INPUTS
% GLMprm: GLM parameters (see which_GLM_MS2.m)
% 
% sub_nm: subject identification full name
%
% OUTPUTS
% con_names: list of contrast names in the order they will be applied at
% the first and second level.
%
% con_vec: contrast vectors to be used at the first level, their order
% follows con_names
%
% See also contrasts_megaconcatenation_MS2_NicoC

%% extract number of regressors of non-interest (movement + tooslow/too fast trials eventually)
% [nb_nonI] = MS2_nreg_noninterest_estimation(sub_nm, GLMprm);

%% contrast index
jCon = 0;

%% extract number of regressors to consider + index for each one of them
[ n_prm, prm_idx ] = LGC_RL_First_level_n_prm( GLMprm, sub_nm );
n_regs = n_prm.all;

%% prepare contrasts for each task
single_cons = fieldnames(prm_idx);
n_single_con = length(single_cons);

con_names   = cell(1,n_single_con*2);
con_vec     = zeros(n_single_con*2, n_regs);

for iSingleCon = 1:n_single_con
    
    % initialize contrast
    curr_vec_nm = single_cons{iSingleCon};
    curr_vec_idx = prm_idx.(curr_vec_nm);
    
    % positive contrast
    jCon = jCon + 1;
    con_names{jCon} = [curr_vec_nm,'_pos'];
    con_vec(jCon, curr_vec_idx) = 1;
    
    % negative contrast
    jCon = jCon + 1;
    con_names{jCon} = [curr_vec_nm,'_neg'];
    con_vec(jCon, curr_vec_idx) = -1;
    
end % regressor loop

end % function