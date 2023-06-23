function [con_names, con_vec] = LGC_RL_load_con(GLMprm)
%[con_names, con_vec] = LGC_RL_load_con(GLMprm)
% LGC_RL_load_con extracts the name (in con_names) and the corresponding
% vectors (in con_vec) of each regressor that will be applied in the first
% and the second level.
%
% INPUTS
% GLMprm: GLM parameters (see which_GLM_MS2.m)
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
[ n_prm, prm_idx ] = LGC_RL_First_level_n_prm( GLMprm );
n_regs = n_prm.all;

%% prepare contrasts for each task
single_cons = fieldnames(prm_idx);
n_single_con = length(single_cons);

con_names   = cell(1,n_single_con*2+9*2);
con_vec     = zeros(n_single_con*2+9*2, n_regs);

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

%% add comparison between runs for the contrasts of interest (PE, SV and RT)
RPE_r3_nm = 'L_mod_fbk_PE_GL_Pairs_run3';
RPE_r2_nm = 'L_mod_fbk_PE_GL_Pairs_run2';
RPE_r1_nm = 'L_mod_fbk_PE_GL_Pairs_run1';
RPE_r3_idx = prm_idx.(RPE_r3_nm);
RPE_r2_idx = prm_idx.(RPE_r2_nm);
RPE_r1_idx = prm_idx.(RPE_r1_nm);

%% PE TE3 - TE2
% positive contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_fbk_PE_GL_Pairs_r3_min_r2';
con_vec(jCon, RPE_r3_idx) = 1;
con_vec(jCon, RPE_r2_idx) = -1;

% negative contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_fbk_PE_GL_Pairs_r2_min_r3';
con_vec(jCon, RPE_r3_idx) = -1;
con_vec(jCon, RPE_r2_idx) = 1;
    
%% PE TE3 - TE1
% positive contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_fbk_PE_GL_Pairs_r3_min_r1';
con_vec(jCon, RPE_r3_idx) = 1;
con_vec(jCon, RPE_r1_idx) = -1;

% negative contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_fbk_PE_GL_Pairs_r1_min_r3';
con_vec(jCon, RPE_r3_idx) = -1;
con_vec(jCon, RPE_r1_idx) = 1;
%% PE TE2 - TE1
% positive contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_fbk_PE_GL_Pairs_r2_min_r1';
con_vec(jCon, RPE_r2_idx) = 1;
con_vec(jCon, RPE_r1_idx) = -1;

% negative contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_fbk_PE_GL_Pairs_r1_min_r2';
con_vec(jCon, RPE_r2_idx) = -1;
con_vec(jCon, RPE_r1_idx) = 1;
%% Reaction times
RT_r3_nm = 'L_mod_stim_RT_GL_Pairs_run3';
RT_r2_nm = 'L_mod_stim_RT_GL_Pairs_run2';
RT_r1_nm = 'L_mod_stim_RT_GL_Pairs_run1';
RT_r3_idx = prm_idx.(RT_r3_nm);
RT_r2_idx = prm_idx.(RT_r2_nm);
RT_r1_idx = prm_idx.(RT_r1_nm);

%% RT TE3 - TE2
% positive contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_RT_GL_Pairs_r3_min_r2';
con_vec(jCon, RT_r3_idx) = 1;
con_vec(jCon, RT_r2_idx) = -1;

% negative contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_RT_GL_Pairs_r2_min_r3';
con_vec(jCon, RT_r3_idx) = -1;
con_vec(jCon, RT_r2_idx) = 1;
%% RT TE3 - TE1
% positive contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_RT_GL_Pairs_r3_min_r1';
con_vec(jCon, RT_r3_idx) = 1;
con_vec(jCon, RT_r1_idx) = -1;

% negative contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_RT_GL_Pairs_r1_min_r3';
con_vec(jCon, RT_r3_idx) = -1;
con_vec(jCon, RT_r1_idx) = 1;
%% RT TE2 - TE1
% positive contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_RT_GL_Pairs_r2_min_r1';
con_vec(jCon, RT_r2_idx) = 1;
con_vec(jCon, RT_r1_idx) = -1;

% negative contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_RT_GL_Pairs_r1_min_r2';
con_vec(jCon, RT_r2_idx) = -1;
con_vec(jCon, RT_r1_idx) = 1;

%% value
SV_r3_nm = 'L_mod_stim_SV_GL_Pairs_run3';
SV_r2_nm = 'L_mod_stim_SV_GL_Pairs_run2';
SV_r1_nm = 'L_mod_stim_SV_GL_Pairs_run1';
SV_r3_idx = prm_idx.(SV_r3_nm);
SV_r2_idx = prm_idx.(SV_r2_nm);
SV_r1_idx = prm_idx.(SV_r1_nm);

%% value TE3 - TE2
% positive contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_SV_GL_Pairs_r3_min_r2';
con_vec(jCon, SV_r3_idx) = 1;
con_vec(jCon, SV_r2_idx) = -1;

% negative contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_SV_GL_Pairs_r2_min_r3';
con_vec(jCon, SV_r3_idx) = -1;
con_vec(jCon, SV_r2_idx) = 1;
    
%% value TE3 - TE1
% positive contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_SV_GL_Pairs_r3_min_r1';
con_vec(jCon, SV_r3_idx) = 1;
con_vec(jCon, SV_r1_idx) = -1;

% negative contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_SV_GL_Pairs_r1_min_r3';
con_vec(jCon, SV_r3_idx) = -1;
con_vec(jCon, SV_r1_idx) = 1;
%% value TE2 - TE1
% positive contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_SV_GL_Pairs_r2_min_r1';
con_vec(jCon, SV_r2_idx) = 1;
con_vec(jCon, SV_r1_idx) = -1;

% negative contrast
jCon = jCon + 1;
con_names{jCon} = 'L_mod_stim_SV_GL_Pairs_r1_min_r2';
con_vec(jCon, SV_r2_idx) = -1;
con_vec(jCon, SV_r1_idx) = 1;
end % function