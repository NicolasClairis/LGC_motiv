function [list_all_GLMprm, n_potReg] = LGC_list_potential_onsets_and_modulators
%[list_all_GLMprm, n_potReg] = LGC_list_potential_onsets_and_modulators
% MS2_list_potential_onsets_and_modulators outputs all the potential onsets
% and regressors that can be used in the first level GLM.
%
% See also which_GLM_MS2, First_level_MS2_megaconcatenation_NicoC_batch and contrasts_megaconcatenation_MS2_NicoC

%% general parameters
list_all_GLMprm.gal = {'grey_mask','add_drv','orth_vars',...
    'onsets_only',...
    'FIR','FIR_dur','FIR_nBins',...
    'zscore'};
n_potReg.gal = length(list_all_GLMprm.gal);

%% Reinforcement-learning
list_all_GLMprm.RL.onsets       = {'o_stim',    'o_answer',     'o_chosen',     'o_fbk',    'o_cross',  'o_missed_trials_stim'};
list_all_GLMprm.RL.durations    = {'dur_stim',  'dur_answer',   'dur_chosen',   'dur_fbk',  'dur_cross','dur_missed_trials_stim'};
n_potReg.RL.onsets = length(list_all_GLMprm.RL.onsets);
list_all_GLMprm.RL.modulators = {'mod_stim','mod_chosen','mod_fbk'};
n_potReg.RL.mods = length(list_all_GLMprm.RL.modulators);
% stimuli on screen modulators
list_all_GLMprm.RL.mod_stim = {'trialN',...
    'mdl_type','mdl_n','SV','dQ','pBest',...
    'ROI_activity_yn','ROI_activity_GLM','ROI_activity_period','ROI_activity_ROI_nm',...
    'RT'};
n_potReg.RL.mod_stim = length(list_all_GLMprm.RL.mod_stim);
% chosen option circled in red modulators
list_all_GLMprm.RL.mod_chosen = {'trialN',...
    'mdl_type','mdl_n','SV','dQ','pBest',...
    'ROI_activity_yn','ROI_activity_GLM','ROI_activity_period','ROI_activity_ROI_nm',...
    'RT'};
n_potReg.RL.mod_chosen = length(list_all_GLMprm.RL.mod_chosen);
% feedback modulators
list_all_GLMprm.RL.mod_fbk = {'trialN','fbk',...
    'mdl_type','mdl_n','PE','PE_bis',...
    'totalGain',...
    'ROI_activity_yn','ROI_activity_GLM','ROI_activity_period','ROI_activity_ROI_nm'};
n_potReg.RL.mod_fbk = length(list_all_GLMprm.RL.mod_fbk);


end % function