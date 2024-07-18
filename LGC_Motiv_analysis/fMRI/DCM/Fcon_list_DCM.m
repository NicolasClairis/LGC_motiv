function[con_names, con_vector] = Fcon_list_DCM(study_nm, sub_nm, GLM, DCM_mode, computer_root, preproc_sm_kernel, condition, biasFieldCorr)
% [con_names, con_vector] = Fcon_list_DCM(study_nm, sub_nm, GLM, DCM_mode, computer_root, preproc_sm_kernel, condition, biasFieldCorr)
% [con_names, con_vector] = LGCM_contrasts(study_nm, sub_nm, GLM,...
%   computer_root, preproc_sm_kernel, condition, biasFieldCorr, DCM_GLM)
% LGCM_contrasts will define the contrast names and contrast vector for the
% subject and study entered in input.
%
% INPUTS
% computer_root: path where data is stored
%
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
%
% sub_nm: subject name
%
% GLM: GLM number
%
% DCM_mode:
% (1) all sessions modeled independently like in a classic univariate GLM
% => hard to manipulate for DCM but could be useful for testing
% session-specific effects or comparing sessions
% (2) sessions pooled within each task (ex: session 1 and 3 of physical
% effort will be concatenated into one single regressor) but each task will
% be modeled separately
% (3) all sessions pooled together
% (4) all trial periods are pooled together across sessions except for
% choice and effort which are modeled independently for each task (but
% pooled across sessions of the same task)
% (5) all trial periods are pooled together across sessions except for
% the effort period which is modeled independently for each task (but
% pooled across sessions of the same task)
%
% preproc_sm_kernel: kernel used in preprocessing for smoothing the data
%
% condition: define subjects and runs to include
% 'fMRI': all subjects where fMRI ok
% 'fMRI_no_move': remove runs with too much movement
%
% biasFieldCorr: binary variable indicating whether results are based on
% bias-field corrected files (1) or not (0)
%
% DCM_GLM: use preprocessing adapted for DCM or not (the main difference is
% that the preprocessing will be including or not a slice-timing correction)
% (0) use basic preprocessing
% (1) use preprocessing including slice-timing correction
%
% OUTPUTS
% con_names: list of contrast names
%
% con_vector: matrix with all contrasts organized in nb_contrasts (rows)*nb_regressors (columns)

%% initialize the variables of interest
[reg_names.Ep, reg_names.Em] = deal( {} );
n_regs.Ep = 0;
n_regs.Em = 0;

end % function

function[] = dispRegFn(text, dispRegs)
% function to display the namer of the contrasts
switch dispRegs
    case 1
        disp(text);
end
end