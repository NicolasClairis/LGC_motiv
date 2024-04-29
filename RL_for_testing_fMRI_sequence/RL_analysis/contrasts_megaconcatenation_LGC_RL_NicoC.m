function [] = contrasts_megaconcatenation_LGC_RL_NicoC(GLM, checking, subject_id)
%% contrasts_megaconcatenation_LGC_RL_NicoC(GLM, checking, subject_id, pc_cluster) 
% concatenates together all tasks in the same GLM
%
% INPUTS
% GLM: GLM number to be used (check which_GLM_MS2.m for more details about each
% GLM specifications)
%
% checking: if 0 or empty, run the GLM with all the subjects, if equal to
% 1, displays the interactive window on only one subject
%
% subject_id: input string subject identification name (but can be left empty)
%
% pc_cluster: 'pc' or 'cluster' to know in which folder you should find the
% files of interest
%
% See also which_GLM_MS2.m, First_level_MS2_megaconcatenation_NicoC_batch

close all; clc;

%% working directories
root = ['E:', filesep, 'study2', filesep,'pilots', filesep,...
    'fMRI_pilots',filesep];
scripts_folder = fullfile('C:','Users','clairis','Desktop','GitHub',...
    'LGC_motiv','RL_for_testing_fMRI_sequence','RL_analysis');

%% by default checking = 0 if not selected
if ~exist('checking','var') || isempty(checking)
    checking = 0;
end

%% GLM parameters
if ~exist('GLM','var') || isempty(GLM)
    GLM = 2;
end
[GLMprm] = which_GLM_LGC_pilot(GLM);
grey_mask = GLMprm.gal.grey_mask;
%% subjects identification
if ~exist('subject_id','var') || isempty(subject_id)
    subject_id = {'fMRI_pilot1_AC'};
end
NS = length(subject_id);

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% preprocessing to use
biasFieldCorr = 1;

switch biasFieldCorr
    case 0
        preproc_folder = 'preproc_sm_8mm';
    case 1
        preproc_folder = 'preproc_sm_8mm_with_BiasFieldCorrection';
end
switch grey_mask
    case 0
        maskThresh = 0.8; % 0.8 by default
    case 2
        maskThresh = 0.1; % lower
end

%% loop through subjects to extract all the regressors
matlabbatch = cell(NS,1);
batch_idx = 0;
for iSub = 1:NS
    
    %% extract current subject name
    sub_nm = subject_id{iSub};
    
    batch_idx = batch_idx + 1;
    
    switch biasFieldCorr
        case 0
            switch grey_mask
                case {0,2}
                    run_foldername = ['GLM',num2str(GLM),'_SPMmask',num2str(maskThresh*100),'percent'];
                case 1
                    run_foldername = ['GLM',num2str(GLM),'_individualGreyMask'];
            end
        case 1
            switch grey_mask
                case {0,2}
                    run_foldername = ['GLM',num2str(GLM),'_SPMmask',num2str(maskThresh*100),'percent_with_BiasFieldCorrection'];
                case 1
                    run_foldername = ['GLM',num2str(GLM),'_individualGreyMask_with_BiasFieldCorrection'];
            end
    end
    matlabbatch{batch_idx}.spm.stats.con.spmmat = {fullfile(root,sub_nm,'fMRI_analysis','functional',...
        preproc_folder,filesep,run_foldername,'SPM.mat')};
    
    %% extract contrasts list (vectors + corresponding names
    [con_names, con_vec] = LGC_RL_load_con(GLMprm);
    n_con = length(con_names);
    
    %% add each contrast to the list
    for iCon = 1:n_con
        matlabbatch{batch_idx}.spm.stats.con.consess{iCon}.tcon.name     = con_names{iCon};
        matlabbatch{batch_idx}.spm.stats.con.consess{iCon}.tcon.weights  = con_vec(iCon,:);
        matlabbatch{batch_idx}.spm.stats.con.consess{iCon}.tcon.sessrep  = 'none';
    end
    
    matlabbatch{batch_idx}.spm.stats.con.delete = 1; % deletes previous contrasts
    
end % subject loop

cd(scripts_folder);

%% display spm batch before running it
%% display spm batch before running it or run it directly
switch checking
    case 0
        spm_jobman('run',matlabbatch);
    case 1
        spm_jobman('interactive',matlabbatch);
end

end % function end