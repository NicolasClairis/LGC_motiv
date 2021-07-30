function[] = LGCM_contrasts_spm(study_nm, GLM, checking)
% [] = LGCM_contrasts_spm(study_nm, GLM, checking)
% LGCM_contrasts_spm prepares the batch for SPM to perform the first level
% contrasts in each individual of the study
%
% INPUTS
% study_nm: study name
%
% GLM: GLM number
%
% checking: display batch before performing it or not? (1 by default)


close all; clc;

%% working directories
computer_root = fullfile('C:','Users','clairis','Desktop');
% computer_root = fullfile('C:','Users','Loco','Downloads');
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        root = fullfile(computer_root,'study1');
    case 'study2_clinical'
        root = fullfile(computer_root,'study2');
end

%% checking by default the batch before launching it
if ~exist('checking','var') || isempty(checking) || ~ismember(checking,[0,1])
    checking = 1;
end

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% define subjects of interest
[subject_id, NS] = LGCM_subject_selection(study_nm);

%% loop through subjects to extract all the regressors
matlabbatch = cell(NS,1);

for iSubject = 1:NS
    %% extract current subject name
    sub_nm = subject_id{iSubject};
    
    %% extract contrasts list (vectors + corresponding names
    [con_names, con_vector] = LGCM_contrasts(study_nm, sub_nm, GLM);
    
    %% define results directory
    matlabbatch{iSubject}.spm.stats.con.spmmat = {fullfile(root,sub_nm,'fMRI_analysis','functional',['GLM',num2str(GLM)],'SPM.mat')};
    %% add each contrast to the list
    for iCon = 1:length(con_names)
        matlabbatch{iSubject}.spm.stats.con.consess{iCon}.tcon.name     = con_names{iCon};
        matlabbatch{iSubject}.spm.stats.con.consess{iCon}.tcon.weights  = con_vector(iCon,:);
        matlabbatch{iSubject}.spm.stats.con.consess{iCon}.tcon.sessrep  = 'none';
    end
    
    matlabbatch{iSubject}.spm.stats.con.delete = 1; % deletes previous contrasts
    
end % subject loop

%% display spm batch before running it
if exist('checking','var') && checking == 1
    spm_jobman('interactive',matlabbatch);
%     spm_jobman('run',matlabbatch);
elseif ~exist('checking','var') || isempty(checking) || checking == 0
    spm_jobman('run',matlabbatch);
end

end % function