function[] = LGC_First_level_RL_NicoC_batch(GLM, checking)%, subject_id)
%[] = LGC_First_level_RL_NicoC_batch(GLM, checking, subject_id)
% Batch script for first level fMRI for fMRI sequence of LGC
% Use of SPM-12
% define subjects manually, ask the GLM you want to use and then displays
% it into SPM GUI ("interactive") or runs it directly ("run") depending on
% the parameter entered in the last line
%
% INPUTS
% GLM: GLM identification number (see which_GLM_MS2.m for details)
%
% checking
% (0) directly runs the first level GLM
% (1) only displays interactive window on one subject to see if everything
% is fine
%
% subject_id: input string subject identification name (but can be left empty)
%
% Makes a global GLM where all three tasks of MS2 (learning, grip and stroop) are pooled together
%
% See also: preprocessing_NicoC_batch, onsets_for_fMRI_MS2,
% group_onsets_for_fMRI_MS2, which_GLM_MS2

%% working directories
root = ['E:', filesep, 'study2', filesep,'pilots', filesep,...
    'fMRI_pilots',filesep];
scripts_folder = fullfile('C:','Users','clairis','Desktop','GitHub',...
    'LGC_motiv','RL_for_testing_fMRI_sequence','RL_analysis');

%% by default checking = 0 if not selected
if ~exist('checking','var') || isempty(checking)
    checking = 0;
end

%% GLM
if ~exist('GLM','var') || isempty(GLM)
    %     GLM = input(sprintf('GLM number? \n '));
    GLM=2;
end
% load specific GLM parameters
[GLMprm] = which_GLM_LGC_pilot(GLM);
grey_mask = GLMprm.gal.grey_mask;
%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% main parameters
% number of subjects
subject_id = {'fMRI_pilot1_AC'};
NS = length(subject_id);

% bias-field corrected preprocessing or not?
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
        maskThresh = 0.1; % lower threshold
    otherwise
        maskThresh = -Inf;
end
learningRuns    = 3;
nbRuns = learningRuns;

% particular sequence parameters
TR = 2.0;

nb_batch_per_subj = 2; % (model + estimate)

matlabbatch = cell(nb_batch_per_subj*NS,1);

for iSub = 1:NS
    
    sub_nm = subject_id{iSub};
    % extract subject number
    if strcmp(sub_nm,'fMRI_pilot1_AC')
        subid = '1';
    else
        error('case where sub_nm different from fMRI_pilot1_AC not ready yet');
    end
    % define working folders
    subj_folder             = [root, filesep, sub_nm];
    subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep];
    subj_anat_folder = [subj_analysis_folder,filesep,'anatomical',filesep];
    subj_scans_folder       = [subj_folder, filesep, 'fMRI_scans' filesep];
    % folder names
    subj_scan_folders_names = ls([subj_scans_folder, filesep, '*_run*']); % takes all functional runs folders
    
    % create folder for storing data for this subject
    switch biasFieldCorr
        case 0
            switch grey_mask
                case {0,2}
                    filename = [subj_analysis_folder 'functional', filesep,...
                        preproc_folder,filesep,...
                        'GLM',num2str(GLM),'_SPMmask',num2str(maskThresh*100),'percent'];
                case 1
                    filename = [subj_analysis_folder 'functional', filesep,...
                        preproc_folder,filesep,...
                        'GLM',num2str(GLM),'_individualGreyMask'];
            end
        case 1 % bias-field corrected files
            switch grey_mask
                case {0,2}
                    filename = [subj_analysis_folder 'functional', filesep,...
                        preproc_folder,filesep,...
                        'GLM',num2str(GLM),'_SPMmask',num2str(maskThresh*100),'percent_with_BiasFieldCorrection'];
                case 1
                    filename = [subj_analysis_folder 'functional', filesep,...
                        preproc_folder,filesep,...
                        'GLM',num2str(GLM),'_individualGreyMask_with_BiasFieldCorrection'];
            end
    end
    mkdir(filename);
    
    %% starting 1st level GLM batch
    sub_idx = nb_batch_per_subj*(iSub-1) + 1 ;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.dir            = {filename};
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.RT      = TR;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    
    % loop through runs and tasks
    for iRun = 1:nbRuns
%         runname = num2str(iRun);
        % erase useless spaces from folder with run name
        n_char = size(subj_scan_folders_names(iRun,:),2);
        for iLetter = 1:n_char
            if strcmp(subj_scan_folders_names(iRun,iLetter),' ') == 0 % erase space
                subj_runFoldername(iLetter) = subj_scan_folders_names(iRun,iLetter);
            end
        end
        
        % load scans in the GLM
        cd([subj_scans_folder filesep subj_runFoldername, filesep, preproc_folder,filesep]); % go to run folder
        switch biasFieldCorr
            case 0
                preprocessed_filenames = cellstr(spm_select('ExtFPList',pwd,'^swr.*\.nii$')); % extracts all the preprocessed swrf files (smoothed, normalized, realigned)
            case 1
                preprocessed_filenames = cellstr(spm_select('ExtFPList',pwd,'^swbr.*\.nii$')); % extracts all the preprocessed swrf files (smoothed, normalized, realigned)
        end
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).scans = preprocessed_filenames;
        
        % identify task corresponding to this run
        cd(scripts_folder);
        matlabbatch = LGC_First_level_RL_prm(matlabbatch, GLMprm,...
            subid, sub_idx, iRun, subj_analysis_folder);

        %% global run parameters (rp movement file, etc.)
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).multi = {''};
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).regress = struct('name', {}, 'val', {});
        mvmtFolder = [subj_scans_folder, filesep, subj_runFoldername, filesep];
        movement_file = ls([mvmtFolder, 'rp*']);
        movement_filePath = [mvmtFolder, movement_file];
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).multi_reg = {movement_filePath};
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).hpf = 128;
        
        % clear folder variable at the end of each run loop to avoid bugs
        clear('subj_runFoldername');
    end % run loop
    
    
    %% global parameters for subject batch
    matlabbatch{sub_idx}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % no derivative
    matlabbatch{sub_idx}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.global = 'None';
    % no grey mask
    matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = maskThresh; % default value
    switch grey_mask
        case {0,2}
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {''};
        case 1
            % find grey matter mask
            mask = ls([subj_anat_folder,'bmwc1*']); % modulated grey matter mask
            mask = [subj_anat_folder,mask];
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {mask};
    end
    matlabbatch{sub_idx}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    
    %% estimate model
    estimate_mdl_rtg_idx = nb_batch_per_subj*iSub;
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{sub_idx}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.write_residuals = 0;
    % classical model
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Classical = 1;
end % subject loop

cd(scripts_folder)

%% display spm batch before running it or run it directly
switch checking
    case 0
%         spm_jobman('interactive',matlabbatch);
        spm_jobman('run',matlabbatch);
    case 1
%         spm_jobman('run',matlabbatch);
        spm_jobman('interactive',matlabbatch);
end

end % function