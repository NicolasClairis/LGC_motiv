
clc; close all;

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% general parameters
% check the batch before launching the script?
checking = 0;
% GLM number
[GLM_nm] = deal([]);
while isempty(GLM_nm)
    % repeat until all questions are answered
    info = inputdlg({'GLM number?'});
    [GLM_nm] = info{1};
end
GLM = str2double(GLM_nm);
GLMprm = which_GLM(GLM);
add_drv = GLMprm.gal.add_drv;
grey_mask = GLMprm.gal.grey_mask;

% repetition time for fMRI
TR = 2.00;

nb_batch_per_subj = 2; % model + estimate

%% working directories
computer_root = fullfile('C:','Users','clairis','Desktop');
scripts_folder = fullfile(computer_root,'GitHub','LGC_motiv','LGC_Motiv_analysis','fMRI');
addpath(scripts_folder);
root = fullfile(computer_root,'study1','fMRI_pilots');

%% list subjects to analyze
subject_id = {'pilot_s2'};%'pilot_s1','pilot_s2'
NS = length(subject_id);
%% loop through subjects
matlabbatch = cell(nb_batch_per_subj*NS,1);
for iS = 1:NS
    sub_nm = subject_id{iS};
    % define working folders
    subj_folder             = [root, filesep, sub_nm];
    subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep];
    subj_anat_folder        = [subj_analysis_folder, filesep, 'anatomical' filesep];
    subj_scans_folder       = [subj_folder, filesep, 'fMRI_scans' filesep];
    subj_behavior_folder    = [subj_folder, filesep, 'behavior' filesep];
    
    % create folder to store the results for the current subject
    resultsFolderName = [subj_analysis_folder 'functional', filesep,...
        'GLM',num2str(GLM)];
    mkdir(resultsFolderName);
    
    %% define number of runs
    switch sub_nm
        case 'pilot_s1'
            nb_runs = 2;
        case 'pilot_s2'
            nb_runs = 1;
        otherwise
            nb_runs = 4;
    end
    
    %% load fMRI data
    subj_scan_folders_names = ls([subj_scans_folder, filesep, '*_run*']); % takes all functional runs folders
    % remove AP/PA corrective runs
    [subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names);
    
    %% starting 1st level GLM batch
    sub_idx = nb_batch_per_subj*(iS-1) + 1 ;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.dir            = {resultsFolderName};
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.RT      = TR;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    % loop through runs and tasks
    for iRun = 1:nb_runs
        run_nm = num2str(iRun);
        % erase useless spaces from folder with run name
        n_char = size(subj_scan_folders_names(iRun,:),2);
        for iLetter = 1:n_char
            if strcmp(subj_scan_folders_names(iRun,iLetter),' ') == 0 % erase space
                subj_runFoldername_tmp(iLetter) = subj_scan_folders_names(iRun, iLetter);
            end
        end
        
        % load scans in the GLM
        cd([subj_scans_folder filesep subj_runFoldername_tmp, filesep]); % go to run folder
        preprocessed_filenames = cellstr(spm_select('ExtFPList',pwd,'^swr.*\.nii$')); % extracts all the preprocessed swrf files (smoothed, normalized, realigned)
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).scans = preprocessed_filenames;
        
        %% load regressors of interest
        % 1) identify which task corresponds to the current run
        currRunBehaviorFileNames = ls([subj_behavior_folder,'*_session',num2str(iRun),'_*_task.mat']);
        if size(currRunBehaviorFileNames,1) > 1
            error(['problem file identification: too many files popping out with run number',run_nm]);
        end
        if strcmp(currRunBehaviorFileNames(16:23),'physical') ||...
                strcmp(currRunBehaviorFileNames(17:24),'physical')
            task_nm = 'physical';
        elseif strcmp(currRunBehaviorFileNames(16:21),'mental') ||...
                strcmp(currRunBehaviorFileNames(17:22),'mental')
            task_nm = 'mental';
        else
            error('problem in identifying task type because file name doesn''t match');
        end
        % perform 1st level
        matlabbatch = First_level_loadRegressors(matlabbatch, GLMprm, sub_idx, iRun,...
            subj_behavior_folder, currRunBehaviorFileNames, task_nm);
        
        %% global run parameters (rp movement file, etc.)
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).multi = {''};
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).regress = struct('name', {}, 'val', {});
        mvmtFolder = [subj_scans_folder, filesep, subj_runFoldername_tmp, filesep];
        movement_file = ls([mvmtFolder, 'rp*']);
        movement_filePath = [mvmtFolder, movement_file];
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).multi_reg = {movement_filePath};
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).hpf = 128;
        
        %% clear run name to avoid bugs across runs
        clear('subj_runFoldername_tmp');
    end % run loop
    
    %% global parameters for subject batch
    matlabbatch{sub_idx}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    
    % add temporal derivative or not
    if add_drv == 0
        matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    elseif add_drv == 1
        matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
    elseif add_drv == 2
        matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
    end
    
    matlabbatch{sub_idx}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.global = 'None';
    switch grey_mask
        case 0 % filter with threshold here if no mask entered
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = 0.8; % default value
        case {1,2,3} % no filter here since the mask will do the job
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = -Inf; % implicitly masks the first level depending on the probability that the voxel is "relevant"
    end
    
    
    % add grey mask or not
    switch grey_mask
        case 0
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {''};
        case {1,2,3}
            % find grey matter mask
            switch grey_mask
                case 1 % grey mask per subject
                    mask_file = ls([subj_anat_folder filesep 'bmwc1s*']); % modulated grey matter mask
                    mask_file_path = [subj_anat_folder, mask_file];
                case 2 % grey matter filter across subs
                    %                     mask_file_path = [root, filesep, 'Second_level', filesep,...
                    %                         'mean_anatomy', filesep,'mean_greyMatter_XX_subjects.nii'];
                    error('not ready yet: you need to average the grey matter of your participants first');
                case 3 % SPM template
                    mask_file_path = [root, filesep, 'Second_level', filesep,...
                        'mean_anatomy', filesep,'bgrey_10perc_SPM_template.nii'];
            end
            % load grey mask (or check if file missing)
            if exist(mask_file_path,'file')
                matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {mask_file_path};
            else
                error('problem: mean anatomy file not found for this number of subjects');
            end
    end
    
    matlabbatch{sub_idx}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %% estimate model
    estimate_mdl_rtg_idx = nb_batch_per_subj*iS;
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{sub_idx}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.write_residuals = 0;
    % classical model
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Classical = 1;
end

%% display spm batch before running it
if exist('checking','var') && checking == 1
    spm_jobman('interactive',matlabbatch);
    %     spm_jobman('run',matlabbatch);
elseif ~exist('checking','var') || isempty(checking) || checking == 0
    spm_jobman('run',matlabbatch);
end