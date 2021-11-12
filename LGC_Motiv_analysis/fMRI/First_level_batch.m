function[] = First_level_batch(study_nm)
% First_level_batch will perform 1st level for LGC motivation fMRI studies.
%
% INPUTS
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
%
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
if ~exist('study_nm','var') || isempty(study_nm)
    study_names = {'fMRI_pilots','study1','study2'};
    study_nm_idx = listdlg('ListString',study_names);
    study_nm = study_names{study_nm_idx};
end
GLM = str2double(GLM_nm);
GLMprm = which_GLM(GLM);
add_drv = GLMprm.gal.add_drv;
grey_mask = GLMprm.gal.grey_mask;

% repetition time for fMRI
TR = 2.00;

nb_batch_per_subj = 2; % model + estimate

%% working directories
computer_root = LGCM_root_paths();
% scripts_folder = fullfile(computer_root,'GitHub','LGC_motiv','LGC_Motiv_analysis','fMRI');
% addpath(scripts_folder);
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        root = fullfile(computer_root,'study1');
    case 'study2'
        root = fullfile(computer_root,'study2');
end

%% list subjects to analyze
[subject_id, NS] = LGCM_subject_selection(study_nm);
%% loop through subjects
matlabbatch = cell(nb_batch_per_subj*NS,1);
for iS = 1:NS
    sub_nm = subject_id{iS};
    % check incompatibility between some GLM and some subjects
    checkGLM_and_subjectIncompatibility(study_nm, sub_nm, GLMprm);
    
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
    switch study_nm
        case 'fMRI_pilots'
            switch sub_nm
                case {'pilot_s1','pilot_s3'}
                    nb_runs = 2;
                case 'pilot_s2'
                    nb_runs = 1;
                otherwise
                    nb_runs = 4;
            end
        case 'study1'
            switch sub_nm
                case 'CID074' % CID074: ignore first run which crashed at the beginning due to high voltage
                    nb_runs = 3;
                otherwise
                    nb_runs = 4;
            end
    end
    
    %% load fMRI data
    subj_scan_folders_names = ls([subj_scans_folder, filesep, '*run*']); % takes all functional runs folders
    % clear files that should not be taken into account in the first level
    % analysis
    [subj_scan_folders_names] = First_level_subRunFilter(study_nm, sub_nm, subj_scan_folders_names);
    
    %% starting 1st level GLM batch
    sub_idx = nb_batch_per_subj*(iS-1) + 1 ;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.dir            = {resultsFolderName};
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.RT      = TR;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    % loop through runs and tasks
    for iRun = 1:nb_runs
        % fix run index if some runs were removed
        [~, jRun] = First_level_subRunFilter(study_nm, sub_nm, [], iRun);
        run_nm = num2str(jRun);
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
        % in case data is not in .nii but in .img & .hdr
        if isempty(preprocessed_filenames{1})
            preprocessed_filenames = cellstr(spm_select('ExtFPList',pwd,'^swr.*\.img$')); % extracts all the preprocessed swrf files (smoothed, normalized, realigned)
        end
        % check if still empty => if yes, stop the script
        if isempty(preprocessed_filenames{1})
            error('problem with format of preprocessed files: impossible to find them. Did you preprocess the data?');
        end
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).scans = preprocessed_filenames;
        
        %% load regressors of interest
        % 1) identify which task corresponds to the current run
        currRunBehaviorFileNames = ls([subj_behavior_folder,'*_session',run_nm,'_*_task.mat']);
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
        matlabbatch = First_level_loadRegressors(matlabbatch, GLMprm, study_nm, sub_idx, iRun,...
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


end % function