function[] = LGC_First_level_RL_NicoC_batch(GLM, checking, subject_id)
%[] = LGC_First_level_RL_NicoC_batch(GLM, checking, subject_id, pc_cluster)
% Batch script for first level fMRI for multiband (1.1s) fMRI sequence of
% MotiScan-2
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
pc_clusterPath = 'pc_home';
% pc_clusterPath = 'pc_lab';
pc_cluster = 'pc';
switch pc_clusterPath
    case 'pc_home'
        root = ['F:', filesep, 'study2', filesep,'pilots', filesep];
        scripts_folder = fullfile('C:','Users','Nicolas Clairis','Documents','GitHub',...
            'Bat-Motiv','Bat-analyse','various',...
            'NicoC_analysis_fMRI','one_seq_batchs');
    case 'pc_lab'
        root = ['E:', filesep, 'study2', filesep,'pilots', filesep];
        scripts_folder = fullfile('C:','Users','clairis','Desktop','GitHub',...
            'Bat-Motiv','Bat-analyse','various',...
            'NicoC_analysis_fMRI','one_seq_batchs');
end
warning('check paths here');

%% by default checking = 0 if not selected
if ~exist('checking','var') || isempty(checking)
    checking = 0;
end

%% GLM
if ~exist('GLM','var') || isempty(GLM)
    GLM = input(sprintf('GLM number? \n '));
end
% load specific GLM parameters
[GLMprm] = which_GLM_LGC_pilot(GLM);

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% main parameters
learningRuns    = 3;
nbRuns = learningRuns;

% particular sequence parameters
TR = 2.0;
warning('PLEASE CHECK TR is ok?');

nb_batch_per_subj = 2; % (model + estimate)

matlabbatch = cell(nb_batch_per_subj*NS,1);

for iSub = 1:NS
    
    sub_nm = subject_id{iSub};
    % extract subject number
    if strcmp(sub_nm(3),'_')
        subid   = sub_nm(2);
    elseif strcmp(sub_nm(3),'_') == 0 && strcmp(sub_nm(4),'_')
        subid   = sub_nm(2:3);
    end
    % define working folders
    subj_folder             = [root, filesep, sub_nm];
    subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep];
    subj_scans_folder       = [subj_folder, filesep, 'fMRI_scans' filesep];
    subj_behavior_folder    = [subj_folder, filesep, 'behavior' filesep];
    % folder names
    warning('DEFINE RUN NAME HERE');
    subj_scan_folders_names = ls([subj_scans_folder, filesep, '*TR1100_2_5iso_PA_RUN*']); % takes all functional runs folders (if TR = 1.10s, for multiband seq in particular)
    
    % create folder for storing data for this subject
    filename = [subj_analysis_folder 'functional', filesep,...
        'GLM',num2str(GLM),'_megaconcatenation'];
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
        runname = num2str(iRun);
        % erase useless spaces from folder with run name
        n_char = size(subj_scan_folders_names(iRun,:),2);
        for iLetter = 1:n_char
            if strcmp(subj_scan_folders_names(iRun,iLetter),' ') == 0 % erase space
                subj_runFoldername(iLetter) = subj_scan_folders_names(iRun,iLetter);
            end
        end
        
        % load scans in the GLM
        cd([subj_scans_folder filesep subj_runFoldername, filesep]); % go to run folder
        preprocessed_filenames = cellstr(spm_select('ExtFPList',pwd,'^swrf.*\.img$')); % extracts all the preprocessed swrf files (smoothed, normalized, realigned)
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).scans = preprocessed_filenames;
        
        % identify task corresponding to this run
        behaviorFile = ls([subj_behavior_folder, filesep,...
            'global_sub_',subid,'_session_',runname,'_*']);
        matlabbatch = LGC_First_level_RL_prm(matlabbatch, GLMprm,...
            subid, sub_idx, iRun, subj_analysis_folder);

        %% global run parameters (rp movement file, etc.)
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).multi = {''};
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).regress = struct('name', {}, 'val', {});
        mvmtFolder = [subj_scans_folder, filesep, subj_runFoldername, filesep];
        switch pc_cluster
            case 'pc'
                movement_file = ls([mvmtFolder, 'rp*']);
            case 'cluster'
                movement_file = cell2mat(get_folder_list('rp*',mvmtFolder));
        end
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
    matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = 0.8; % default value
    matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {''};    
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