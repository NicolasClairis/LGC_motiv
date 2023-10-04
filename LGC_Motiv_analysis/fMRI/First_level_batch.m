function[] = First_level_batch(GLM, checking, condition, study_nm, subject_id, NS, biasFieldCorr)
%[] = First_level_batch(GLM, checking, condition, study_nm, subject_id, NS, biasFieldCorr)
% First_level_batch will perform 1st level for LGC motivation fMRI studies.
%
% INPUTS
% GLM: GLM number
%
% checking:
% (0) launch 1st level directly
% (1) display SPM batch before launching to be able to check the inputs
%
% condition:
% 'fMRI': all fMRI compatible data
% 'fMRI_no_move': remove runs with too much movement
%
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
% will be asked if not entered in the inputs
%
% subject_id: list of subject (determined automatically if not defined in
% the inputs)
%
% NS: number of subjects (determined automatically if not defined in
% the inputs)
%
% biasFieldCorr: use bias-field corrected images (1) or not (0)? By default
% will not use bias-field corrected images
%
% See also which_GLM.m, 
% First_level_loadEachCondition.m, First_level_loadRegressors.m,
% First_level_subRunFilter.m, LGCM_contrasts.m and Second_level_batch.m

clc; close all;

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% general parameters

% by default launch the script unless specified otherwise
% if checking = 1, then batch will be displayed before launching
if ~exist('checking','var') || isempty(checking)
    checking = 0;
end
% GLM number
if ~exist('GLM','var') || isempty(GLM) || ~isnumeric(GLM) || GLM < 0
    [GLM_nm] = deal([]);
    while isempty(GLM_nm)
        % repeat until all questions are answered
        info = inputdlg({'GLM number?'});
        [GLM_nm] = info{1};
    end
    GLM = str2double(GLM_nm);
end
if ~exist('study_nm','var') || isempty(study_nm)
    %     study_names = {'study1','study2','fMRI_pilots'};
    %     study_nm_idx = listdlg('ListString',study_names);
    %     study_nm = study_names{study_nm_idx};
    study_nm = 'study1'; % define by default
end
GLMprm = which_GLM(GLM);
add_drv = GLMprm.gal.add_drv;
grey_mask = GLMprm.gal.grey_mask;
autocorrel = GLMprm.gal.autocorrel;

% value of the smoothing during preprocessing?
preproc_sm_kernel = 8;

% use bias-field corrected images or not?
if ~exist('biasFieldCorr','var') || ~ismember(biasFieldCorr,[0,1])
    biasFieldCorr = 0;
end
switch biasFieldCorr
    case 0
        prefix = 'swr';
    case 1
        prefix = 'swbr';
end

% repetition time for fMRI
TR = 2.00;

nb_batch_per_subj = 2; % model + estimate

%% working directories
% computerRoot = LGCM_root_paths();
computerRoot = ['E:',filesep];
% scripts_folder = fullfile(computer_root,'GitHub','LGC_motiv','LGC_Motiv_analysis','fMRI');
% addpath(scripts_folder);
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computerRoot,'fMRI_pilots');
    case 'study1'
        root = fullfile(computerRoot,'study1');
    case 'study2'
        root = fullfile(computerRoot,'study2');
end

%% list subjects to analyze
if ~exist('condition','var') ||...
        ~strcmp(condition(1:4),'fMRI')
    condition = subject_condition;
end
if ~exist('subject_id','var') || ~exist('NS','var') ||...
        isempty(subject_id) || isempty(NS)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
end
%% loop through subjects
matlabbatch = cell(nb_batch_per_subj*NS,1);
for iS = 1:NS
    sub_nm = subject_id{iS};
    % check incompatibility between some GLM parameters and subject
    checkGLM_and_subjectIncompatibility(study_nm, sub_nm, condition, GLMprm);
    
    % define working folders
    subj_folder             = [root, filesep, 'CID',sub_nm];
    subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep];
    subj_anat_folder        = [subj_analysis_folder, filesep, 'anatomical' filesep];
    subj_scans_folder       = [subj_folder, filesep, 'fMRI_scans' filesep];
    subj_behavior_folder    = [subj_folder, filesep, 'behavior' filesep];
    
    % create folder to store the results for the current subject
    switch biasFieldCorr
        case 0
            sm_folderName = [subj_analysis_folder 'functional', filesep,...
                'preproc_sm_',num2str(preproc_sm_kernel),'mm',filesep];
        case 1
            sm_folderName = [subj_analysis_folder 'functional', filesep,...
                'preproc_sm_',num2str(preproc_sm_kernel),'mm_with_BiasFieldCorrection',filesep];
    end
    if ~exist(sm_folderName,'dir')
        mkdir(sm_folderName);
    end
    [resultsFolderName] = fMRI_subFolder(sm_folderName, GLM, condition);
    if ~exist(resultsFolderName,'dir')
        mkdir(resultsFolderName);
    else
        if checking == 0
            rmdir(resultsFolderName,'s');
            mkdir(resultsFolderName);
            warning(['First level folder ',resultsFolderName,' already existed. ',...
                'It was deleted and recreated for CID',sub_nm,'.']);
        elseif checking == 1
            warning(['First level folder ',resultsFolderName,' already exists ',...
                'for CID',sub_nm,'. Please check as much as you want ',...
                'but don''t run this script before deleting the folder.']);
        end
    end
    
    %% define number of runs
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    %% load fMRI data
    subj_scan_folders_names = ls([subj_scans_folder, filesep, '*run*']); % takes all functional runs folders
    % check if extraction worked
    if ~exist('subj_scan_folders_names','var') ||...
            isempty(subj_scan_folders_names)
        error(['The fMRI files could not be extracted. ',...
            'Check the name of the files maybe there is something wrong there.']);
    end
    % clear files that should not be taken into account in the first level
    % analysis
    [subj_scan_folders_names] = First_level_subRunFilter(study_nm, sub_nm,...
        subj_scan_folders_names, [], condition);
    
    %% starting 1st level GLM batch
    sub_idx = nb_batch_per_subj*(iS-1) + 1 ;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.dir            = {resultsFolderName};
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.RT      = TR;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    % loop through runs and tasks
    for iRun = 1:n_runs
        % fix run index if some runs were removed
        [~, jRun] = First_level_subRunFilter(study_nm, sub_nm, [], iRun, condition);
        run_nm = num2str(jRun);
        % erase useless spaces from folder with run name
        n_char = size(subj_scan_folders_names(iRun,:),2);
        for iLetter = 1:n_char
            if strcmp(subj_scan_folders_names(iRun,iLetter),' ') == 0 % erase space
                subj_runFoldername_tmp(iLetter) = subj_scan_folders_names(iRun, iLetter);
            end
        end
        
        % load scans in the GLM
        switch biasFieldCorr
            case 0
                runPath = [subj_scans_folder filesep subj_runFoldername_tmp, filesep,...
                    'preproc_sm_',num2str(preproc_sm_kernel),'mm',filesep]; % run folder
            case 1
                runPath = [subj_scans_folder filesep subj_runFoldername_tmp, filesep,...
                    'preproc_sm_',num2str(preproc_sm_kernel),'mm_with_BiasFieldCorrection',filesep]; % run folder
        end
        preprocessed_filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'.*\.nii$'])); % extracts all the preprocessed swrf files (smoothed, normalized, realigned)
        % in case data is not in .nii but in .img & .hdr
        if isempty(preprocessed_filenames{1})
            preprocessed_filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'.*\.img$'])); % extracts all the preprocessed swrf files (smoothed, normalized, realigned)
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
        % check that task type matches the one predicted by
        % runs_definition.m script to be sure that all works ok
        if (strcmp(task_nm,'physical') && ~strcmp(runs.tasks(iRun),'Ep')) ||...
            (strcmp(task_nm,'mental') && ~strcmp(runs.tasks(iRun),'Em'))
            error(['problem with run task type for subject ',sub_nm,' and run ',num2str(jRun)]);
        end
        % perform 1st level
        matlabbatch = First_level_loadRegressors(matlabbatch, GLMprm, study_nm, sub_nm, sub_idx, iRun, jRun,...
            subj_behavior_folder, currRunBehaviorFileNames, task_nm, computerRoot);
        
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
    
    %% add temporal derivative or not
    if add_drv == 0
        matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    elseif add_drv == 1
        matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
    elseif add_drv == 2
        matlabbatch{sub_idx}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
    end
    %% other parameters
    matlabbatch{sub_idx}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{sub_idx}.spm.stats.fmri_spec.global = 'None';
    switch grey_mask
        case 0 % filter with threshold here if no mask entered
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = 0.8; % default value
        case 4 % filter with threshold here if no mask entered
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = 0.5; % try lower value to include ventral striatum
        case 5 % filter with threshold here if no mask entered
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = 0.3; % try lower value to include ventral striatum
        case 6 % filter with threshold here if no mask entered
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = 0.1; % try lower value to include ventral striatum
        case {1,2,3,7} % no filter here since the mask will do the job
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mthresh = -Inf; % implicitly masks the first level depending on the probability that the voxel is "relevant"
    end
    
    
    %% add grey mask or not
    switch grey_mask
        case {0,4,5,6}
            matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {''};
        case {1,2,3,7}
            maskProbaThreshold_nm = num2str(GLMprm.gal.mask_probaThreshold);
            % find grey matter mask
            switch grey_mask
                case 1 % grey mask per subject
                    error('individual grey matter not extracted yet');
                case 2 % SPM template
                    mask_file_path = [root, filesep, 'grey_matter_mask', filesep,...
                        'bgrey_10perc_SPM_template.nii'];
                    error('copy-paste SPM template file first and then remove this line');
                case 3 %  grey matter filter across subs
                    switch checking
                        case 0
                            mask_file_path = [root, filesep, 'grey_matter_mask', filesep,...
                                'bmean_greyM_',num2str(NS),'_subjects_',condition,'_',maskProbaThreshold_nm,'percentGreyM.nii'];
                        case 1 % to avoid bugs when checking
                            mask_file_path = [root, filesep, 'grey_matter_mask', filesep,...
                                'bmean_greyM_63_subjects_',condition,'_',maskProbaThreshold_nm,'percentGreyM.nii'];
                    end
                case 7 % grey + white matter filter across subs
                    switch checking
                        case 0
                            mask_file_path = [root, filesep, 'greyAndWhite_matter_mask', filesep,...
                                'bmean_greyAndWhiteM_',num2str(NS),'_subjects_',condition,'_',maskProbaThreshold_nm,'percentM.nii'];
                        case 1 % to avoid bugs when checking
                            mask_file_path = [root, filesep, 'greyAndWhite_matter_mask', filesep,...
                                'bmean_greyAndWhiteM_63_subjects_',condition,'_',maskProbaThreshold_nm,'percentM.nii'];
                    end
            end
            % load grey mask (or check if file missing)
            if exist(mask_file_path,'file')
                matlabbatch{sub_idx}.spm.stats.fmri_spec.mask = {mask_file_path};
            else
                error('problem: mean anatomy file not found for this set of subjects');
            end
    end
    
    %% temporal autocorrelation
    switch autocorrel
        case 0
            matlabbatch{sub_idx}.spm.stats.fmri_spec.cvi = 'AR(1)';
        case 1
            matlabbatch{sub_idx}.spm.stats.fmri_spec.cvi = 'FAST';
    end
    
    %% estimate model
    estimate_mdl_rtg_idx = nb_batch_per_subj*iS;
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File',...
        substruct('.','val', '{}',{sub_idx}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.write_residuals = 0;
    % classical model
    matlabbatch{estimate_mdl_rtg_idx}.spm.stats.fmri_est.method.Classical = 1;
end

%% display spm batch before running it or run it directly
if exist('checking','var') && checking == 1
    spm_jobman('interactive',matlabbatch);
    %     spm_jobman('run',matlabbatch);
elseif ~exist('checking','var') || isempty(checking) || checking == 0
    spm_jobman('run',matlabbatch);
end


end % function