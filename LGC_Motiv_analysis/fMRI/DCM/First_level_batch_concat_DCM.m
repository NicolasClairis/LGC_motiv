function[] = First_level_batch_concat_DCM(GLM, checking, condition, study_nm, subject_id, NS, biasFieldCorr)
%[] = First_level_batch_concat_DCM(GLM, checking, condition, study_nm, subject_id, NS, biasFieldCorr)
% First_level_batch_concat_DCM will perform 1st level for LGC motivation fMRI studies.
% The main difference with Firt_level_batch.m is that
% First_level_batch_concat_DCM will concatenate all the sessions together
% in one single session in order to have only one single session to use for
% DCM. You then need to call the spm_fmri_concatenate function to account
% for the session-specific effects and, in principle, you are then good to
% go with the DCM analysis.
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
% will not use bias-field corrected images (=0)
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

% determine how sessions will be concatenated
[DCM_mode] = which_DCM_mode_for_GLM;

% value of the smoothing during preprocessing?
preproc_sm_kernel = 8;

% use bias-field corrected images or not?
if ~exist('biasFieldCorr','var') || ~ismember(biasFieldCorr,[0,1])
    biasFieldCorr = 0;
end

% adapt files if bias-field corrected preprocessing or not
switch biasFieldCorr
    case 0
        prefix = 'swar'; % smoothed(s)-normalized(w)-slice-timed(a)-realigned(r) files
    case 1
        error('case with bias-field correction and DCM slice-timing preprocessing not ready yet');
end

% repetition time for fMRI
TR = 2.00;

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
        isempty(condition) ||...
        ~strcmp(condition(1:4),'fMRI')
    condition = subject_condition;
end
if ~exist('subject_id','var') || ~exist('NS','var') ||...
        isempty(subject_id) || isempty(NS)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
end

%% quick check condition and GLM are compatible
isGLMok_bis(GLMprm, condition);

%% loop through subjects
matlabbatch1 = cell(NS,1);
for iS = 1:NS
    sub_nm = subject_id{iS};
    % check incompatibility between some GLM parameters and subject
    checkGLM_and_subjectIncompatibility(study_nm, sub_nm, condition, GLMprm);
    
    % define working folders
    subj_folder             = [root, filesep, 'CID',sub_nm];
    subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep];
%     subj_anat_folder        = [subj_analysis_folder, filesep, 'anatomical' filesep];
    subj_scans_folder       = [subj_folder, filesep, 'fMRI_scans' filesep];
    subj_behavior_folder    = [subj_folder, filesep, 'behavior' filesep];
    
    % create folder to store the results for the current subject
    switch biasFieldCorr
        case 0
            sm_folderName = [subj_analysis_folder 'functional', filesep,...
                'preproc_sm_',num2str(preproc_sm_kernel),'mm_DCM',filesep];
        case 1
            error('case with bias-field correction and DCM slice-timing preprocessing not ready yet');
    end
    if ~exist(sm_folderName,'dir')
        mkdir(sm_folderName);
    end
    [resultsFolderName] = fMRI_subFolder_DCM(sm_folderName, GLM, condition, DCM_mode);
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
    
    %% initialize number of scans for each run (important for the concatenation required for DCM)
    n_scansPerRun.(['CID',sub_nm]) = NaN(1,n_runs);
    
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
    matlabbatch1{iS}.spm.stats.fmri_spec.dir            = {resultsFolderName};
    matlabbatch1{iS}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch1{iS}.spm.stats.fmri_spec.timing.RT      = TR;
    matlabbatch1{iS}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch1{iS}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    %% load scans in the batch (concatenate all sessions together)
    % loop through runs and tasks
    for iRun = 1:n_runs
        % erase useless spaces from folder with run name
        subj_runFoldername_tmp = strrep(subj_scan_folders_names(iRun, :),' ','');
        
        % load scans in the GLM
        switch biasFieldCorr
            case 0
                runPath = [subj_scans_folder filesep subj_runFoldername_tmp, filesep,...
                    'preproc_sm_',num2str(preproc_sm_kernel),'mm_DCM',filesep]; % run folder
            case 1
                error('case with bias-field correction and DCM slice-timing preprocessing not ready yet');
        end
        preprocessed_filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'.*\.nii$'])); % extracts all the preprocessed swar files (smoothed, normalized, slice-timed, realigned)
        % in case data is not in .nii but in .img & .hdr
        if isempty(preprocessed_filenames{1})
            preprocessed_filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'.*\.img$'])); % extracts all the preprocessed swar files (smoothed, normalized, slice-timed, realigned)
        end
        % check if still empty => if yes, stop the script
        if isempty(preprocessed_filenames{1})
            error('problem with format of preprocessed files: impossible to find them. Did you preprocess the data?');
        end
        % extract number of scans for current run
        n_scansPerRun.(['CID',sub_nm])(iRun) = size(preprocessed_filenames,1);
        % load scans in matlabbatch
        switch iRun
            case 1
                matlabbatch1{iS}.spm.stats.fmri_spec.sess(1).scans = preprocessed_filenames;
            otherwise
                matlabbatch1{iS}.spm.stats.fmri_spec.sess(1).scans = [matlabbatch1{iS}.spm.stats.fmri_spec.sess(1).scans;...
                    preprocessed_filenames]; % add the files of the current run on each iteration
        end
        %% clear run name to avoid bugs across runs
        clear('subj_runFoldername_tmp');
    end % run loop for concatenating scans
    
    %% load regressors of interest. Group sessions depending on the DCM_mode selected
    % extract onsets + regressors
    matlabbatch1 = First_level_loadRegressors_DCM(matlabbatch1, GLMprm, study_nm, sub_nm, iS, runs, n_runs,...
        subj_behavior_folder, computerRoot, n_scansPerRun.(['CID',sub_nm]), TR, DCM_mode);
    
    %% global run parameters (rp movement file, etc.)
    % 1) concatenate all rp files together
    [mvmtFolder, mvmt_file_nm] = concatenate_rp_files(subj_scans_folder, subj_scan_folders_names, n_scansPerRun.(['CID',sub_nm]));
    movement_filePath = [mvmtFolder, mvmt_file_nm];
    
    % 2) load concatenated rp files inside the matlabbatch
    matlabbatch1{iS}.spm.stats.fmri_spec.sess(1).multi = {''};
    matlabbatch1{iS}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    matlabbatch1{iS}.spm.stats.fmri_spec.sess(1).multi_reg = {movement_filePath};
    matlabbatch1{iS}.spm.stats.fmri_spec.sess(1).hpf = 128;
    
    %% global parameters for subject batch
    matlabbatch1{iS}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    
    %% add temporal derivative or not
    if add_drv == 0
        matlabbatch1{iS}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    elseif add_drv == 1
        matlabbatch1{iS}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
    elseif add_drv == 2
        matlabbatch1{iS}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
    end
    
    %% other parameters
    matlabbatch1{iS}.spm.stats.fmri_spec.volt = 1;
    matlabbatch1{iS}.spm.stats.fmri_spec.global = 'None';
    switch grey_mask
        case 0 % filter with threshold here if no mask entered
            matlabbatch1{iS}.spm.stats.fmri_spec.mthresh = 0.8; % default value
        case 4 % filter with threshold here if no mask entered
            matlabbatch1{iS}.spm.stats.fmri_spec.mthresh = 0.5; % try lower value to include ventral striatum
        case 5 % filter with threshold here if no mask entered
            matlabbatch1{iS}.spm.stats.fmri_spec.mthresh = 0.3; % try lower value to include ventral striatum
        case 6 % filter with threshold here if no mask entered
            matlabbatch1{iS}.spm.stats.fmri_spec.mthresh = 0.1; % try lower value to include ventral striatum
        case {1,2,3,7} % no filter here since the mask will do the job
            matlabbatch1{iS}.spm.stats.fmri_spec.mthresh = -Inf; % implicitly masks the first level depending on the probability that the voxel is "relevant"
    end
    
    %% add grey mask or not
    switch grey_mask
        case {0,4,5,6}
            matlabbatch1{iS}.spm.stats.fmri_spec.mask = {''};
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
                matlabbatch1{iS}.spm.stats.fmri_spec.mask = {mask_file_path};
            else
                error('problem: mean anatomy file not found for this set of subjects');
            end
    end
    
    %% temporal autocorrelation
    switch autocorrel
        case 0
            matlabbatch1{iS}.spm.stats.fmri_spec.cvi = 'AR(1)';
        case 1
            matlabbatch1{iS}.spm.stats.fmri_spec.cvi = 'FAST';
    end
    
end % loop through subjects to load matlabbatch

%% checking => display or no checking => launch spm_fmri_concatenate + 1st level estimation
switch checking
    case 0 % not checking => run 1st level + spm_fmri_concatenate + 1st level estimation
        % run spm batch for 1st level
        spm_jobman('run',matlabbatch1);
        
        %% apply spm_fmri_concatenate on each SPM.m to correct for pooling all sessions into one
        % spm_fmri_concatenate will create one constant/session + it will adjust the high-pass filter and non-sphericity
        % estimates as if sessions were separate
        for iS = 1:NS
            sub_nm = subject_id{iS};
            
            % define working folders
            subj_folder             = [root, filesep, 'CID',sub_nm];
            subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep];
            % create folder to store the results for the current subject
            switch biasFieldCorr
                case 0
                    sm_folderName = [subj_analysis_folder 'functional', filesep,...
                        'preproc_sm_',num2str(preproc_sm_kernel),'mm_DCM',filesep];
                case 1
                    error('case with bias-field correction and DCM slice-timing preprocessing not ready yet');
            end
            if ~exist(sm_folderName,'dir')
                mkdir(sm_folderName);
            end
            [subPath] = fMRI_subFolder_DCM(sm_folderName, GLM, condition, DCM_mode);
            
            % apply spm_fmri_concatenate => will rewrite SPM.mat file
            spm_fmri_concatenate([subPath,'SPM.mat'], n_scansPerRun.(['CID',sub_nm]));
        end % loop through subjects
        
        %% estimate the model AFTER applying spm_fmri_concatenate
        % (important to do in that order but see https://www.fil.ion.ucl.ac.uk/spm/docs/wikibooks/Concatenation/
        % for more details as to why)
        matlabbatch2_estimate = cell(NS,1);
        for iS = 1:NS
            sub_nm = subject_id{iS};
            % extract path for folder where 1st level SPM.mat is stored
            switch biasFieldCorr
                case 0
                    sm_folderName = [root, filesep, 'CID',sub_nm, filesep,...
                        'fMRI_analysis',filesep 'functional', filesep,...
                        'preproc_sm_',num2str(preproc_sm_kernel),'mm_DCM',filesep];
                case 1
                    error('case with bias-field correction and DCM slice-timing preprocessing not ready yet');
            end
            [resultsFolderName] = fMRI_subFolder_DCM(sm_folderName, GLM, condition, DCM_mode);
            
            % load path into batch
            matlabbatch2_estimate{iS}.spm.stats.fmri_est.spmmat = {[resultsFolderName,'SPM.mat']};
            matlabbatch2_estimate{iS}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch2_estimate{iS}.spm.stats.fmri_est.method.Classical = 1; % perform classical model
        end
        
        % run spm batch for model estimation
        spm_jobman('run',matlabbatch2_estimate);
        
    case 1 % checking => display 1st level batch and stop there (avoid launching spm_fmri_concatenate + 1st level estimation, only for display: careful to not launch)
        % run spm batch for 1st level
        spm_jobman('interactive',matlabbatch1);
        warning('This is only for checking and display. Do NOT launch the batch as you need to apply spm_fmri_concatenate + 1st level estimation afterwards');
end % checking

end % function