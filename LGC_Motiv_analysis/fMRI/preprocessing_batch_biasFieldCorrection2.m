function[] = preprocessing_batch_biasFieldCorrection2(study_nm, sub_nm)
%[] = preprocessing_batch_biasFieldCorrection2(study_nm, sub_nm)
% preprocessing for fMRI data similar to preprocessing_batch.m but will
% also include a step for correcting bias field inhomogeneities (the goal
% being to compensate for the fact that some regions have very low
% intensity of signal and are therefore dismissed by SPM filters, in
% particular the ventral striatum, globus pallidus and occipital cortex
% have much less signal at 7T. The bias field correction script is adapted from
% Antoine Lutti, Dominik Bach and Guillaume Flandin.
% enter subject identification in 'subject_id' (XXX corresponding to CID number)
%
% all data is stored in the fMRI_scans folders + the script also creates
% fMRI_analysis folder where the anatomical files are copied for the
% preprocessing + functional folder created for 1st level analysis
%
% Preprocessing entails:
% 1) realignement
% 2) bias-field correction of the EPI images based on the 1st EPI image of
% each run
% 3) co-registration
% 4) segmentation
% 5) normalisation in MNI space
% 6) spatial smoothing
%
% INPUTS
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
%
% sub_nm: subject name of the participant you want to preprocess. If left
% empty, by default the script will preprocess all the individuals of the
% study defined in 'study_nm' input.
%
% See also First_level_batch, contrasts_batch and
% Second_level_batch

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% define subjects and working directories
% [computerRoot, spmFolderPath] = LGCM_root_paths();
computerRoot = ['E:',filesep];
spmFolderPath = fullfile('C:','Users','clairis','Desktop');
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm_List = {'study1','study2','fMRI_pilots','study2_pilots'};
    study_nm_idx = listdlg('ListString',study_nm_List);
    study_nm = study_nm_List{study_nm_idx};
end
% switch study_nm
%     case 'fMRI_pilots' % pilots
%         root = [fullfile(computerRoot,'fMRI_pilots'),filesep];
%     case 'study1'
%         root = [fullfile(computerRoot,'study1'),filesep];
%     case 'study2'
%         root = [fullfile(computerRoot,'study2'),filesep];
%     case 'study2_pilots'
%         root = [fullfile(computerRoot,'study2','pilots','fMRI_pilots'),filesep];
% end
root = [fullfile(computerRoot,'test'),filesep];

if ~exist('sub_nm','var') || isempty(sub_nm)
    subject_id = LGCM_subject_selection(study_nm);
else
    subject_id = {sub_nm};
end

%% smoothing kernel
smKernel = 8;

%% remove from the list participants for which preprocessing was already made
% [subject_id, NS] = preproc_already_made(computerRoot, study_nm, subject_id, smKernel);
NS = length(subject_id);
if NS >= 1
    %% give path for anatomical template
    spmTemplatePath = fullfile(spmFolderPath,'spm12','spm12','tpm','TPM.nii');
    
    
    %% define number of preprocessing steps
    nb_preprocessingSteps_step1 = 1;
    % nb of preprocessing steps per subject:
    % 1)realign-reslice
    
    % Second step:
    % 1) compute bias field for each session of each subject
    
    nb_preprocessingSteps_step3 = 5;
    % Third step:
    % 1) coregistration functional - anatomical
    % 2) segmentation
    % 3) normalization functional
    % 4) normalization anatomical
    % 5) smoothing functional
    
    
    %% initialize first batch
    matlabbatch_realign = cell(nb_preprocessingSteps_step1*NS,1);
    
    for iS = 1:NS % loop through subjects
        disp(['loading re-alignement for preprocessing subject ',num2str(iS),'/',num2str(NS)]);
        sub_nm = subject_id{iS};
        switch study_nm
            case {'study1','study2'}
                sub_fullNm = ['CID',sub_nm];
            case {'fMRI_pilots','study2_pilots'}
                sub_fullNm = sub_nm;
        end
        subFolder = [root, sub_fullNm, filesep];
        % create working directories and copy anat. file inside
        % \fMRI_analysis\anatomical folder
        subj_analysis_folder = [subFolder,'fMRI_analysis',filesep];
        if ~exist(subj_analysis_folder,'dir')
            mkdir(subj_analysis_folder);
        end
        newAnatFolder = [subj_analysis_folder,'anatomical',filesep];
        if ~exist(newAnatFolder,'dir')
            mkdir(newAnatFolder);
        end
        fMRI_analysis_folder = [subj_analysis_folder,'functional',filesep];
        if ~exist(fMRI_analysis_folder,'dir')
            mkdir(fMRI_analysis_folder);
        end
        subj_scans_folder = [root, sub_fullNm, filesep,'fMRI_scans', filesep];
        anat_folder = ls([subj_scans_folder,'*UNI-DEN*']);
        copyfile([subj_scans_folder,anat_folder],...
            newAnatFolder); % copies the contempt of the anatomical folder
        
        %% extract folders where functional runs are stored
        subj_scan_folders_names = ls([subj_scans_folder,'*run*']); % takes all functional runs folders
        % remove AP/PA top-up corrective runs when they were performed (only 2
        % first pilots)
        if strcmp(study_nm,'fMRI_pilots') && ismember(sub_nm,{'pilot_s1','pilot_s2'})
            [subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names);
        end
        %% define number of sessions to analyze
        if strcmp(study_nm,'fMRI_pilots') &&...
                ismember(sub_nm, {'pilot_s1','pilot_s2','pilot_s3'}) % only 2 sessions for these pilots
            n_runs = 2;
        elseif strcmp(study_nm,'study1') &&...
                ismember(sub_nm,{'040'}) % fMRI had to be crashed during run 3
            n_runs = 2;
        else
            n_runs = 4;
        end
        for iRun = 1:n_runs % loop through runs
            runPath = [subj_scans_folder, subj_scan_folders_names(iRun,:)]; % extract run folder
            [filenames] = preproc_run_name_extraction(study_nm, sub_nm, runPath);
            runFileNames.(['run_',num2str(iRun)]) = filenames;
        end % run loop
        
        %% realignement
        preproc_step = 1; % first step of preprocessing
        realign_step = nb_preprocessingSteps_step1*(iS-1) + preproc_step; % number depends also on the number of subjects
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.data = cell(n_runs, 1);
        for iRun = 1:n_runs
            matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.data{iRun,1} = runFileNames.(['run_',num2str(iRun)]);
        end
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch_realign{realign_step}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    end % subject loop
    spm_jobman('run',matlabbatch_realign);
    
    % display info that this step is over
    disp(['Realignement of ',num2str(NS),' subjects has been correctly performed.']);
    
    %% apply bias field correction to all the re-aligned files
    for iS = 1:NS % loop through subjects
        disp(['performing bias-field correction for preprocessing subject ',num2str(iS),'/',num2str(NS)]);
        sub_nm = subject_id{iS};
        switch study_nm
            case {'study1','study2'}
                sub_fullNm = ['CID',sub_nm];
            case {'fMRI_pilots','study2_pilots'}
                sub_fullNm = sub_nm;
        end
        subj_scans_folder = [root, sub_fullNm, filesep,'fMRI_scans',filesep];
        % define number of sessions to check
        if strcmp(study_nm,'fMRI_pilots') &&...
                ismember(sub_nm, {'pilot_s1','pilot_s2','pilot_s3'}) % only 2 sessions for these pilots
            n_runs = 2;
        elseif strcmp(study_nm,'study1') &&...
                ismember(sub_nm,{'040'}) % fMRI had to be crashed during run 3
            n_runs = 2;
        else
            n_runs = 4;
        end
        
        for iRun = 1:n_runs
            % re-initialize the matlab batch for each run
            matlabbatch_biasField = cell(1,1);
            %% extract folders where functional runs are stored
            subj_scan_folders_names = ls([subj_scans_folder,'*run*']); % takes all functional runs folders
            % remove AP/PA top-up corrective runs when they were performed (only 2
            % first pilots)
            if strcmp(study_nm,'fMRI_pilots') && ismember(sub_nm,{'pilot_s1','pilot_s2'})
                [subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names);
            end
            % select re-aligned files for the current run
            runPath = [subj_scans_folder, subj_scan_folders_names(iRun,:)]; % extract run folder
            [realigned_filenames] = preproc_run_name_extraction(study_nm, sub_nm, runPath,'r');
            
            %% bias field estimation for correction using SPM segmentation (applied on the first image of each fMRI session)
            matlabbatch_biasField{1}.spm.tools.preproc8.channel.vols = realigned_filenames(1);
            matlabbatch_biasField{1}.spm.tools.preproc8.channel.write = [1 0];
            spm_jobman('run', matlabbatch_biasField);
            
            %% perform bias-field correction
            % extract bias field path
            bf = fullfile(strrep(runPath,' ',''), filesep, spm_select('List', runPath, '^BiasField.*\.nii$'));
            bfv = spm_vol(bf);
            BF = double(spm_read_vols(bfv));
            % perform correction and create 1 new file for each
            for iFile = 1:numel(realigned_filenames)
                % read file
                fn = [realigned_filenames{iFile}];
                V = spm_vol(fn);
                Y = spm_read_vols(V);
                % apply bias field
                Y = BF.*Y;
                % save file
                [pth, fnn, ext] = fileparts(fn);
                nfn = fullfile(pth, ['b', fnn, ext]);
                V.fname = strrep(nfn,',1','');
                spm_write_vol(V,Y);
            end
        end % run loop
    end % subject loop
    
    % display info that this step is over
    disp(['Bias-field correction of ',num2str(NS),' subjects has been correctly performed.']);
    
    %% final part of preprocessing
    for iS = 1:NS
        sub_nm = subject_id{iS};
        switch study_nm
            case {'study1','study2'}
                sub_fullNm = ['CID',sub_nm];
            case {'fMRI_pilots','study2_pilots'}
                sub_fullNm = sub_nm;
        end
        %% define number of sessions to analyze
        if strcmp(study_nm,'fMRI_pilots') &&...
                ismember(sub_nm, {'pilot_s1','pilot_s2','pilot_s3'}) % only 2 sessions for these pilots
            n_runs = 2;
        elseif strcmp(study_nm,'study1') &&...
                ismember(sub_nm,{'040'}) % fMRI had to be crashed during run 3
            n_runs = 2;
        else
            n_runs = 4;
        end
        subFolder = [root, sub_fullNm, filesep];
        subj_analysis_folder = [subFolder,'fMRI_analysis',filesep];
        newAnatFolder = [subj_analysis_folder,'anatomical', filesep];
        subj_scans_folder = [subFolder,'fMRI_scans',filesep];
        subj_scan_folders_names = ls([subj_scans_folder,'*run*']); % takes all functional runs folders
        % remove AP/PA top-up corrective runs when they were performed (only 2
        % first pilots)
        if strcmp(study_nm,'fMRI_pilots') && ismember(sub_nm,{'pilot_s1','pilot_s2'})
            [subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names);
        end
        
        matlabbatch_finalPreproc = cell(nb_preprocessingSteps_step3*NS,1);
        %% coregistration
        preproc_step = 1;
        coreg_step = nb_preprocessingSteps_step3*(iS-1) + preproc_step;
        % extract mean functional which will serve for co-registration of
        % the anatomical scan
        run1Path = [subj_scans_folder, strrep(subj_scan_folders_names(1,:),' ',''),filesep]; % in general mean file automatically saved inside run 1 folder
        meanFile = ls([run1Path,'mean*.nii']);
        matlabbatch_finalPreproc{coreg_step}.spm.spatial.coreg.estimate.ref = {[run1Path,meanFile]};
        if strcmp(study_nm,'fMRI_pilots') && ismember(sub_nm,{'pilot_s1'})
            anat_file = ls([newAnatFolder,'LGCM_*.nii']);
        elseif strcmp(study_nm,'fMRI_pilots') && ismember(sub_nm,{'pilot_s3'})
            anat_file = ls([newAnatFolder,'ABNC_*.img']);
        elseif strcmp(study_nm,'study2_pilots') && ismember(sub_nm,{'fMRI_pilot1_AC'})
            anat_file = ls([newAnatFolder,'AC*UNI-DEN.nii']);
        else
            %         anat_file = ls('mp2rage_*.nii');
            anat_file = ls([newAnatFolder,'CID*.nii']);
        end
        matlabbatch_finalPreproc{coreg_step}.spm.spatial.coreg.estimate.source = {[newAnatFolder, anat_file]};
        matlabbatch_finalPreproc{coreg_step}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch_finalPreproc{coreg_step}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch_finalPreproc{coreg_step}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch_finalPreproc{coreg_step}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch_finalPreproc{coreg_step}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        %% segmentation
        preproc_step = 2;
        segm_step = nb_preprocessingSteps_step3*(iS-1) + preproc_step;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.channel.vols = {[newAnatFolder,anat_file]};
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.channel.write = [1 1];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(1).tpm = {[spmTemplatePath,',1']};
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(1).warped = [1 1]; % normalise grey matter and record modulated and unmodulated warped tissue
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(2).tpm = {[spmTemplatePath,',2']};
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(3).tpm = {[spmTemplatePath,',3']};
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(4).tpm = {[spmTemplatePath,',4']};
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(5).tpm = {[spmTemplatePath,',5']};
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(6).tpm = {[spmTemplatePath,',6']};
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch_finalPreproc{segm_step}.spm.spatial.preproc.warp.write = [1 1];
        %% normalization
        preproc_step = 3;
        normf_step = nb_preprocessingSteps_step3*(iS-1) + preproc_step;
        matlabbatch_finalPreproc{normf_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations',...
            substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','fordef', '()',{':'}));
        % extract run files
        rfiles_to_normalise_to_MNI = [];
        for iRun = 1:n_runs % select re-aligned files for the current run
            runPath = [subj_scans_folder, subj_scan_folders_names(iRun,:)]; % extract run folder
            [realign_biasFieldCorr_run_files_tmp] = preproc_run_name_extraction(study_nm, sub_nm, runPath,'br');
            rfiles_to_normalise_to_MNI = [rfiles_to_normalise_to_MNI; realign_biasFieldCorr_run_files_tmp];
        end % run loop
        rfiles_to_normalise_to_MNI = [rfiles_to_normalise_to_MNI; [run1Path,meanFile]];
        matlabbatch_finalPreproc{normf_step}.spm.spatial.normalise.write.subj.resample = rfiles_to_normalise_to_MNI;
        matlabbatch_finalPreproc{normf_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch_finalPreproc{normf_step}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch_finalPreproc{normf_step}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch_finalPreproc{normf_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
        %% normalization anatomical scan
        preproc_step = 4;
        norma_step = nb_preprocessingSteps_step3*(iS-1) + preproc_step;
        matlabbatch_finalPreproc{norma_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations',...
            substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','fordef', '()',{':'}));
        matlabbatch_finalPreproc{norma_step}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Segment: Bias Corrected (1)',...
            substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        matlabbatch_finalPreproc{norma_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch_finalPreproc{norma_step}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
        matlabbatch_finalPreproc{norma_step}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch_finalPreproc{norma_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
        %% smoothing
        preproc_step = 5;
        smooth_step = nb_preprocessingSteps_step3*(iS-1) + preproc_step;
        matlabbatch_finalPreproc{smooth_step}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)',...
            substruct('.','val', '{}',{normf_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        matlabbatch_finalPreproc{smooth_step}.spm.spatial.smooth.fwhm = [smKernel smKernel smKernel];
        matlabbatch_finalPreproc{smooth_step}.spm.spatial.smooth.dtype = 0;
        matlabbatch_finalPreproc{smooth_step}.spm.spatial.smooth.im = 0;
        matlabbatch_finalPreproc{smooth_step}.spm.spatial.smooth.prefix = 's';
        
        cd(root);
    end % subject loop
    % display spm batch before running it
    spm_launch_or_display = 'run'; % 'run' or 'interactive'
    spm_jobman(spm_launch_or_display,matlabbatch_finalPreproc);
    
    % display info that this step is over
    if strcmp(spm_launch_or_display,'run')
        disp(['Last part of pre-processing (co-registration, normalization in ',...
            'MNI space and smoothing) of ',num2str(NS),' subjects has been correctly performed.']);
    end
    
    %% move files output from the last step
    if strcmp(spm_launch_or_display,'run')
        biasFieldCorr = 1;
        move_preproc_files_to_saveFolder(root, study_nm,...
            subject_id, NS, smKernel, biasFieldCorr);
    end % move files
else
    disp(['All subjects have already been preprocessed with smoothing ',...
        'kernel of ',num2str(smKernel),'mm.']);
end % if all subjects have already been preprocessed, don't do anything

end % function