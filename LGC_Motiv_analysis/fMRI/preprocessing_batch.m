function[] = preprocessing_batch(study_nm, sub_nm)
%[] = preprocessing_batch(study_nm, sub_nm)
% preprocessing for fMRI data
% enter subject identification in 'subject_id' (XXX corresponding to CID number)
%
% all data is stored in the fMRI_scans folders + the script also creates
% fMRI_analysis folder where the anatomical files are copied for the
% preprocessing + functional folder created for 1st level analysis
%
% Preprocessing entails 1) realignement, 2) co-registration, 3)
% segmentation, 4-5) normalisation in MNI space, 6) spatial smoothing
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
switch study_nm
    case 'fMRI_pilots' % pilots
        root = [fullfile(computerRoot,'fMRI_pilots'),filesep];
    case 'study1'
        root = [fullfile(computerRoot,'study1'),filesep];
    case 'study2'
        root = [fullfile(computerRoot,'study2'),filesep];
    case 'study2_pilots'
        root = [fullfile(computerRoot,'study2','pilots','fMRI_pilots'),filesep];
end

if ~exist('sub_nm','var') || isempty(sub_nm)
    subject_id = LGCM_subject_selection(study_nm);
else
    subject_id = {sub_nm};
end

%% smoothing kernel
smKernel = 8;

%% remove from the list participants for which preprocessing was already made
[subject_id, NS] = preproc_already_made(computerRoot, study_nm, subject_id, smKernel);

if NS >= 1
    %% give path for anatomical template
    spmTemplatePath = fullfile(spmFolderPath,'spm12','spm12','tpm','TPM.nii');
    
    
    %% define number of preprocessing steps
    nb_preprocessingSteps = 6;
    % nb of preprocessing steps per subject:
    % 1)realign-reslice
    % 2)coregistration functional - anatomical
    % 3) segmentation
    % 4) normalization functional
    % 5) normalization anatomical
    % 6) smoothing functional
    
    
    %% initialize batch
    matlabbatch = cell(nb_preprocessingSteps*NS,1);
    
    for iS = 1:NS % loop through subjects
        disp(['loading batch for preprocessing subject ',num2str(iS),'/',num2str(NS)]);
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
        subj_analysis_folder = [subFolder,'fMRI_analysis'];
        if exist(subj_analysis_folder,'dir') ~= 7
            mkdir(subj_analysis_folder);
        end
        newAnatFolder = [subFolder,'anatomical'];
        if exist(newAnatFolder,'dir') ~= 7
            mkdir(newAnatFolder);
        end
        fMRI_analysis_folder = [subFolder,'functional'];
        if exist(fMRI_analysis_folder,'dir') ~= 7
            mkdir(fMRI_analysis_folder);
        end
        subj_scans_folder = [root, sub_fullNm, filesep,'fMRI_scans', filesep];
        anat_folder = ls([subj_scans_folder,'*UNI-DEN*']);
        copyfile([subj_scans_folder,anat_folder], newAnatFolder); % copies the contempt of the anatomical folder
        
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
        end
        
        %% realignement
        preproc_step = 1; % first step of preprocessing
        realign_step = nb_preprocessingSteps*(iS-1) + preproc_step; % number depends also on the number of subjects
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.data = cell(n_runs, 1);
        for iRun = 1:n_runs
            matlabbatch{realign_step}.spm.spatial.realign.estwrite.data{iRun,1} = runFileNames.(['run_',num2str(iRun)]);
        end
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.weight = '';
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.which = [2 1];
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.interp = 4;
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.mask = 1;
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
        %% coregistration
        preproc_step = 2;
        coreg_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{coreg_step}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image',...
            substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
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
        matlabbatch{coreg_step}.spm.spatial.coreg.estimate.source = {[newAnatFolder, anat_file]};
        matlabbatch{coreg_step}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        %% segmentation
        preproc_step = 3;
        segm_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{segm_step}.spm.spatial.preproc.channel.vols = {[newAnatFolder,anat_file]};
        matlabbatch{segm_step}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{segm_step}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{segm_step}.spm.spatial.preproc.channel.write = [1 1];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).tpm = {[spmTemplatePath,',1']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).warped = [1 1]; % normalise grey matter and record modulated and unmodulated warped tissue
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).tpm = {[spmTemplatePath,',2']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).tpm = {[spmTemplatePath,',3']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).tpm = {[spmTemplatePath,',4']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).tpm = {[spmTemplatePath,',5']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).tpm = {[spmTemplatePath,',6']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{segm_step}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{segm_step}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{segm_step}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{segm_step}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{segm_step}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{segm_step}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{segm_step}.spm.spatial.preproc.warp.write = [1 1];
        %% normalization
        preproc_step = 4;
        normf_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{normf_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
        for iRun = 1:n_runs
            matlabbatch{normf_step}.spm.spatial.normalise.write.subj.resample(iRun) = cfg_dep(['Realign: Estimate & Reslice: Resliced Images (Sess ',num2str(iRun),')'], substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{iRun}, '.','rfiles'));
        end
        matlabbatch{normf_step}.spm.spatial.normalise.write.subj.resample(n_runs+1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
        matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
        %% normalization anatomical scan
        preproc_step = 5;
        norma_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{norma_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
        matlabbatch{norma_step}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
        matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
        %% smoothing
        preproc_step = 6;
        smooth_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{smooth_step}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{normf_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
        matlabbatch{smooth_step}.spm.spatial.smooth.fwhm = [smKernel smKernel smKernel];
        matlabbatch{smooth_step}.spm.spatial.smooth.dtype = 0;
        matlabbatch{smooth_step}.spm.spatial.smooth.im = 0;
        matlabbatch{smooth_step}.spm.spatial.smooth.prefix = 's';
        
        cd(root);
    end
    
    % display spm batch before running it
%     spm_jobman('interactive',matlabbatch);
    spm_jobman('run',matlabbatch);
    
    %% move files output from the last step
    for iS = 1:NS
        sub_nm = subject_id{iS};
        switch study_nm
            case {'study1','study2'}
                sub_fullNm = ['CID',sub_nm];
            case {'fMRI_pilots','study2_pilots'}
                sub_fullNm = sub_nm;
        end
        subj_scans_folder = [root, sub_fullNm, filesep,'fMRI_scans'];
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
        
        for iRun = 1:n_runs % loop through runs for 3 ratings, 3 choices 1D, 3 choices 2D runs
            runPath = [subj_scan_folders_names(iRun,:),filesep]; % go to run folder
            preproc_newFolder_nm = ['preproc_sm_',num2str(smKernel),'mm'];
            if ~exist([runPath,preproc_newFolder_nm],'dir')
                mkdir([runPath,preproc_newFolder_nm]);
            else
                error(['preprocessing folder ',preproc_newFolder_nm,' already exists for subject ',sub_fullNm,' run ',num2str(iRun)]);
                % note: something should be done above to avoid re-doing the
                % preprocessing for the subjects where it was already done
            end
            
            if strcmp(study_nm,'fMRI_pilots')
                if ismember(sub_nm,{'pilot_s1'})
                    filenames = ls([runPath,'*swrLGCM*.nii']);
                elseif ismember(sub_nm,{'pilot_s2'})
                    filenames = ls([runPath,'*swrrun*.nii']);
                elseif ismember(sub_nm,{'pilot_s3'})
                    filenames = ls([runPath,'*swrABNC*.img']);
                    filenames = [filenames; ls([runPath,'*swrABNC*.hdr'])];
                else
                    filenames = ls([runPath,'*swrCID*.nii']);
                end
            elseif strcmp(study_nm,'study2_pilots')
                if ismember(sub_nm,'fMRI_pilot1_AC')
                    filenames = ls([runPath,'*swrAC*.nii']);
                end
            else
                filenames = ls([runPath,'*swrCID*.nii']);
                %             error('please check the format (nii/img) and the start of the name of each run because it has to be stabilized now...');
            end
            
            % move files
            for iFile = 1:length(filenames)
                movefile([runPath,filenames(iFile,:)],...
                    [runPath,preproc_newFolder_nm]);
            end
        end % run loop
    end % subject loop
    
else
    disp(['All subjects have already been preprocessed with smoothing ',...
        'kernel of ',num2str(smKernel),'mm.']);
end % if all subjects have already been preprocessed, don't do anything

end % function