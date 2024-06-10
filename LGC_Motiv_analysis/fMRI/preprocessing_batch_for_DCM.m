function[] = preprocessing_batch_for_DCM(study_nm, sub_nm, checking)
%[] = preprocessing_batch_for_DCM(study_nm, sub_nm, checking)
% preprocessing for fMRI data. Similar to preprocessing_batch, but adding
% slice timing procedure to adapt the data for DCM processing.
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
% checking: (0) runs the batch immediately, (1) shows you the batch for testing
% but avoid running it as it will not move the resulting files in that case
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

if (~exist('sub_nm','var') || isempty(sub_nm))
    switch checking
        case 0
            condition = subject_condition;
            subject_id = LGCM_subject_selection(study_nm, condition);
        case 1 % only 1 subject for testing
            subject_id = {'001'};
    end
else
    subject_id = {sub_nm};
end

%% smoothing kernel
smKernel = 8;

%% acquisition parameters (important for slice-timing in particular)
nslices = 63; % number of slices during acquisition
TR = 2.00; % repetition time in seconds (ie time between first slice of current scan and first slice of next scan)
TA = TR - (TR/nslices); % time between first and last slice within one scan
% "We need the time to acquire all but the last slice. SPM calls this TA
% (time of acquisition). This odd parameter comes from deep in the history of the SPM slice-timing routine."
% https://bic-berkeley.github.io/psych-214-fall-2016/spm_slice_timing_solution.html
refslice_sliceT = floor(nslices/2); % reference slice to use (generally the middle one)

% define order of acquisition of the slices in one volume
% slice_order = [1:2:nslices, 2:2:nslices]; % interleaved ascending order
% slice_order = [nslices:(-2):1, (nslices-1):(-2):1]; % interleaved descending order
slice_order = NaN(1,nslices); % interleaved middle-top, ie. scanning top->middle and middle->bottom in parallel
for iSlice = 1:nslices
    slice_order(iSlice) = round((nslices - iSlice)/2 + (rem((nslices - iSlice),2)*(nslices - 1)/2)) + 1;
end % loop over slices

%% decide whether you want to display the preprocessing script at the end 
% or nor run it directly
switch checking
    case 0
        spm_launch_or_display = 'run'; % 'run' or 'interactive'
    case 1
        spm_launch_or_display = 'interactive'; % 'run' or 'interactive'
end % checking data

%% remove from the list participants for which preprocessing was already made
switch spm_launch_or_display
    case 'run'
        [subject_id, NS] = preproc_already_made(computerRoot, study_nm, subject_id, smKernel);
    otherwise
        warning(['Careful spm_launch_or_display is set to "interactive" ',...
            'and no filter has been done to remove subjects already ',...
            'preprocessed, so be cautious before launching it, but feel free ',...
            'to explore as much as you want.']);
        NS = length(subject_id);
end

if NS >= 1
    %% give path for anatomical template
    spmTemplatePath = fullfile(spmFolderPath,'spm12','spm12','tpm','TPM.nii');
    
    
    %% define number of preprocessing steps
    nb_preprocessingSteps = 7;
    % nb of preprocessing steps per subject:
    % 1) realign-reslice
    % 2) slice-timing
    % 3) coregistration functional - anatomical
    % 4) segmentation
    % 5) normalization functional
    % 6) normalization anatomical
    % 7) smoothing functional
    
    
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
        subj_analysis_folder = [subFolder,'fMRI_analysis',filesep];
        if exist(subj_analysis_folder,'dir') ~= 7
            mkdir(subj_analysis_folder);
        end
        newAnatFolder = [subj_analysis_folder,'anatomical',filesep];
        if exist(newAnatFolder,'dir') ~= 7
            mkdir(newAnatFolder);
        end
        fMRI_analysis_folder = [subj_analysis_folder,'functional',filesep];
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
        %% slice-timing (mandatory step for DCM)
        preproc_step = 2;
        sliceT_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        for iRun = 1:n_runs
            matlabbatch{sliceT_step}.spm.temporal.st.scans{iRun}(1) = cfg_dep(['Realign: Estimate: Realigned Images (Sess ',num2str(iRun),')'],...
                substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
                substruct('.','sess', '()',{iRun}, '.','cfiles'));
        end % run loop
        matlabbatch{sliceT_step}.spm.temporal.st.nslices = nslices;
        matlabbatch{sliceT_step}.spm.temporal.st.tr = TR;
        matlabbatch{sliceT_step}.spm.temporal.st.ta = TA;
        matlabbatch{sliceT_step}.spm.temporal.st.so = slice_order;
        matlabbatch{sliceT_step}.spm.temporal.st.refslice = refslice_sliceT;
        matlabbatch{sliceT_step}.spm.temporal.st.prefix = 'a';
        
        %% coregistration of anatomical file to functional space
        preproc_step = 3;
        coreg_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{coreg_step}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image',...
            substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','rmean'));
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
        %% segmentation anatomical scan
        preproc_step = 4;
        segm_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{segm_step}.spm.spatial.preproc.channel.vols = {[newAnatFolder,anat_file]};
        matlabbatch{segm_step}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{segm_step}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{segm_step}.spm.spatial.preproc.channel.write = [1 1];
        % grey matter
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).tpm = {[spmTemplatePath,',1']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).warped = [1 1]; % normalise grey matter and record modulated and unmodulated warped tissue
        % white matter
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).tpm = {[spmTemplatePath,',2']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).warped = [1 1]; % normalise white matter and record modulated and unmodulated warped tissue
        % CSF
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).tpm = {[spmTemplatePath,',3']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).warped = [0 0];
        % dura matter
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).tpm = {[spmTemplatePath,',4']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).warped = [0 0];
        % skull
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).tpm = {[spmTemplatePath,',5']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).warped = [0 0];
        % air and any other abnormal tissue
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).tpm = {[spmTemplatePath,',6']};
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).warped = [0 0];
        % other parameters
        matlabbatch{segm_step}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{segm_step}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{segm_step}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{segm_step}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{segm_step}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{segm_step}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{segm_step}.spm.spatial.preproc.warp.write = [1 1]; % record inverse and forward deformation fields
        %% normalization functional images in MNI
        preproc_step = 5;
        normf_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{normf_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations',...
            substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','fordef', '()',{':'}));
        for iRun = 1:n_runs
            matlabbatch{normf_step}.spm.spatial.normalise.write.subj.resample(iRun) = cfg_dep(['Slice Timing: Slice Timing Corr. Images (Sess ',num2str(iRun),')'],...
                substruct('.','val', '{}',{sliceT_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
                substruct('.','sess', '()',{iRun}, '.','files'));
        end
        matlabbatch{normf_step}.spm.spatial.normalise.write.subj.resample(n_runs+1) = cfg_dep('Realign: Estimate & Reslice: Mean Image',...
            substruct('.','val', '{}',{sliceT_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','rmean'));
        matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
        %% normalization anatomical scan in MNI
        preproc_step = 6;
        norma_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{norma_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations',...
            substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{norma_step}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Segment: Bias Corrected (1)',...
            substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
        matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
        %% smoothing functional images
        preproc_step = 7;
        smooth_step = nb_preprocessingSteps*(iS-1) + preproc_step;
        matlabbatch{smooth_step}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)',...
            substruct('.','val', '{}',{normf_step}, '.',...
            'val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('()',{1}, '.','files'));
        matlabbatch{smooth_step}.spm.spatial.smooth.fwhm = [smKernel smKernel smKernel];
        matlabbatch{smooth_step}.spm.spatial.smooth.dtype = 0;
        matlabbatch{smooth_step}.spm.spatial.smooth.im = 0;
        matlabbatch{smooth_step}.spm.spatial.smooth.prefix = 's';
        
        cd(root);
    end % subject loop
    
    % display spm batch before running it?
    spm_jobman(spm_launch_or_display,matlabbatch);
    
    %% move files output from the last step (if you launched the script)
    if strcmp(spm_launch_or_display,'run')
        move_preproc_files_to_saveFolder(root, study_nm,...
            subject_id, NS, smKernel, 0);
    end % move files
else
    disp(['All subjects have already been preprocessed with smoothing ',...
        'kernel of ',num2str(smKernel),'mm.']);
end % if all subjects have already been preprocessed, don't do anything

end % function