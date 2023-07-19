function[] = redo_segmentation(study_nm, sub_nm)
% script aiming at relaunching all the segmentation for study 1 but this
% time also extracting the white matter files and not just the grey matter
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% sub_nm: subject name

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
    condition = subject_condition;
    subject_id = LGCM_subject_selection(study_nm, condition);
else
    subject_id = {sub_nm};
end

%% smoothing kernel
smKernel = 8;

%% decide whether you want to display the preprocessing script at the end 
% or nor run it directly
spm_launch_or_display = 'interactive'; % 'run' or 'interactive'

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
        
        
        %% coregistration
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
        %% segmentation
        preproc_step = 1;
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
        matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).warped = [1 1]; % normalise white matter and record modulated and unmodulated warped tissue
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
        
        cd(root);
    end
    
    % display spm batch before running it?
    spm_jobman(spm_launch_or_display,matlabbatch);
    
    %% move files where they need to be sent
else
    disp(['All subjects have already been preprocessed with smoothing ',...
        'kernel of ',num2str(smKernel),'mm.']);
end % if all subjects have already been preprocessed, don't do anything

end % function