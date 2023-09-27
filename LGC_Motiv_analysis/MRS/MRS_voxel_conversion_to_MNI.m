function[] = MRS_voxel_conversion_to_MNI(condition, study_nm)
% MRS_voxel_conversion_to_MNI converts MRS voxel from native space to MNI
% space
%
% INPUTS
% condition: to know which subjects to include (will be asked if left
% empty)
%
% study_nm: which study? ('study1' by default)


%% subject selection
% condition
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
% study name
if ~exist('study_nm','var')
    study_nm = 'study1';
end
% list of subjects
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
study_path = fullfile('E:',study_nm);

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% main parameters
MRS_ROIs = {'dmPFC','ai'};
nROIs = length(MRS_ROIs);
% initialize batch
matlabbatch = cell(1,1);

%% subject loop
jStep = 0;
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_fullNm = ['CID',sub_nm];
    for iROI = 1:nROIs
        MRS_ROI_nm = MRS_ROIs{iROI};
        
        %% extract forward (native space=>MNI) deformation field
        switch sub_nm
            case {'021','056','088'}
                switch MRS_ROI_nm
                    case 'dmPFC'
                        sub_MRI_folder = [fullfile(study_path, sub_fullNm,'MRS','MRI','dmPFC_MRI'),filesep];
                    case 'ai'
                        sub_MRI_folder = [fullfile(study_path, sub_fullNm,'MRS','MRI','ai_MRI'),filesep];
                end
            otherwise
                sub_MRI_folder = [fullfile(study_path, sub_fullNm,'MRS','MRI'),filesep];
        end
        forward_defField_file = ls([sub_MRI_folder,'y_CID*UNI-DEN.nii']);
        
        %% extract MRS ROI
        sub_ROI_folder = [fullfile(study_path, sub_fullNm,'MRS','MRS_voxels'),filesep];
        switch MRS_ROI_nm
            case 'dmPFC'
                ROI_file = 'dmpfc.nii';
            case 'ai'
                ROI_file = 'ai.nii';
        end
        ROI_fullfile = [sub_ROI_folder, ROI_file];
        
        %% launch normalization if the voxel exists
        if exist(ROI_fullfile,'file')
            jStep = jStep + 1;
            matlabbatch{jStep}.spm.spatial.normalise.write.subj.def(1) = {[sub_MRI_folder,forward_defField_file]};
            matlabbatch{jStep}.spm.spatial.normalise.write.subj.resample(1) = {ROI_fullfile};
            matlabbatch{jStep}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                78 76 85];
            matlabbatch{jStep}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
            matlabbatch{jStep}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{jStep}.spm.spatial.normalise.write.woptions.prefix = 'w';
        end
    end % ROI loop
end % subject loop

%% launch batch
spm_jobman('run',matlabbatch);

end % function