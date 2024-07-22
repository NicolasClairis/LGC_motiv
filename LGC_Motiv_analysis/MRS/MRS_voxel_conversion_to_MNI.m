function[] = MRS_voxel_conversion_to_MNI(condition, study_nm)
% MRS_voxel_conversion_to_MNI converts MRS voxels from native space to MNI
% space using SPM normalisation batch with MRS MRI deformation field and
% MRS voxel as input. MRS_voxel_conversion_to_MNI needs
% MRS_MRI_conversion_to_MNI.m to be runned first in order to extract the
% information relative to the deformation field from native space to MNI
% space.
% MRS_voxel_conversion_to_MNI will produce both a probabilistic map 
% (w-VoxelName.nii) and a binary mask (bw-VoxelName.nii) (based on 
% bThreshold threshold used to consider a voxel to be part of the MRS voxel
% or not).
%
% INPUTS
% condition: to know which subjects to include (will be asked if left
% empty)
%
% study_nm: which study? ('study1' by default)
%
% See also MRS_MRI_conversion_to_MNI.m

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
% binarization threshold
bThreshold = 0.1;
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
        ROI_binary_file = [sub_ROI_folder, 'bw',ROI_file];
        
        %% filter case when ROI not extracted (because signal too bad or else)
        if exist(ROI_fullfile,'file')
            %% launch normalization if the voxel exists
            jStep = jStep + 1;
            matlabbatch{jStep}.spm.spatial.normalise.write.subj.def(1) = {[sub_MRI_folder,forward_defField_file]};
            matlabbatch{jStep}.spm.spatial.normalise.write.subj.resample(1) = {ROI_fullfile};
            matlabbatch{jStep}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                78 76 85];
            matlabbatch{jStep}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
            matlabbatch{jStep}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{jStep}.spm.spatial.normalise.write.woptions.prefix = 'w';
            
            %% binarize mask
            jStep = jStep + 1;
            matlabbatch{jStep}.spm.util.imcalc.input(1) = cfg_dep('Image Calculator: ImCalc Computed Image: output',...
                substruct('.','val', '{}',{jStep-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
            matlabbatch{jStep}.spm.util.imcalc.output           = ROI_binary_file;
            matlabbatch{jStep}.spm.util.imcalc.outdir           = {sub_ROI_folder};
            matlabbatch{jStep}.spm.util.imcalc.expression       = ['i1>',num2str(bThreshold)]; % keep only voxels with signal intensity higher than bThreshold
            matlabbatch{jStep}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
            matlabbatch{jStep}.spm.util.imcalc.options.dmtx     = 0;
            matlabbatch{jStep}.spm.util.imcalc.options.mask     = 0;
            matlabbatch{jStep}.spm.util.imcalc.options.interp   = 1;
            matlabbatch{jStep}.spm.util.imcalc.options.dtype    = 4;
        end % ROI existing filter
    end % ROI loop
end % subject loop

%% launch batch
spm_jobman('run',matlabbatch);

end % function