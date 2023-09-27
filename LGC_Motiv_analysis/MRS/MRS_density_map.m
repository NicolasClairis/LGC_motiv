function[] = MRS_density_map(study_nm, condition)

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
[subject_id, NS_total] = LGCM_subject_selection(study_nm, condition);

%% working directories
study_path = fullfile('E:',study_nm);
MRS_voxel_results_folder = fullfile(study_path,'MRS_mask');
% create folder to store the results
if ~exist(MRS_voxel_results_folder,'dir')
    mkdir(MRS_voxel_results_folder);
end
%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% main parameters
switch study_nm
    case 'study1'
        MRS_ROIs = {'dmPFC','ai'};
    otherwise
        error([study_nm,' not ready yet']);
end
nROIs = length(MRS_ROIs);
% initialize batch
matlabbatch = cell(nROIs*2,1);

batch_idx = 0;
for iROI = 1:nROIs
    ROI_nm = MRS_ROIs{iROI};
    ROI_binary_nm = ['bw',ROI_nm,'.nii'];
    
    %% extract individual file names
    NS.(ROI_nm) = NS_total;
    MRS_voxel_perSub = cell(NS.(ROI_nm),1);
    % add all the MRS voxel files for all the subjects
    jS = 0;
    for iS = 1:NS_total
        sub_fullNm = ['CID',subject_id{iS}];
        sub_MRS_voxels_folder = [fullfile(study_path,sub_fullNm,...
            'MRS','MRS_voxels'),filesep];
        sub_MRS_voxel_file_fullNm = [sub_MRS_voxels_folder, ROI_binary_nm];
        if exist(sub_MRS_voxel_file_fullNm,'file')
            jS = jS + 1;
            MRS_voxel_perSub(jS) = {sub_MRS_voxel_file_fullNm};
        else % remove one subject if no ROI
            NS.(ROI_nm) = NS.(ROI_nm) - 1;
            MRS_voxel_perSub(end) = [];
        end
    end
    
    NS_str_ROI = num2str(NS.(ROI_nm));
    
    %% 1) take mean MRS voxel across the selected participants
    batch_idx = batch_idx + 1;
    mean_filename = ['mean_',ROI_nm,'_MRS_voxel__',condition,'_',NS_str_ROI,'_subjects'];
    matlabbatch{batch_idx}.spm.util.imcalc.input = MRS_voxel_perSub;
    matlabbatch{batch_idx}.spm.util.imcalc.output = mean_filename;
    matlabbatch{batch_idx}.spm.util.imcalc.outdir = {MRS_voxel_results_folder};
    matlabbatch{batch_idx}.spm.util.imcalc.expression = 'mean(X)';
    matlabbatch{batch_idx}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{batch_idx}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{batch_idx}.spm.util.imcalc.options.mask = 0;
    matlabbatch{batch_idx}.spm.util.imcalc.options.interp = 1;
    matlabbatch{batch_idx}.spm.util.imcalc.options.dtype = 4; % read as 4D image (necessary to be able to use the 'mean(X)' function of SPM)
    
    %% 2) compute density map (ie % of subjects present for each voxel)
    batch_idx = batch_idx + 1;
    density_filename = ['density_map_',ROI_nm,'_MRS_voxel__',condition,'_',NS_str_ROI,'_subjects'];
    matlabbatch{batch_idx}.spm.util.imcalc.input = MRS_voxel_perSub;
    matlabbatch{batch_idx}.spm.util.imcalc.output = density_filename;
    matlabbatch{batch_idx}.spm.util.imcalc.outdir = {MRS_voxel_results_folder};
    matlabbatch{batch_idx}.spm.util.imcalc.expression = ['sum(X)./',NS_str_ROI];
    matlabbatch{batch_idx}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{batch_idx}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{batch_idx}.spm.util.imcalc.options.mask = 0;
    matlabbatch{batch_idx}.spm.util.imcalc.options.interp = 1;
    matlabbatch{batch_idx}.spm.util.imcalc.options.dtype = 4; % read as 4D image (necessary to be able to use the 'sum(X)' function of SPM)
end % loop over ROIs

%% launch SPM script
spm_jobman('run',matlabbatch);
end % function