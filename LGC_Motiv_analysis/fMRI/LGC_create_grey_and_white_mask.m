function[]=LGC_create_grey_and_white_mask(study_nm, condition)
% []=LGC_create_grey_and_white_mask(study_nm, condition)
% LGC_create_grey_and_white_mask will create a grey + white matter mask for
% study defined in study_nm and subjects based on the condition entered. If
% those parameters are not defined, they will be asked.
%
% INPUTS
% study_nm: study name
%
% condition: condition
%
%

%% initiate SPM
spm('defaults','fmri');
spm_jobman('initcfg');

%% main parameters
% % define smoothing kernel for the last step
% sm_kernel = 8;
% define threshold to use for probability to be in the grey matter
proba_threshold = 0.00; % for 50%, write 0.50
%% check inputs
if ~exist('study_nm','var') || isempty(study_nm)
    %    error('study name not defined');
    study_nm = 'study1';
end

if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
%% working directories
rootPath = [fullfile('E:',study_nm),filesep];
outDir = fullfile(rootPath,'greyAndWhite_matter_mask');
if ~exist(outDir,'dir')
    mkdir(outDir);
end
greyM_folder = fullfile(rootPath,'grey_matter_mask');
whiteM_folder = fullfile(rootPath,'white_matter_mask');

%% select subjects
[~, NS] = LGCM_subject_selection(study_nm,condition);

%% prepare main files of interest
% input files
greyM_binary_filename = [greyM_folder,filesep,...
    'bmean_greyM_',...
    num2str(NS),'_subjects_',condition,...
    '_',num2str(proba_threshold*100),'percentGreyM.nii'];
whiteM_binary_filename = [whiteM_folder,filesep,...
    'bmean_whiteM_',...
    num2str(NS),'_subjects_',condition,...
    '_',num2str(proba_threshold*100),'percentWhiteM.nii'];

% resulting files
whiteM_greyM_pool_filename = ['pool_greyAndWhiteM_',...
    num2str(NS),'_subjects_',condition,...
    '_',num2str(proba_threshold*100),'percentM'];
binary_filename = ['bmean_greyAndWhiteM_',...
    num2str(NS),'_subjects_',condition,...
    '_',num2str(proba_threshold*100),'percentM'];

%% pool grey and white matter together (overlapping voxels will have weird non-binary values)
batch_avg = 1;
matlabbatch{batch_avg}.spm.util.imcalc.input            = {[greyM_binary_filename,',1'];...
    [whiteM_binary_filename,',1']};
matlabbatch{batch_avg}.spm.util.imcalc.output           = whiteM_greyM_pool_filename;
matlabbatch{batch_avg}.spm.util.imcalc.outdir           = {outDir};
matlabbatch{batch_avg}.spm.util.imcalc.expression       = 'sum(X)';
matlabbatch{batch_avg}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
matlabbatch{batch_avg}.spm.util.imcalc.options.dmtx     = 1; % data matrix set to 1 allows to pool all files together and average all voxels together using the expression mean(X)
matlabbatch{batch_avg}.spm.util.imcalc.options.mask     = 0;
matlabbatch{batch_avg}.spm.util.imcalc.options.interp   = 1;
matlabbatch{batch_avg}.spm.util.imcalc.options.dtype    = 4;
%% produce binary mask where any voxel which is either part of 
% grey matter or of white matter with a probability higher than 
% proba_thresholdis included
batch_binarize = 2;
% matlabbatch{2}.spm.util.imcalc.input(1) = cfg_dep('Image Calculator: ImCalc Computed Image: pool_greyAndWhiteM_65_subjects_fMRI_noSatTaskSub_noSatRun_20percentWhiteM', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{batch_binarize}.spm.util.imcalc.input(1) = cfg_dep('Image Calculator: ImCalc Computed Image: output',...
    substruct('.','val', '{}',{batch_avg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{batch_binarize}.spm.util.imcalc.output           = binary_filename;
matlabbatch{batch_binarize}.spm.util.imcalc.outdir           = {outDir};
% keep voxels with proba of being grey OR white matter > proba_threshold
% threshold here is at zero because the inputs are already binarized within
% each dimension (grey and white matter)
matlabbatch{batch_binarize}.spm.util.imcalc.expression       = 'i1>0'; 
matlabbatch{batch_binarize}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
matlabbatch{batch_binarize}.spm.util.imcalc.options.dmtx     = 0;
matlabbatch{batch_binarize}.spm.util.imcalc.options.mask     = 0;
matlabbatch{batch_binarize}.spm.util.imcalc.options.interp   = 1;
matlabbatch{batch_binarize}.spm.util.imcalc.options.dtype    = 4;
% %% smooth a bit the resulting mask to avoid spurious voxels
% batch_smooth = 3;
% matlabbatch{batch_smooth}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: output', substruct('.','val', '{}',{batch_binarize}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
% matlabbatch{batch_smooth}.spm.spatial.smooth.fwhm = [sm_kernel sm_kernel sm_kernel];
% matlabbatch{batch_smooth}.spm.spatial.smooth.dtype = 0;
% matlabbatch{batch_smooth}.spm.spatial.smooth.im = 0;
% matlabbatch{batch_smooth}.spm.spatial.smooth.prefix = 's';

%% perform the batch
% spm_jobman('interactive',matlabbatch);
spm_jobman('run',matlabbatch);

end % function