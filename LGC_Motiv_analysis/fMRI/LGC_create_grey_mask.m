function[]=LGC_create_grey_mask(study_nm, condition)
% []=LGC_create_grey_mask(study_nm, condition)
% LGC_create_grey_mask will create a grey mask for study defined in
% study_nm and subjects based on the condition entered. If those
% parameters are not defined, they will be asked.
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
rootPath = [fullfile('E:',study_nm),filesep];
outDir = fullfile(rootPath,'grey_matter_mask');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% select subjects
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% prepare anatomy
anatFiles = cell(NS,1);
for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    subFolder = fullfile(rootPath, sub_nm, 'fMRI_analysis','anatomical');
    anatFile_tmp = ls([subFolder,filesep,'mwc1*.nii']);
    anatFiles(iS,:) = {[subFolder,filesep,anatFile_tmp]};
end % subject loop
mean_filename = ['mean_greyM_proba_',num2str(NS),'_subjects_',condition];
binary_filename = ['bmean_greyM_',...
    num2str(NS),'_subjects_',condition,...
    '_',num2str(proba_threshold*100),'percentGreyM'];

%% average all anatomical files together
batch_avg = 1;
matlabbatch{batch_avg}.spm.util.imcalc.input            = anatFiles;
matlabbatch{batch_avg}.spm.util.imcalc.output           = mean_filename;
matlabbatch{batch_avg}.spm.util.imcalc.outdir           = {outDir};
matlabbatch{batch_avg}.spm.util.imcalc.expression       = 'mean(X)';
matlabbatch{batch_avg}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
matlabbatch{batch_avg}.spm.util.imcalc.options.dmtx     = 1; % data matrix set to 1 allows to pool all files together and average all voxels together using the expression mean(X)
matlabbatch{batch_avg}.spm.util.imcalc.options.mask     = 0;
matlabbatch{batch_avg}.spm.util.imcalc.options.interp   = 1;
matlabbatch{batch_avg}.spm.util.imcalc.options.dtype    = 4;
%% binarize resulting grey mask
batch_binarize = 2;
matlabbatch{batch_binarize}.spm.util.imcalc.input(1) = cfg_dep('Image Calculator: ImCalc Computed Image: output', substruct('.','val', '{}',{batch_avg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{batch_binarize}.spm.util.imcalc.output           = binary_filename;
matlabbatch{batch_binarize}.spm.util.imcalc.outdir           = {outDir};
matlabbatch{batch_binarize}.spm.util.imcalc.expression       = ['i1>',num2str(proba_threshold)]; % keep only voxels with probability of being in the grey matter higher than proba threshold%
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