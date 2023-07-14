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

error(['Plutot que de regrouper grey + white individually, ',...
    'create one mask for each at the group level and then check ',...
    'if proba of grey OR white matter is > threshold => include ',...
    'while regrouping would imply to multiply = average = interaction']);

%% initiate SPM
spm('defaults','fmri');
spm_jobman('initcfg');

%% main parameters
% define smoothing kernel for the last step
sm_kernel = 8;
% define threshold to use for probability to be in the grey matter
proba_threshold = 0.50;
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

%% select subjects
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% prepare anatomy
[greyM_anatFiles, whiteM_anatFiles,...
    anatFiles] = deal(cell(NS,1));
for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    subFolder = fullfile(rootPath, sub_nm, 'fMRI_analysis','anatomical');
    % extract individual grey matter mask
    greyM_anatFile_tmp = ls([subFolder,filesep,'mwc1*.nii']);
    greyM_anatFiles(iS,:) = {[subFolder,filesep,greyM_anatFile_tmp]};
    % extract individual white matter mask
    whiteM_anatFile_tmp = ls([subFolder,filesep,'mwc2*.nii']);
    whiteM_anatFiles(iS,:) = {[subFolder,filesep,whiteM_anatFile_tmp]};
    
    % extract individual grey matter mask
    anatFiles(iS,:) = {['mw',sub_nm,'_grey_and_white_matter.nii']};
    
    %% group grey and white matter together for each subject
    matlabbatch{batch_avg}.spm.util.imcalc.input            = [greyM_anatFiles(iS,:);...
        whiteM_anatFiles(iS,:)];
    matlabbatch{batch_avg}.spm.util.imcalc.output           = anatFiles(iS,:);
    matlabbatch{batch_avg}.spm.util.imcalc.outdir           = {subFolder};
    matlabbatch{batch_avg}.spm.util.imcalc.expression       = 'mean(X)';
    matlabbatch{batch_avg}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
    matlabbatch{batch_avg}.spm.util.imcalc.options.dmtx     = 1; % data matrix set to 1 allows to pool all files together and average all voxels together using the expression mean(X)
    matlabbatch{batch_avg}.spm.util.imcalc.options.mask     = 0;
    matlabbatch{batch_avg}.spm.util.imcalc.options.interp   = 1;
    matlabbatch{batch_avg}.spm.util.imcalc.options.dtype    = 4;
end % subject loop
mean_filename = ['mean_greyAndWhiteM_proba_',num2str(NS),'_subjects_',condition];
binary_filename = ['bmean_greyAndWhiteM_',...
    num2str(NS),'_subjects_',condition,...
    '_',num2str(proba_threshold*100),'percentGreyWhiteM'];

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