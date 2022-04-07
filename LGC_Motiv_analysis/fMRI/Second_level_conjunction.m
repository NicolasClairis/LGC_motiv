%% conjunction betwnneen two contrasts defined by the user

%% clear
clear all; close all; clc;

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% working directories
% computer_root = LGCM_root_paths();
computerRoot = ['E:',filesep];
if ~exist('study_nm','var') || isempty(study_nm)
    %     study_names = {'study1','study2','fMRI_pilots'};
    %     study_nm_idx = listdlg('ListString',study_names);
    %     study_nm = study_names{study_nm_idx};
    study_nm = 'study1'; % by default
end
switch study_nm
    case 'fMRI_pilots'
        studyRoot = fullfile(computerRoot,'fMRI_pilots');
    case 'study1'
        studyRoot = fullfile(computerRoot,'study1');
    case 'study2'
        studyRoot = fullfile(computerRoot,'study2');
end

%% define GLM to check
GLM = spm_input('GLM number?',1,'e');
GLM_str = num2str(GLM);

%% define preprocessing kernel to use
preproc_sm_kernel = 8;
% preproc_sm_kernel = spm_input('smoothing kernel to use?',1,'e','8');

%% extract subjects
[subject_id, NS] = LGCM_subject_selection('study1');
NS_str = num2str(NS);

%% define contrasts on which conjunction will be operated
[con_names, con_vector] = LGCM_contrasts(study_nm, subject_id{iS}, GLM, computerRoot, preproc_sm_kernel);

selectedCon = NaN(1,2);
for iCon = 1:2
    selectedCon(iCon) = spm_input();
end % contrast loop

%% extract folder of interest
GLM_folder = [computerRoot, filesep, 'GLM', GLM_str];

%% initialize
batch_idx = 0;

%% mean anatomical file
% extract parameters for mean calculation
mean_filename = ['mean_anat_',NS_str,'_subjects'];
wms_anat = cell(NS,1);
% add all the anat files for all the subjects
% extract anat EPI
for iS = 1:NS
    wms_anat_name = ls([studyRoot,subject_id{iS}, filesep, 'fMRI_analysis',filesep,'anatomical',filesep 'wms*']);
    wms_anat(iS) = {[studyRoot,subject_id{iS}, filesep, 'fMRI_analysis', filesep, 'anatomical', filesep, wms_anat_name]};
end

% spm calculation of the mean
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.util.imcalc.input = wms_anat;
matlabbatch{batch_idx}.spm.util.imcalc.output = mean_filename;
matlabbatch{batch_idx}.spm.util.imcalc.outdir = {taskResults_folder(1:end-1)};
matlabbatch{batch_idx}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{batch_idx}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{batch_idx}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{batch_idx}.spm.util.imcalc.options.mask = 0;
matlabbatch{batch_idx}.spm.util.imcalc.options.interp = 1;
matlabbatch{batch_idx}.spm.util.imcalc.options.dtype = 4;

%% 2nd level concatenation of EPIs
%% enter contrasts in spm
batch_idx = batch_idx + 1;

% start second level:
matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {GLM_folder};
% list of inputs
conlist = cell(NS,con_length); % 1 con per EPI-subject
% extract MBBjune2016 contrasts
for iContrast = 1:con_length
    con_str = num2str(conOfInterest(iContrast));
    for iS = 1:NS        
        subject_folder = [sub_firstLevel_folder, 'GLM',GLM_str,'_megaconcatenation', filesep];
        if conOfInterest(iContrast) < 10
            conlist(iS,iContrast) = {[subject_folder,s,'con_000',con_str,'.nii,1']};
        elseif conOfInterest(iContrast) >= 10 && conOfInterest(iContrast) < 100
            conlist(iS,iContrast) = {[subject_folder,s,'con_00',con_str,'.nii,1']};
        elseif conOfInterest(iContrast) >= 100 && conOfInterest(iContrast) < 1000
            conlist(iS,iContrast) = {[subject_folder,s,'con_0',con_str,'.nii,1']};
        end
    end
end

%% Be careful t2.scans (for two-sample t.test)
if which_technique == 1
    
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.scans1 = conlist(:,1);
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.scans2 = conlist(:,2);
    
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.dept = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.variance = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.gmsca = 0;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.ancova = 0;
    matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
elseif which_technique == 2
    allConList = [conlist(:,1); conlist(:,2)];
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.mreg.scans = allConList;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.mreg.incint = 0;
    
    % one covariate per contrast of interest
    % 0 everywhere except for the con of interest = 1
    for iContrast = 1:con_length
        tmp = zeros(NS*con_length,1);
        tmp(1+(iContrast-1)*NS:iContrast*NS) = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov(iContrast).c = tmp;
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov(iContrast).cname = list_con{conOfInterest(iContrast)};
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov(iContrast).iCFI = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov(iContrast).iCC = 5;
    end
end

matlabbatch{batch_idx}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{batch_idx}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.glonorm = 1;

%% model estimation
batch_model_rtg = batch_idx;
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{batch_model_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;

%% t.contrast
batch_estm_rtg = batch_idx;
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{batch_estm_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

for iContrast = 1:con_length
    matlabbatch{batch_idx}.spm.stats.con.consess{iContrast}.tcon.name = list_con{conOfInterest(iContrast)};
    matlabbatch{batch_idx}.spm.stats.con.consess{iContrast}.tcon.weights = [zeros(1, iContrast - 1), 1, zeros(1,con_length - iContrast)];
    matlabbatch{batch_idx}.spm.stats.con.consess{iContrast}.tcon.sessrep = 'none';
    matlabbatch{batch_idx}.spm.stats.con.delete = 0;
end

%% f.contrast
fIdentityConne = eye(con_length);
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{batch_estm_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.name = 'Fcon_conjunction';
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.weights = fIdentityConne;
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.delete = 0;

%%
cd(results_folder);

%% display spm batch before running it
spm_jobman('interactive',matlabbatch);
% spm_jobman('run',matlabbatch);