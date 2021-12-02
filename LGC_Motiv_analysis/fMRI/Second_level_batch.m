function[] = Second_level_batch(study_nm, GLM)
% script to launch second level on LGC Motivation studies
%
% INPUTS
% study_nm : name of the study to analyze
%
% GLM: number of the GLM


%% clear workspace
close all; clc;

%% initiate SPM
spm('defaults','fmri');
spm_jobman('initcfg');

%% check the batch before launching the script?
checking = 0;

%% value of the smoothing during preprocessing?
preproc_sm_kernel = 6;

%% define study and list of subjects to include
% define study
if ~exist('study_nm','var') || isempty(study_nm)
    study_names = {'fMRI_pilots','study1','study2'};
    study_nm_idx = listdlg('ListString',study_names);
    study_nm = study_names{study_nm_idx};
end

% define subjects
[subject_id, NS] = LGCM_subject_selection(study_nm);
NS_str = num2str(NS);

%% working directories
computer_root = LGCM_root_paths();
% scripts_folder = fullfile(computer_root,'GitHub','LGC_motiv','LGC_Motiv_analysis','fMRI');
% addpath(scripts_folder);
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        root = fullfile(computer_root,'study1');
    case 'study2'
        root = fullfile(computer_root,'study2');
end

%% define GLM
if ~exist('GLM','var') || isempty(GLM) || GLM <= 0
    GLM = spm_input('GLM number?',1,'e',[]);
end
GLM_str = num2str(GLM);
GLMprm = which_GLM(GLM);

%% create results folder
results_folder = [root,filesep,'Second_level',filesep,...
    'GLM',GLM_str,'_',NS_str,'subs_preprocSm',num2str(preproc_sm_kernel),'mm',filesep];
if exist(results_folder,'dir') ~= 7
    mkdir(results_folder);
end

%% 1) take mean anatomy across participants
batch_idx = 0;

% extract parameters for mean calculation
mean_filename = ['mean_anat_',NS_str,'_subjects'];
wms_anat = cell(NS,1);
% add all the anat files for all the subjects
% extract anat EPI
for iS = 1:NS
    sub_anat_folder = [root,filesep,'CID',subject_id{iS},filesep,...
        'fMRI_analysis',filesep,'anatomical',filesep];
    wms_anat_name = ls([sub_anat_folder,'wm*']);
    wms_anat(iS) = {[sub_anat_folder, wms_anat_name]};
end

% spm calculation of the mean
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.util.imcalc.input = wms_anat;
matlabbatch{batch_idx}.spm.util.imcalc.output = mean_filename;
matlabbatch{batch_idx}.spm.util.imcalc.outdir = {results_folder(1:end-1)};
matlabbatch{batch_idx}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{batch_idx}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{batch_idx}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{batch_idx}.spm.util.imcalc.options.mask = 0;
matlabbatch{batch_idx}.spm.util.imcalc.options.interp = 1;
matlabbatch{batch_idx}.spm.util.imcalc.options.dtype = 4;

%% perform the second level

% load all contrasts of interest
con_names = LGCM_contrasts(study_nm, subject_id{1}, GLM, computer_root);
n_con = length(con_names);

% loop over contrasts
for iCon = 1:n_con
    con_str = num2str(iCon);

    % directory for concatenated contrast
    conFolder_nm = strrep(con_names{iCon},' ','_');
    conFolder_nm = strrep(conFolder_nm,':','');
    conFolder_nm = [results_folder, conFolder_nm];
    mkdir(conFolder_nm);

    % start second level:
    batch_idx = batch_idx + 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {conFolder_nm};
    % list of inputs
    conlist = cell(NS,1); % 1 con per EPI-subject
    % extract contrasts per subject
    for iS = 1:NS
        sub_nm = subject_id{iS};

        % check incompatibility between some GLM and some subjects to see
        % if you need to redefine the list of subjects included in the
        % analysis
        checkGLM_and_subjectIncompatibility(study_nm, sub_nm, GLMprm);

        subject_folder = [root,filesep,'CID',sub_nm, filesep, 'fMRI_analysis' filesep,...
            'functional' filesep, 'preproc_sm_',num2str(preproc_sm_kernel),'mm',filesep...
            'GLM',GLM_str, filesep];
        if iCon < 10
            conlist(iS) = {[subject_folder,'con_000',con_str,'.nii,1']};
        elseif iCon >= 10 && iCon < 100
            conlist(iS) = {[subject_folder,'con_00',con_str,'.nii,1']};
        elseif iCon >= 100 && iCon < 1000
            conlist(iS) = {[subject_folder,'con_0',con_str,'.nii,1']};
        elseif iCon >= 1000 && iCon < 10000
            conlist(iS) = {[subject_folder,'con_',con_str,'.nii,1']};
        end
    end

    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t1.scans = conlist;

    % default parameters
    matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{batch_idx}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{batch_idx}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{batch_idx}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.globalm.glonorm = 1;

    % model estimation
    batch_model_rtg = batch_idx;
    batch_idx = batch_idx + 1;
    matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File',...
        substruct('.','val', '{}',{batch_model_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;
    % contrast
    batch_estm_rtg = batch_idx;
    batch_idx = batch_idx + 1;
    matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File',...
        substruct('.','val', '{}',{batch_estm_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.name = con_names{iCon};
    matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{batch_idx}.spm.stats.con.delete = 0;

end % loop over contrasts

    cd(results_folder);

%% display spm batch before running it or run it directly
if exist('checking','var') && checking == 1
    spm_jobman('interactive',matlabbatch);
    %     spm_jobman('run',matlabbatch);
elseif ~exist('checking','var') || isempty(checking) || checking == 0
    spm_jobman('run',matlabbatch);
end


end % function