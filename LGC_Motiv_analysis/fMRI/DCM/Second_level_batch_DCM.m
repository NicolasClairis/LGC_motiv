function[] = Second_level_batch_DCM(GLM, condition, gender, biasFieldCorr)
% Second_level_batch_DCM(GLM, condition, gender, biasFieldCorr)
% script to launch second level on LGC Motivation studies for DCM first
% level.
%
% INPUTS
% GLM: number of the GLM
%
% condition: define subjects and runs to include
%
% gender:
% 'all': all subjects by default
% 'males': remove females from the list
% 'females': remove males from the list
%
% biasFieldCorr: use bias-field corrected images (1) or not (0)? By default
% will not use bias-field corrected images
%

%% clear workspace
close all; clc;

%% initiate SPM
spm('defaults','fmri');
spm_jobman('initcfg');

%% check the batch before launching the script?
checking = 0;

%% value of the smoothing during preprocessing?
preproc_sm_kernel = 8;

%% use bias-field corrected files or not?
if ~exist('biasFieldCorr','var') || ~ismember(biasFieldCorr,[0,1])
    biasFieldCorr = 0;
end

%% which DCM mode to use? will define how sessions are grouped
[DCM_mode] = which_DCM_mode_for_GLM;

%% define study and list of subjects to include
% define study
if ~exist('study_nm','var')
    %     study_names = {'fMRI_pilots','study1','study2'};
    %     study_nm_idx = listdlg('ListString',study_names);
    %     study_nm = study_names{study_nm_idx};
    study_nm = 'study1'; % by default
end

% define subjects
if ~exist('condition','var') || ~strcmp(condition(1:4),'fMRI')
    condition = subject_condition;
end
if ~exist('gender','var') ||...
        isempty(gender) ||...
        ~ismember(gender,{'all','males','females'})
    gender = 'all';
end
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, gender);
NS_str = num2str(NS);

%% working directories
% computer_root = LGCM_root_paths();
computer_root = ['E:',filesep];
% scripts_folder = fullfile(computer_root,'GitHub','LGC_motiv','LGC_Motiv_analysis','fMRI');
% addpath(scripts_folder);
switch study_nm
    case 'fMRI_pilots'
        studyRoot = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        studyRoot = fullfile(computer_root,'study1');
    case 'study2'
        studyRoot = fullfile(computer_root,'study2');
end

%% define GLM
if ~exist('GLM','var') || isempty(GLM) || GLM <= 0
    GLM = spm_input('GLM number?',1,'e',[]);
end
GLM_str = num2str(GLM);
GLMprm = which_GLM(GLM);

%% create results folder
switch biasFieldCorr
    case 0
        biasField_sufix = '';
    case 1
        biasField_sufix = '_with_BiasFieldCorrection';
end
mainPath = [studyRoot,filesep,'Second_level',biasField_sufix,'_DCM',filesep];
if ~exist(mainPath,'dir')
    mkdir(mainPath);
end
switch condition
    case {'fMRI'}
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noSatRunSub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatRunSubs_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noSatTaskSub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskSubs_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noSatTaskSub_noSatRun'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskSubs_noSatRun_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noMoveSub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noFullMvmtSub_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noMoveSub_bis'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtSubLenient_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noMoveSub_ter'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtSubStringent_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noSatTaskSub_noMove_bis_Sub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskNoMvmtSub_DCM_mode',num2str(DCM_mode),filesep];
    case {'fMRI_noSatTask','fMRI_noSatTask_bayesianMdl'}
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTask_DCM_mode',num2str(DCM_mode),filesep];
    case {'fMRI_noSatRun','fMRI_noSatRun_bayesianMdl'}
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatRun_DCM_mode',num2str(DCM_mode),filesep];
    case {'fMRI_noSatRun_choiceSplit_Elvl','fMRI_noSatTaskSub_noSatRun_choiceSplit_Elvl'}
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatRun_Elvl_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noMove_bis'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtRunLenient_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noMove_ter'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtRunStringent_DCM_mode',num2str(DCM_mode),filesep];
    case {'fMRI_noSatTask_noMove_bis','fMRI_noSatTask_noMove_bis_bayesianMdl'}
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskNoMvmtRun_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noSatTaskSub_noMoveSub_noSatRun'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatNoMoveTaskSubs_noSatRun_DCM_mode',num2str(DCM_mode),filesep];
    case 'fMRI_noSatTaskSub_noMoveSub_noSatRun_noMoveRun'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatNoMoveTaskSubs_noSatNoMoveRun_DCM_mode',num2str(DCM_mode),filesep];
    otherwise
        error(['folder not ready yet for the condition ',condition]);
end

% create folder to store the results
if ~exist(results_folder,'dir')
    mkdir(results_folder);
else
    error(['Folder with the name ',results_folder,' already exists.']);
end

%% 1) take mean anatomy across participants
batch_idx = 0;

% extract parameters for mean calculation
mean_filename = ['mean_anat_',NS_str,'_subjects'];
wms_anat = cell(NS,1);
% add all the anat files for all the subjects
% extract anat EPI
for iS = 1:NS
    sub_anat_folder = [studyRoot,filesep,'CID',subject_id{iS},filesep,...
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
switch condition
    case 'fMRI_noSatRun_choiceSplit_Elvl_bis' % in this case, subject_id{1} does not include Em => need to adapt to include Em
        fcon_names_perSub = Fcon_list_DCM(study_nm, subject_id{16}, GLM, DCM_mode, computer_root, preproc_sm_kernel, condition, biasFieldCorr);
    otherwise
        fcon_names_perSub = Fcon_list_DCM(study_nm, subject_id{1}, GLM, DCM_mode, computer_root, preproc_sm_kernel, condition, biasFieldCorr);
end
n_con = length(fcon_names_perSub);

% loop over contrasts
for iCon = 1:n_con
    current_fcon_nm = fcon_names_perSub{iCon};

    % directory for concatenated contrast
    conFolder_nm = strrep(current_fcon_nm,' ','_');
    conFolder_nm = strrep(conFolder_nm,':','');
    conFolder_nm = strrep(conFolder_nm,'_|','_abs_');
    conFolder_nm = strrep(conFolder_nm,'|','');
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
        checkGLM_and_subjectIncompatibility(study_nm, sub_nm, condition, GLMprm);
        subject_main_folder = [studyRoot,filesep,'CID',sub_nm, filesep,...
            'fMRI_analysis' filesep, 'functional' filesep,...
            'preproc_sm_',num2str(preproc_sm_kernel),'mm',biasField_sufix,'_DCM',filesep];
        subject_main_folder = fMRI_subFolder_DCM(subject_main_folder, GLM, condition);
        if isempty(subject_main_folder)
            error(['condition ',condition,' not planned yet. Please add it.']);
        end
        
        %% adapt contrast index since some conditions and contrasts are missing for some subjects
        fcon_names_perSub = Fcon_list_DCM(study_nm, sub_nm, GLM,...
            DCM_mode, computer_root, preproc_sm_kernel, condition, biasFieldCorr);
        if sum(strcmp(current_fcon_nm, fcon_names_perSub)) > 0
            % extract index (for this subject) of the current contrast
            jCon = find(strcmp(current_fcon_nm, fcon_names_perSub));
            
            % extract name for this particular subject of the contrast of
            % interest
            [con_str] = conNumber2conName(jCon);
            conlist(iS) = {[subject_main_folder,'con_',con_str,'.nii,1']};
        end % in case contrast exists for the current subject
    end % subject loop
    
    % remove empty subjects from the list for the current contrast
    conlist(cellfun(@isempty,conlist)) = [];

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

    %% model estimation
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
    matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.name = current_fcon_nm;
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