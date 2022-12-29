function[] = Second_level_batch_metabolites_comparison(GLM, condition, gender)
% Second_level_batch_metabolites_comparison(GLM, condition, gender)
% script to launch second level on LGC Motivation studies
%
% INPUTS
% GLM: number of the GLM
%
% condition: define subjects and runs to include (you can also use
% subject_condition.mat to select the condition)
%
% gender:
% 'all': all subjects by default
% 'males': remove females from the list
% 'females': remove males from the list

%% clear workspace
close all; clc;

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

%% ask which metabolite to focus on?
[low_met_subs, high_met_subs, metabolite_nm, MRS_ROI_nm] = medSplit_metabolites(study_nm, subject_id);
NS_low = sum(low_met_subs);
NS_high = sum(high_met_subs);
%% initiate SPM
spm('defaults','fmri');
spm_jobman('initcfg');

%% check the batch before launching the script?
checking = 0;

%% value of the smoothing during preprocessing?
preproc_sm_kernel = 8;

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
mainPath = [studyRoot,filesep,'Second_level',filesep];
MRS_split_nm = [MRS_ROI_nm,'_',metabolite_nm,'MedSplit'];
switch condition
    case {'fMRI'}
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_',MRS_split_nm,filesep];
    case 'fMRI_noSatRunSub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatRunSubs_',MRS_split_nm,filesep];
    case 'fMRI_noSatTaskSub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskSubs_',MRS_split_nm,filesep];
    case 'fMRI_noMoveSub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noFullMvmtSub_',MRS_split_nm,filesep];
    case 'fMRI_noMoveSub_bis'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtSubLenient_',MRS_split_nm,filesep];
    case 'fMRI_noMoveSub_ter'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtSubStringent_',MRS_split_nm,filesep];
    case 'fMRI_noSatTaskSub_noMove_bis_Sub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskNoMvmtSub_',MRS_split_nm,filesep];
    case 'fMRI_noSatTask'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTask_',MRS_split_nm,filesep];
    case 'fMRI_noSatRun'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatRun_',MRS_split_nm,filesep];
    case 'fMRI_noMove_bis'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtRunLenient_',MRS_split_nm,filesep];
    case 'fMRI_noMove_ter'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtRunStringent_',MRS_split_nm,filesep];
    case 'fMRI_noSatTask_noMove_bis'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskNoMvmtRun_',MRS_split_nm,filesep];
    otherwise
        error(['folder not ready yet for the condition ',condition]);
end

% create folder to store the results
if ~exist(results_folder,'dir')
    mkdir(results_folder);
else
    error(['run folder with the name ',results_folder,' already exists.']);
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
con_names = LGCM_contrasts(study_nm, subject_id{1}, GLM,...
    computer_root, preproc_sm_kernel, condition);
n_con = length(con_names);

% loop over contrasts
for iCon = 1:n_con
    current_con_nm = con_names{iCon};

    % directory for concatenated contrast
    conFolder_nm = strrep(current_con_nm,' ','_');
    conFolder_nm = strrep(conFolder_nm,':','');
    conFolder_nm = [results_folder, conFolder_nm];
    mkdir(conFolder_nm);

    % start second level:
    batch_idx = batch_idx + 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {conFolder_nm};
    % list of inputs
    conlist_low = cell(NS_low,1); % 1 con per EPI-subject
    conlist_high = cell(NS_high,1); % 1 con per EPI-subject
    % re-initialize counters
    jS_low = 0;
    jS_high = 0;
    % extract contrasts per subject
    for iS = 1:NS
        sub_nm = subject_id{iS};
        
        % check incompatibility between some GLM and some subjects to see
        % if you need to redefine the list of subjects included in the
        % analysis
        checkGLM_and_subjectIncompatibility(study_nm, sub_nm, condition, GLMprm);
        subject_main_folder = [studyRoot,filesep,'CID',sub_nm, filesep, 'fMRI_analysis' filesep,...
            'functional' filesep, 'preproc_sm_',num2str(preproc_sm_kernel),'mm',filesep];
        subject_main_folder = fMRI_subFolder(subject_main_folder, GLM, condition);
        if isempty(subject_main_folder)
            error(['condition ',condition,' not planned yet. Please add it.']);
        end
        
        %% adapt contrast index since some conditions and contrasts are missing for some subjects
        con_names_perSub = LGCM_contrasts(study_nm, sub_nm, GLM,...
            computer_root, preproc_sm_kernel, condition);
        if sum(strcmp(current_con_nm, con_names_perSub)) > 0
            % extract index (for this subject) of the current contrast
            jCon = find(strcmp(current_con_nm, con_names_perSub));
            
            % extract name for this particular subject of the contrast of
            % interest
            [con_str] = conNumber2conName(jCon);
            if low_met_subs(iS) == 1 && high_met_subs(iS) == 0
                jS_low = jS_low + 1;
                conlist_low(jS_low) = {[subject_main_folder,'con_',con_str,'.nii,1']};
            elseif low_met_subs(iS) == 0 && high_met_subs(iS) == 1
                jS_high = jS_high + 1;
                conlist_high(jS_high) = {[subject_main_folder,'con_',con_str,'.nii,1']};
            else
                error('Who''s that dude with weird metabolites?');
            end
        end % in case contrast exists for the current subject
    end % subject loop
    
    % remove empty subjects from the list for the current contrast
    conlist_low(cellfun(@isempty,conlist_low)) = [];
    conlist_high(cellfun(@isempty,conlist_high)) = [];

    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.scans1 = conlist_low;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.scans2 = conlist_high;

    % default parameters
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.dept = 0;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.variance = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.gmsca = 0;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.ancova = 0;
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
    matlabbatch{batch_idx}.spm.stats.con.consess{1}.tcon.name = current_con_nm;
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