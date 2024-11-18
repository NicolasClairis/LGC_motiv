% conjunction between two contrasts defined by the user

%% clear
clear; close all; clc;

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

%% use bias-field corrected maps?
biasFieldCorr = 0;

%% use preprocessing optimized for DCM (with slice-timing correction)?
DCM_GLM = 0;

%% extract subjects
[condition] = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);
NS_str = num2str(NS);

%% define contrasts on which conjunction will be operated
switch condition
    case 'fMRI_noSatRun_choiceSplit_Elvl_bis' % in this case, subject_id{1} does not include Em => need to adapt to include Em
        con_names = LGCM_contrasts(study_nm, subject_id{16}, GLM,...
            computerRoot, preproc_sm_kernel, condition, biasFieldCorr, DCM_GLM);
    otherwise
        con_names = LGCM_contrasts(study_nm, subject_id{1}, GLM,...
            computerRoot, preproc_sm_kernel, condition, biasFieldCorr, DCM_GLM);
end
n_cons = length(con_names);

nConsForConj = 2;
selectedCon = NaN(1, nConsForConj);
[selectedConNames, selectedConShortNames] = deal(cell(1,nConsForConj));
for iCon = 1:2
    selectedCon(iCon) = spm_input(['Please select contrast ',num2str(iCon),...
        ' for conjunction:'],...
        1,'m',con_names,1:n_cons,0);
    selectedConNames{iCon} = con_names{selectedCon(iCon)};
    selectedConShortNames{iCon} = spm_input(['Short name for: ',...
        con_names{selectedCon(iCon)}],1,'s');
end % contrast loop
% create name for resulting conjunction
conj_name = [selectedConShortNames{1},'_CONJ_',selectedConShortNames{2}];
% replace spaces by '_' if some spaces were left
conj_name = strrep(conj_name,' ','_');

%% extract folder of interest
switch biasFieldCorr
    case 0
        resultsFolder = [studyRoot, filesep,'Second_level',filesep,...
            'GLM', GLM_str,'_conjunction_',NS_str,'subs',filesep];
    case 1
        resultsFolder = [studyRoot, filesep,'Second_level_with_BiasFieldCorrection',filesep,...
            'GLM', GLM_str,'_conjunction_',NS_str,'subs',filesep];
end
if ~exist(resultsFolder,'dir')
    mkdir(resultsFolder);
end
conResultsFolder = [resultsFolder,filesep,conj_name];
if ~exist(conResultsFolder,'dir')
    mkdir(conResultsFolder);
else
    error([conResultsFolder,' already exists. ',...
        'Please rename folders to avoid confusion.']);
end
%% initialize
batch_idx = 0;

%% mean anatomical file
% extract parameters for mean calculation
mean_filename = ['mean_anat_',NS_str,'_subjects'];
wms_anat = cell(NS,1);
% add all the anat files for all the subjects
% extract anat EPI
for iS = 1:NS
    sub_nm = subject_id{iS};
    subAnatPath = [studyRoot,filesep,'CID',sub_nm, filesep,...
        'fMRI_analysis',filesep,'anatomical',filesep];
    wms_anat_name = ls([subAnatPath, 'wm*']);
    wms_anat(iS) = {[subAnatPath, wms_anat_name]};
end

% spm calculation of the mean
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.util.imcalc.input = wms_anat;
matlabbatch{batch_idx}.spm.util.imcalc.output = mean_filename;
matlabbatch{batch_idx}.spm.util.imcalc.outdir = {resultsFolder(1:end-1)};
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
matlabbatch{batch_idx}.spm.stats.factorial_design.dir = {conResultsFolder};
% list of inputs
conlist_allSubs = cell(NS, nConsForConj); % 1 con per EPI-subject
% extract contrasts for each subject
for iCon = 1:nConsForConj
    current_con_nm = selectedConNames{iCon};
    for iS = 1:NS
        sub_nm = subject_id{iS};
        % working directory
        switch biasFieldCorr
            case 0
                subfMRIPath = [studyRoot,filesep,'CID',sub_nm, filesep,...
                    'fMRI_analysis',filesep,'functional',filesep,...
                    'preproc_sm_',num2str(preproc_sm_kernel),'mm',filesep];
            case 1
                subfMRIPath = [studyRoot,filesep,'CID',sub_nm, filesep,...
                    'fMRI_analysis',filesep,'functional',filesep,...
                    'preproc_sm_',num2str(preproc_sm_kernel),'mm_with_BiasFieldCorrection',filesep];
        end
        subject_main_folder = fMRI_subFolder(subfMRIPath, GLM, condition);
        
        % extract contrast names for current subject
        con_names_perSub = LGCM_contrasts(study_nm, sub_nm, GLM,...
            computerRoot, preproc_sm_kernel, condition, biasFieldCorr);
        if sum(strcmp(current_con_nm, con_names_perSub)) > 0
            % extract index (for this subject) of the current contrast
            jCon = find(strcmp(current_con_nm, con_names_perSub));
            % extract name for this particular subject of the contrast of
            % interest
            [con_str] = conNumber2conName(jCon);
            conlist_allSubs(iS,iCon) = {[subject_main_folder,'con_',con_str,'.nii,1']};
        end
    end
end

%% remove subjects where one of the contrast is missing
okSubs = false(1,NS);
for iS = 1:NS
    okSubs(iS) = (~isempty(conlist_allSubs{iS,1}).*~isempty(conlist_allSubs{iS,2})) == 1;
end % loop through subjects
conlist = conlist_allSubs(okSubs,:);
NS_ok = sum(okSubs);
if NS_ok < NS
    warning(['Only ',num2str(NS_ok),'/',num2str(NS),' subjects with ok data for conjunction']);
end

%% Be careful t2.scans (for two-sample t.test)
which_technique = 1;
if which_technique == 1 % 2-sample t.test (unpaired) allowing to group the contrasts for the conjunction together
    
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.scans1 = conlist(:,1);
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.scans2 = conlist(:,2);
    
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.dept = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.variance = 1;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.gmsca = 0;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.t2.ancova = 0;
    matlabbatch{batch_idx}.spm.stats.factorial_design.cov = struct('c', {},...
        'cname', {}, 'iCFI', {}, 'iCC', {});
elseif which_technique == 2
    allConList = [conlist(:,1); conlist(:,2)];
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.mreg.scans = allConList;
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});
    matlabbatch{batch_idx}.spm.stats.factorial_design.des.mreg.incint = 0;
    
    % one covariate per contrast of interest
    % 0 everywhere except for the con of interest = 1
    for iCon = 1:con_length
        tmp = zeros(NS*n_cons,1);
        tmp(1+(iCon-1)*NS:iCon*NS) = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov(iCon).c = tmp;
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov(iCon).cname = selectedConNames{iCon};
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov(iCon).iCFI = 1;
        matlabbatch{batch_idx}.spm.stats.factorial_design.cov(iCon).iCC = 5;
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
matlabbatch{batch_idx}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File',...
    substruct('.','val', '{}',{batch_model_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{batch_idx}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{batch_idx}.spm.stats.fmri_est.method.Classical = 1;

%% t.contrast
batch_estm_rtg = batch_idx;
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File',...
    substruct('.','val', '{}',{batch_estm_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

for iCon = 1:nConsForConj
    matlabbatch{batch_idx}.spm.stats.con.consess{iCon}.tcon.name = con_names{selectedCon(iCon)};
    matlabbatch{batch_idx}.spm.stats.con.consess{iCon}.tcon.weights = [zeros(1, iCon - 1), 1, zeros(1,nConsForConj - iCon)];
    matlabbatch{batch_idx}.spm.stats.con.consess{iCon}.tcon.sessrep = 'none';
    matlabbatch{batch_idx}.spm.stats.con.delete = 0;
end

%% f.contrast
fIdentityConne = eye(nConsForConj);
batch_idx = batch_idx + 1;
matlabbatch{batch_idx}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File',...
    substruct('.','val', '{}',{batch_estm_rtg}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.name = 'Fcon_conjunction';
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.weights = fIdentityConne;
matlabbatch{batch_idx}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{batch_idx}.spm.stats.con.delete = 0;

%%
cd(resultsFolder);

%% display spm batch before running it
% spm_jobman('interactive',matlabbatch);
spm_jobman('run',matlabbatch);