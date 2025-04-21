function[] = Fcontrasts_DCM(GLM, checking, condition, study_nm, subject_id, NS, biasFieldCorr)
% [] = LGCM_contrasts_spm(GLM, checking, condition, study_nm, subject_id, NS, biasFieldCorr)
% LGCM_contrasts_spm prepares the batch for SPM to perform the first level
% contrasts in each individual of the study
%
% INPUTS
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
%
% GLM: GLM number
%
% checking: display batch before performing it or not? (1 by default)
%
% condition: condition to select subjects and sessions (will be asked if
% not defined by default)
%
% subject_id: list of subject (determined automatically if not defined in
% the inputs)
%
% NS: number of subjects (determined automatically if not defined in
% the inputs)
%
% biasFieldCorr: use bias-field corrected images (1) or not (0)? By default
% will not use bias-field corrected images
%

close all; clc;

%% working directories
% computer_root = LGCM_root_paths();
computer_root = ['E:',filesep];
if ~exist('study_nm','var') || isempty(study_nm)
    %     study_names = {'study1','study2','fMRI_pilots'};
    %     study_nm_idx = listdlg('ListString',study_names);
    %     study_nm = study_names{study_nm_idx};
    study_nm = 'study1'; % by default
end
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        root = fullfile(computer_root,'study1');
    case 'study2'
        root = fullfile(computer_root,'study2');
end

%% GLM
if ~exist('GLM','var') || isempty(GLM)
   GLM_info = inputdlg({'GLM number?'});
   GLM = str2double(GLM_info{1});
end

%% preprocessing smoothing kernel to consider
preproc_sm_kernel = 8;

%% use bias-field corrected files or not?
if ~exist('biasFieldCorr','var') || ~ismember(biasFieldCorr,[0,1])
    biasFieldCorr = 0;
end

%% which DCM mode to use? will define how sessions are grouped
[DCM_mode] = which_DCM_mode_for_GLM;

%% checking by default the batch before launching it
if ~exist('checking','var') ||...
        isempty(checking) ||...
        ~ismember(checking,[0,1])
    checking = 1;
end

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% define subjects of interest
if ~exist('condition','var') || isempty(condition) ||...
        ~strcmp(condition(1:4),'fMRI')
    condition = subject_condition;
end
gender = 'all';
if ~exist('subject_id','var') || ~exist('NS','var') ||...
        isempty(subject_id) || isempty(NS)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition, gender);
end

%% indication of GLM
disp(['GLM',num2str(GLM),' DCM_mode ',num2str(DCM_mode),' launching contrasts']);

%% loop through subjects to extract all the regressors
matlabbatch = cell(NS,1);

for iSubject = 1:NS
    %% extract current subject name
    sub_nm = subject_id{iSubject};
    
    %% extract contrasts list (vectors + corresponding names
    [fcon_names, fcon_vector] = Fcon_list_DCM(study_nm, sub_nm, GLM, DCM_mode, computer_root, preproc_sm_kernel, condition, biasFieldCorr);
    
    %% define results directory
    switch biasFieldCorr
        case 0
            mainPath = [fullfile(root,['CID',sub_nm],...
                'fMRI_analysis','functional',...
                ['preproc_sm_',num2str(preproc_sm_kernel),'mm_DCM']), filesep];
        case 1
            error('(Bias field = 1) AND (DCM model) not possible yet. Update scripts accordingly please.')
    end
    [resultsFolderName] = fMRI_subFolder_DCM(mainPath, GLM, condition, DCM_mode);
    matlabbatch{iSubject}.spm.stats.con.spmmat = {fullfile(resultsFolderName,'SPM.mat')};

    %% add each contrast to the list
    for iCon = 1:length(fcon_names)
        con_nm = ['con_',conNumber2conName(iCon)];
        matlabbatch{iSubject}.spm.stats.con.consess{iCon}.fcon.name     = fcon_names{iCon};
        matlabbatch{iSubject}.spm.stats.con.consess{iCon}.fcon.weights  = fcon_vector.(con_nm);
        matlabbatch{iSubject}.spm.stats.con.consess{iCon}.fcon.sessrep  = 'none'; % replicate contrast across sessions: no need as sessions are pooled
    end
    
    matlabbatch{iSubject}.spm.stats.con.delete = 1; % deletes previous contrasts if set to 1
end % subject loop

%% display spm batch before running it
if exist('checking','var') && checking == 1
    spm_jobman('interactive',matlabbatch);
%     spm_jobman('run',matlabbatch);
elseif ~exist('checking','var') || isempty(checking) || checking == 0
    spm_jobman('run',matlabbatch);
end

end % function