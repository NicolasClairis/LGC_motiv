function[badSubList] = find_bad_subs_for_GLM(GLM, condition)
%[badSubList] = find_bad_subs_for_GLM(GLM, condition)
% find_bad_subs_for_GLM aims at identifying the subjects for which the
% first level and/or the contrasts failed for a given GLM (specified in the
% inputs).
%
% INPUTS
% GLM: GLM number
%
% condition: define subjects and runs to include
% 'fMRI': all subjects where fMRI ok
% 'fMRI_no_move': remove runs with too much movement
%
% OUTPUTS
% badSubList: structure with names of the subjects where the GLM was not
% performed either for first level and/or for contrasts
%
% Written by Nicolas Clairis - april 2022

%% define study
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% define working directories
% computer_root = LGCM_root_paths();
computerRoot = ['E:',filesep];
switch study_nm
    case 'fMRI_pilots'
        studyRoot = fullfile(computerRoot,'fMRI_pilots');
    case 'study1'
        studyRoot = fullfile(computerRoot,'study1');
    case 'study2'
        studyRoot = fullfile(computerRoot,'study2');
end

%% define preprocessing smoothing kernel
preproc_sm_kernel = 8;
if ~exist('preproc_sm_kernel','var') || isempty(preproc_sm_kernel)
   preproc_sm_kernel_cell = inputdlg('What smoothing kernel to use for preprocessing?');
   preproc_sm_kernel = str2double(preproc_sm_kernel_cell{1});
end

%% extract GLM if not defined yet
if ~exist('GLM','var') || isempty(GLM)
   GLM_cell = inputdlg('GLM number?');
   GLM = str2num(GLM_cell{1});
end
GLM_folder = [filesep, 'fMRI_analysis',filesep,'functional',filesep,...
    'preproc_sm_',num2str(preproc_sm_kernel),'mm',filesep];
GLM_path = fMRI_subFolder(GLM_folder,GLM,condition);
% extract regressors
[~, n_regsPerTask] = GLM_details(GLM);

%% get full subject list
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

% prepare the vectors of participants where GLM failed
[badFirstLevel_idx, badContrast_idx] = deal(NaN(1,NS));

%% loop on all subjects and check who does not have a first level and/or a constrast
for iS = 1:NS
    sub_nm = subject_id{iS};
    subFolder = [studyRoot, filesep, 'CID',sub_nm, GLM_path];
    % get the theoretical number of runs
    runs = runs_definition(study_nm, sub_nm, condition);
    
    if ~exist(subFolder,'dir') % no first level performed for this subject
        badFirstLevel_idx(iS) = true;
        badContrast_idx(iS) = true;
    else % first level folder exists => is there data? were contrasts performed?
        % get the theoretical number of contrasts
        con_names = LGCM_contrasts(study_nm, sub_nm, GLM,...
            computerRoot, preproc_sm_kernel, condition);
        n_cons = size(con_names,2);
        
        % get number of regressors/task*number of runs/task + 1 constant/run
        n_totalRegs = n_regsPerTask.Ep*(runs.nb_runs.Ep) +...
            n_regsPerTask.Em*(runs.nb_runs.Em) +...
            runs.nb_runs.Ep + runs.nb_runs.Em;
        
        
        betaFiles = ls([subFolder,'beta_*.nii']);
        n_1stLevelBetas = size(betaFiles,1);
        if ~exist([subFolder,'SPM.mat'],'file') || n_1stLevelBetas < n_totalRegs
            badFirstLevel_idx(iS) = true;
            badContrast_idx(iS) = true;
        else % 1st level seemed to work but contrasts could have failed
            badFirstLevel_idx(iS) = false;
            conFiles = ls([subFolder,'spmT_*.nii']);
            n_conFiles = size(conFiles,1);
            if n_conFiles < n_cons
                badContrast_idx(iS) = true;
            else
                badContrast_idx(iS) = false;
            end % contrasts
        end % 1st level
    end % check for first level
end % subject loop

%% extract name of bad subjects
[badSubList.firstLevel, badSubList.contrast] = deal({});
if sum(badFirstLevel_idx) > 0
    for iS = 1:NS
        sub_nm = subject_id{iS};
        if badFirstLevel_idx(iS) == true
            badSubList.firstLevel = [badSubList.firstLevel, sub_nm];
        end % filter bad first level
    end % subject loop
end % if there are some subs with bad first level

if sum(badContrast_idx) > 0
    for iS = 1:NS
        sub_nm = subject_id{iS};
        if badContrast_idx(iS) == true
            badSubList.contrast = [badSubList.contrast, sub_nm];
        end % filter bad contrast
    end % subject loop
end % if there are some subs with bad contrast

end % function