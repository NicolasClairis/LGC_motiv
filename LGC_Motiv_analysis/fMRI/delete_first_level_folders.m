function[] = delete_first_level_folders(study_nm)
% function to order first level files coming from one main sequence in a
% folder in order to make new analysis
% function by Nicolas Clairis - june 2016

baseline_folder = pwd;

kernels = {'8','6','4'};
sm_kernel_idx = listdlg('PromptString','What smoothing kernel do you want to delete?',...
    'ListString',kernels);
sm_kernel = str2double(kernels{sm_kernel_idx});

%% GLM choice
GLM_str = inputdlg('What First level GLM do you want to delete?');
GLM = str2double(GLM_str);

if ~exist('study_nm','var') || isempty(study_nm)
    study_names = {'study1','study2','fMRI_pilots'};
    study_nm_idx = listdlg('ListString',study_names);
    study_nm = study_names{study_nm_idx};
end


% check if sure to delete before doing so
del_check = listdlg('PromptString',...
    ['Are you sure you want to delete the first level folders of GLM ',num2str(GLM),' for',study_nm,'?'],...
    'ListString',{'yes','no'},'SelectionMode','single');
if del_check == 2
    disp('You didn''t want to delete these files finally huh?');
    return;
end

%% define subject list
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

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

% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subFullNm = ['CID',sub_nm];
    subj_folder = [root, filesep, subFullNm];
    subj_analysis_folder = [subj_folder,filesep,'fMRI_analysis',filesep,...
        'functional',filesep,'preproc_sm_',num2str(sm_kernel),'mm', filesep];
    if exist(subj_analysis_folder,'dir')
        cd(subj_analysis_folder);
        
        [resultsFolderName, resultsFolderShortName] = fMRI_subFolder(subj_analysis_folder, GLM, condition);
        if exist(resultsFolderName,'dir')
            rmdir(resultsFolderShortName,'s');
            disp([resultsFolderShortName, ' correctly removed for ',sub_nm]);
        else
            disp([resultsFolderShortName, ' not found for ',sub_nm]);
        end
    else
        disp([subj_analysis_folder, ' not found for ',sub_nm]);
    end
end

%% go back to initial folder
if exist(baseline_folder,'dir') % folder might have been removed by delete_first_level_folders.m
    cd(baseline_folder);
end

end % function