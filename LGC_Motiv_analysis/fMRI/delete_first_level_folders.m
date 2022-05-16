function[] = delete_first_level_folders(study_nm)
% function to order first level files coming from one main sequence in a
% folder in order to make new analysis
% function by NicoC june 2016

baseline_folder = pwd;

sm_kernel = input(sprintf('What smoothing kernel do you want to delete (4/6/8)? \n'));

GLM = input(sprintf('What First level GLM do you want to delete? \n'));

if ~exist('study_nm','var') || isempty(study_nm)
    study_names = {'fMRI_pilots','study1','study2'};
    study_nm_idx = listdlg('ListString',study_names);
    study_nm = study_names{study_nm_idx};
end


% check if sure to delete before doing so
del_check = input(sprintf(['Are you sure you want to delete the first level folders of GLM ',num2str(GLM),' for',study_nm,'?',...
    '\n If yes, press 1. If no, press whatever else. \n']));
if del_check ~= 1
    disp('You didn''t want to delete these files finally huh?');
    return;
end

condition = 'fMRI';
warning('should add a filter here to select fMRI/fMRI_noMove/etc.');

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
    subFullNm = ['CID',subject_id{iS}];
    subj_folder = [root, filesep, subFullNm];
    subj_analysis_folder = [subj_folder,filesep,'fMRI_analysis',filesep,...
        'functional',filesep,'preproc_sm_',num2str(sm_kernel),'mm', filesep];
    cd(subj_analysis_folder);
    
    run_folders = {['GLM',num2str(GLM)]};
    for iRun = 1:length(run_folders)
        if exist([subj_analysis_folder, run_folders{iRun}],'dir')
            rmdir(run_folders{iRun},'s');
            disp([run_folders{iRun}, ' correctly removed for ',subject_id{iS}]);
        else
            disp([run_folders{iRun}, ' not found for ',subject_id{iS}]);
        end
    end
end

%% go back to initial folder
cd(baseline_folder);

end % function