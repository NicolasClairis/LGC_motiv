function[] = copy_behavior_to_server()
% script aimed at moving behavioral data from pc to server to be used by
% future interns.

%% general parameters
study_nm = 'study1';
condition = 'behavior';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
dataRootPath = fullfile('E:','study1');
targetPath = [fullfile('M:','Nicolas_Clairis','behavior_extracted'),...
    filesep];
    
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    % source folder
    subPath = [dataRootPath, filesep, 'CID',sub_nm,filesep,...
        'behavior',filesep];
    cd(subPath);
    % target folders
    subNewBehaviorFolder = [targetPath,'CID',sub_nm];
    mkdir(subNewBehaviorFolder);
    % copy files
    behavioralFiles = ls(['CID',sub_nm,'*.mat']);
    lFiles = size(behavioralFiles,1);
    for iFileP = 1:lFiles
        copyfile(behavioralFiles(iFileP,:),subNewBehaviorFolder);
    end % pulse file loop
    disp(['subject ',sub_nm,' - ',...
        num2str(iS),'/',num2str(NS),' - done']);
end % subject loop

end % function