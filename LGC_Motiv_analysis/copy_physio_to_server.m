function[]=copy_physio_to_server()
% script aimed at moving physiological data from pc to server to be used by
% future interns.

%% general parameters
study_nm = 'study1';
condition = 'behavior';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
dataRootPath = fullfile('E:','study1');
targetPath = fullfile('M:','Nicolas_Clairis','physio_extracted');
pulseTargetPath = [fullfile(targetPath, 'pulse'),filesep];
respiTargetPath = [fullfile(targetPath, 'respiration'),filesep];
    
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    % source folder
    subPath = [dataRootPath, filesep, 'CID',sub_nm,filesep,...
        'physiologs',filesep];
    cd(subPath);
    % target folders
    subPulseFolder = [pulseTargetPath,'CID',sub_nm];
    subRespiFolder = [respiTargetPath,'CID',sub_nm];
    mkdir(subPulseFolder);
    mkdir(subRespiFolder);
    % copy pulse files
    pulseFiles = ls(['CID',sub_nm,'*.puls']);
    lPulse = size(pulseFiles,1);
    for iFileP = 1:lPulse
        copyfile(pulseFiles(iFileP,:),subPulseFolder);
    end % pulse file loop
    % copy respiration files
    respiFiles = ls(['CID',sub_nm,'*.resp']);
    lRespi = size(respiFiles,1);
    for iFileR = 1:lRespi
        copyfile(respiFiles(iFileR,:),subRespiFolder);
    end % respiration file loop
    disp(['subject ',sub_nm,' - ',...
        num2str(iS),'/',num2str(NS),' - done']);
end % subject loop

end % function