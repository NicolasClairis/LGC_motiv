% Script to copy first level data from local pc to server.

%% subject selection
study_nm = 'study1';
rootPath = fullfile('E:',study_nm);
serverPath = fullfile('M:','human_data_private','analyzed_data',study_nm);
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% GLM to copy
GLM_str = inputdlg('GLM to copy?');
GLM = str2double(GLM_str{1});

for iS = 1:NS
    sub_fullNm = ['CID',subject_id{iS}];
    % define source data
    subRootGlobalFolder = [rootPath, filesep, sub_fullNm, filesep,...
        'fMRI_analysis',filesep,'functional',filesep,'preproc_sm_8mm'];
    subRootFolder = fMRI_subFolder(subRootGlobalFolder, GLM, condition);
    
    % define target folder
    subTargetGlobalFolder = [serverPath, filesep, sub_fullNm, filesep,...
        'fMRI_analysis',filesep,'functional',filesep,'preproc_sm_8mm'];
    subTargetFolder = fMRI_subFolder(subTargetGlobalFolder, GLM, condition);
    
    if exist(subRootFolder,'dir') && ~exist(subTargetFolder,'dir')
        status = copyfile(subRootFolder, subTargetFolder);
        if status == 1
            disp(['Sub ',num2str(iS),'/',num2str(NS),' (',sub_fullNm,') - copied.']);
        end
    end
end % sub loop