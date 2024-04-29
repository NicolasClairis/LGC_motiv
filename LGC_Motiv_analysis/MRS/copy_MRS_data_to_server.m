% Script to copy MRS mask from local pc to server.

%% subject selection
study_nm = 'study1';
rootPath = fullfile('E:',study_nm);
serverPath = fullfile('M:','human_data_private','analyzed_data',study_nm);
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

for iS = 1:NS
    sub_fullNm = ['CID',subject_id{iS}];
    % define source data
    subRootFolder = [rootPath, filesep, sub_fullNm, filesep,...
        'MRS',filesep,'MRS_voxels'];
    
    % define target folder
    subTargetFolder = [serverPath, filesep, sub_fullNm, filesep,...
        'MRS',filesep,'MRS_voxels'];
    
    if exist(subRootFolder,'dir') && ~exist(subTargetFolder,'dir')
        status = copyfile(subRootFolder, subTargetFolder);
        if status == 1
            disp(['Sub ',num2str(iS),'/',num2str(NS),' (',sub_fullNm,') - copied.']);
        end
    end
end % sub loop