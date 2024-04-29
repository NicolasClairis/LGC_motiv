% delete_MRS_data_from_pc will delete all the folders containing MRS data
% inside the pc.


%% list of subjects
% study name
study_nm = 'study1';
condition = 'fullList';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directory
rootDir = [fullfile('E:',study_nm),filesep];

%% loop through subjects
no_MRS_subs = [];
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_fullNm = ['CID',sub_nm];
    subFolder = [rootDir, sub_fullNm, filesep];
    MRS_folder = [subFolder, 'MRS',filesep];
    if exist(MRS_folder,'dir')
        rmdir(MRS_folder,'s');
        disp(['MRS folder removed for subject ',sub_nm]);
    else
        disp(['MRS folder missing for subject ',sub_nm]);
        no_MRS_subs = [no_MRS_subs, sub_nm];
    end
end % subject loop