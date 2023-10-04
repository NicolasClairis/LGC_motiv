function[] = delete_fMRI_raw()
% delete_fMRI_raw will delete the 'fMRI_raw' folders within each subject.
% Since these folders only contain all the raw data before the conversion
% to make them usable nifti files with mriconvert.exe, they are just taking
% much space for nothing => delete_fMRI_raw will delete them to free some
% space from your hard drive.

%% working directory
rootPath = LGCM_root_paths;
study_nm = 'study1';
study_path = [fullfile(rootPath, study_nm), filesep];
%% subject selection
condition = 'fullList';
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% loop through subjects
for iS = 1:NS
    sub_CID_nm = ['CID',subject_id{iS}];
    % extract folder name
    subFolder = [study_path, sub_CID_nm, filesep,'fMRI_raw'];
    %% remove folder (if it exists)
    if exist(subFolder,'dir')
       rmdir(subFolder,'s');
       disp([sub_CID_nm,' - ',num2str(iS),'/',num2str(NS),' removed']);
    else
        disp([sub_CID_nm,' no fMRI_raw folder found']);
    end
end % subject loop

end % function