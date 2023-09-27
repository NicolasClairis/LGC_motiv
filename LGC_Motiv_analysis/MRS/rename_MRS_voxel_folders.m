function[] = rename_MRS_voxel_folders()
% rename_MRS_voxel_folders will rename individual folders which were called
% "voxel_nativeSpace" into "MRS_voxels". The old name was framed in the
% idea that native-space and MNI-space data would be in different folders
% but it's actually easier to keep them all in the same folder, hence the
% folder name change.

%% subject selection
condition = 'fullList';
study_nm = 'study1';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);
%% working dir
pcStoragePath = fullfile('E:',study_nm);

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_fullNm = ['CID',sub_nm];
    sub_MRS_folder = fullfile(pcStoragePath, sub_fullNm, 'MRS');
    original_folder_name = fullfile(sub_MRS_folder,'voxel_nativeSpace');
    new_folder_name = fullfile(sub_MRS_folder,'MRS_voxels');
    movefile(original_folder_name, new_folder_name);
end % subject loop

end % function