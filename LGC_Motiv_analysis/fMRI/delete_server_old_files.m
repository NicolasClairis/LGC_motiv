% Ugly script to brutally remove any unwanted GLM done in the past in order
% to clear a bit the server.

%% subject selection
study_nm = 'study1';
serverPath = fullfile('M:','human_data_private','analyzed_data',study_nm);
[subject_id, NS] = LGCM_subject_selection(study_nm,'fullList');

%% launch analysis
for iS = 1:NS
    sub_fullNm = ['CID',subject_id{iS}];
    sub_fMRI_folder = fullfile(serverPath, sub_fullNm, 'fMRI_analysis','functional');
%     %% 4mm smoothing kernel preprocessing
%     low_preproc_4mm_dir = [sub_fMRI_folder,filesep,'preproc_sm_4mm'];
%     if exist(low_preproc_4mm_dir,'dir')
%         rmdir(low_preproc_4mm_dir,'s');
%         disp([sub_fullNm,' preproc 4mm - removed']);
%     end
    
%     %% 6mm smoothing kernel preprocessing
%     low_preproc_6mm_dir = [sub_fMRI_folder,filesep,'preproc_sm_6mm'];
%     if exist(low_preproc_6mm_dir,'dir')
%         rmdir(low_preproc_6mm_dir,'s');
%         disp([sub_fullNm,' preproc 6mm - removed']);
%     end
    
    %% old GLM
    preproc_8mm_dir = [sub_fMRI_folder,filesep,'preproc_sm_8mm'];
    if exist(preproc_8mm_dir,'dir')
        for iGLM = 1:100
            GLM_to_clean = [preproc_8mm_dir,filesep,'GLM',num2str(iGLM)];
            if exist(GLM_to_clean,'dir')
                rmdir(GLM_to_clean,'s');
                disp([sub_fullNm,' GLM',num2str(iGLM),' - removed']);
            end
        end % GLM loop
    end
end % subject loop