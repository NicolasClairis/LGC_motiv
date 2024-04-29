% copy_MRS_anat_from_server will copy MRS original anatomical data + the
% anterior insula and dmPFC mask in the native space of the subject from
% the server to the pc.


%% subject selection
% condition
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
% study name
study_nm = 'study1';
% list of subjects
if ~exist('subject_id','var') || ~exist('NS','var') ||...
        isempty(subject_id) || isempty(NS)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
end

%% working directories
serverPath = fullfile('M:','human_data_private','analyzed_data',study_nm);
pcStoragePath = fullfile('E:',study_nm);

%% loop through subjects to obtain data
[anat_missing, aINS_missing, dmPFC_missing] = deal({});
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_fullNm = ['CID',sub_nm];
    
    %% source folder and files
    switch sub_nm
        case {'005','019','085'} % first visit can be ignored for those subjects
            MRS_folder = 'MRS_v2';
        otherwise
            MRS_folder = 'MRS';
    end
    subServerFolder = [fullfile(serverPath,...
        sub_fullNm,MRS_folder,[sub_fullNm,'_LGCMot']),filesep];
    aINS_file = [subServerFolder, 'ai',filesep,'ai.nii'];
    dmPFC_file = [subServerFolder, 'dmpfc',filesep,'dmpfc.nii'];
    % anatomical file
    switch sub_nm
        case {'021','056','088'}% one MRI for ai and one MRI for dmPFC for subject 021 => need special
            % treatment
            switch sub_nm
                case '021'
                    anat_folder_nm_dmPFC = '1_005_mp2rage_FN600b_FatNav_1mm_UNI-DEN_20220519';
                    anat_folder_nm_ai = '1_028_mp2rage_FN600b_FatNav_1mm_UNI-DEN_20220519';
                case '056'
                    anat_folder_nm_dmPFC = '1_005_mp2rage_FN600b_FatNav_1mm_UNI-DEN_20220218';
                    anat_folder_nm_ai = '1_026_mp2rage_FN600b_FatNav_1mm_UNI-DEN_20220218';
                case '088'
                    anat_folder_nm_dmPFC = '1_005_mp2rage_FN600b_FatNav_1mm_UNI-DEN_20220512';
                    anat_folder_nm_ai = '3_005_mp2rage_FN600b_FatNav_1mm_UNI-DEN_20220516'; % acquired on a different date
            end
            anat_folder_dmPFC = [subServerFolder, anat_folder_nm_dmPFC, filesep];
            anat_folder_ai = [subServerFolder, anat_folder_nm_ai, filesep];
            anat_file_nm_dmPFC = ls([anat_folder_dmPFC, 'CID*UNI-DEN.nii']);
            anat_file_nm_ai = ls([anat_folder_ai, 'CID*UNI-DEN.nii']);
            anat_file_dmPFC = [anat_folder_dmPFC, anat_file_nm_dmPFC];
            anat_file_ai = [anat_folder_ai, anat_file_nm_ai];
        otherwise
            anat_folder_nm = ls([subServerFolder, '*mp2rage_FN600b_FatNav_1mm_UNI-DEN*']);
            anat_folder = [subServerFolder, anat_folder_nm, filesep];
            anat_file_nm = ls([anat_folder, 'CID*UNI-DEN.nii']);
            anat_file = [anat_folder, anat_file_nm];
    end
    
    %% target folder and files
    sub_MRS_folder = fullfile(pcStoragePath, sub_fullNm, 'MRS');
    if ~exist(sub_MRS_folder,'dir')
        mkdir(sub_MRS_folder);
    end
    sub_MRS_MRI_folder = fullfile(sub_MRS_folder,'MRI');
    if ~exist(sub_MRS_MRI_folder,'dir')
        mkdir(sub_MRS_MRI_folder);
    end
    sub_MRS_voxels_folder = fullfile(sub_MRS_folder,'MRS_voxels');
    if ~exist(sub_MRS_voxels_folder,'dir')
        mkdir(sub_MRS_voxels_folder);
    end
    
    %% copy data
    % anatomical file
    switch sub_nm
        case {'021','056','088'} % different anatomy for ai and dmPFC for those subjects
            % dmPFC
            sub_MRS_voxels_folder_dmPFC = [sub_MRS_MRI_folder,filesep,'dmPFC_MRI'];
            if ~exist(sub_MRS_voxels_folder_dmPFC,'dir')
                mkdir(sub_MRS_voxels_folder_dmPFC);
            end
            copyfile(anat_file_dmPFC, sub_MRS_voxels_folder_dmPFC);
            
            % anterior insula
            sub_MRS_voxels_folder_ai = [sub_MRS_MRI_folder,filesep,'ai_MRI'];
            if ~exist(sub_MRS_voxels_folder_ai,'dir')
                mkdir(sub_MRS_voxels_folder_ai);
            end
            copyfile(anat_file_ai, sub_MRS_voxels_folder_ai);
        otherwise
            if exist(anat_file,'file')
                copyfile(anat_file, sub_MRS_MRI_folder);
            else
                disp([sub_fullNm,' anat file not found']);
                anat_missing = [anat_missing, sub_nm];
            end
    end
    
    % anterior insula file
    if exist(aINS_file,'file')
        copyfile(aINS_file, sub_MRS_voxels_folder);
    else
        disp([sub_fullNm,' aINS file not found']);
        aINS_missing = [aINS_missing, sub_nm];
    end
    
    % dmPFC file
    if exist(dmPFC_file,'file')
        copyfile(dmPFC_file, sub_MRS_voxels_folder);
    else
        disp([sub_fullNm,' dmPFC file not found']);
        dmPFC_missing = [dmPFC_missing, sub_nm];
    end
    
    %% message when worked fine
    if exist(anat_file,'file') && exist(aINS_file,'file') && exist(dmPFC_file,'file')
        disp([sub_fullNm,' copied']);
    end
    
end % subject loop