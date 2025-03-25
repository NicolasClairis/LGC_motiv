% script to copy-paste LGC Motivation project behavior + MRS files for
% Arthur to work on it while ignoring heavy fMRI files

%% working directories
source_path = [filesep,filesep,fullfile('sv-nas1.rcp.epfl.ch',...
    'Sandi-lab','human_data_private')];
harddrive_path = input('nom du disque dur? (Format ''X:/'' svp');
target_path = fullfile(harddrive_path,'LGCMot');
% create target path if not already created
if ~exist(target_path,'dir')
    mkdir(target_path);
end

%% list of subjects
subjects = {'CID001','CID002','CID003','CID004','CID005','CID008','CID009',...
            'CID011','CID012','CID013','CID015','CID017','CID018','CID019',...
            'CID020','CID021','CID022','CID024','CID027','CID029',...
            'CID030','CID032','CID034','CID035','CID036','CID038','CID039',...
            'CID040','CID042','CID043','CID044','CID045','CID046','CID047','CID048','CID049',...
            'CID050','CID052','CID053','CID054','CID055','CID056','CID058','CID059',...
            'CID060','CID061','CID062','CID064','CID065','CID068','CID069',...
            'CID071','CID072','CID073','CID074','CID075','CID076','CID078','CID079',...
            'CID080','CID081','CID082','CID083','CID085','CID086','CID087','CID088',...
            'CID090','CID091','CID093','CID094','CID095','CID097','CID099','CID100'};

%% copy-paste all the data
for iS = 1:NS 
    sub_nm = subjects{iS};
    subj_targetPath = fullfile(target_path,sub_nm);
    % create folder with current subject path in the hard drive
    if ~exist(subj_targetPath,'dir')
        mkdir(subj_targetPath);
    end

    % copy each folder
    %% behavior
    bhv_sourcePath = fullfile(source_path,'raw_data_subject',...
        'study1',sub_nm,'behavior');
    [bhv_copy_success] = copyfile(bhv_sourcePath,subj_targetPath,'f');
    
    %% MRS data
    % raw data
    MRS_raw_sourcePath = fullfile(source_path,'raw_data_subject',...
        'study1',sub_nm,'MRS');
    MRS_raw_data_targetPath = fullfile(subj_targetPath,'MRS','raw_data');
    if ~exist(MRS_raw_data_targetPath,'dir')
        mkdir(MRS_raw_data_targetPath);
    end
    [MRS_raw_copy_success] = copyfile(MRS_raw_sourcePath,MRS_raw_data_targetPath,'f');
    %% analyzed MRS data
    MRS_analyzed_sourcePath = fullfile(source_path,'analyzed_data',...
        'study1',sub_nm,'MRS');
    MRS_anal_data_targetPath = fullfile(subj_targetPath,'MRS','raw_data');
    if ~exist(MRS_anal_data_targetPath,'dir')
        mkdir(MRS_anal_data_targetPath);
    end
    [MRS_anal_copy_success] = copyfile(MRS_analyzed_sourcePath,MRS_anal_data_targetPath,'f');

    %% check where we are and that all went well
    if bhv_copy_success && MRS_raw_copy_success && MRS_anal_copy_success
        disp(['Subject ',sub_nm,': ',num2str(iS),'/',num2str(NS),' copy done.']);
    else
        warning(['Subject ',sub_nm,': problem with copy-pasting.']);
    end
end % subject loop

%% copy-paste also the general folders
[summary_success] = copyfile(fullfile(source_path,'Summary'),target_path,'f');
if summary_success
    disp('Summary folder copy done.');
else
    warning('Summary folder: problem with copy-pasting.');
end