function[] = preprocessing_batch(study_nm)
%[] = preprocessing_batch(study_nm)
% preprocessing for fMRI data
% enter subject identification in 'subject_id' (sXX_ddMMyy dd: day,
% MM: month, yy: year) and preprocessing number in 'preproc'
% preproc: 0: "raw" preprocessing (grey matter saved during segmentation)
%          1: added AP-PA correction
%
% all data is stored in the fMRI_scans folders + the script also creates
% fMRI_analysis folder where the anatomical files are copied for the
% preprocessing + functional folder created for 1st level analysis
%
% Preprocessing entails 1) realignement, 2) co-registration, 3)
% segmentation, 4-5) normalisation in MNI space, 6) spatial smoothing
%
% INPUTS
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
%
% See also First_level_batch, contrasts_batch and
% Second_level_batch

% clear;

% %% preprocessing number in case you want to have several versions of it
% preproc = 0; % see intro

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% define subjects and working directories
[computerRoot, spmFolderPath] = LGCM_root_paths();
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm_List = {'study1','study2','fMRI_pilots'};
    study_nm_idx = listdlg('ListString',study_nm_List);
    study_nm = study_nm_List{study_nm_idx};
end
switch study_nm
    case 'fMRI_pilots' % pilots
        root = [fullfile(computerRoot,'fMRI_pilots'),filesep];
    case 'study1'
        root = [fullfile(computerRoot,'study1'),filesep];
    case 'study2'
        root = [fullfile(computerRoot,'study2'),filesep];
end
% subject_id = {'pilot_s2'}; % 'pilot_s1','pilot_s2'
% NS = length(subject_id); % nber of subjects
[subject_id, NS] = LGCM_subject_selection(study_nm);

% give path for anatomical template
spmTemplatePath = fullfile(spmFolderPath,'spm12','spm12','tpm','TPM.nii');


%% define number of preprocessing steps
nb_preprocessingSteps = 6;
 % nb of preprocessing steps per subject:
 % 1)realign-reslice
 % 2)coregistration functional - anatomical
 % 3) segmentation
 % 4) normalization functional
 % 5) normalization anatomical
 % 6) smoothing functional
 
 %% initialize batch
 matlabbatch = cell(nb_preprocessingSteps*NS,1);

for iS = 1:NS % loop through subjects
    disp(['loading batch for preprocessing subject ',num2str(iS),'/',num2str(NS)]);
    sub_nm = subject_id{iS};
    suf_fullNm = ['CID',sub_nm];
    
    % create working directories and copy anat. file inside
    % \fMRI_analysis\anatomical folder
    cd([root, suf_fullNm]);
    if exist('fMRI_analysis','dir') ~= 7
        mkdir fMRI_analysis;
    end
    subj_analysis_folder = [root, suf_fullNm,filesep,'fMRI_analysis'];
    cd(subj_analysis_folder);
    if exist('anatomical','dir') ~= 7
        mkdir anatomical;
    end
    if exist('functional','dir') ~= 7
        mkdir functional;
    end
    subj_scans_folder = [root, suf_fullNm, filesep,'fMRI_scans'];
    cd(subj_scans_folder);
    anat_folder = ls('*UNI-DEN*');
    newAnatFolder = [subj_analysis_folder, filesep,'anatomical',filesep];
    copyfile(anat_folder, newAnatFolder); % copies the contempt of the anatomical folder
    
    %% extract folders where functional runs are stored
    cd(subj_scans_folder);
    subj_scan_folders_names = ls('*run*'); % takes all functional runs folders
    % remove AP/PA top-up corrective runs when they were performed (only 2
    % first pilots)
    if strcmp(study_nm,'fMRI_pilots') && ismember(sub_nm,{'pilot_s1','pilot_s2'})
        [subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names);
    end
    %%
    cd(subj_analysis_folder)
    %% define number of sessions to analyze
    if strcmp(study_nm,'fMRI_pilots') &&...
            ismember(sub_nm, {'pilot_s1','pilot_s2','pilot_s3'}) % only 2 sessions for these pilots
        nb_runs = 2;
    else
        nb_runs = 4;
    end
    cd(subj_scans_folder);
    for iRun = 1:nb_runs % loop through runs for 3 ratings, 3 choices 1D, 3 choices 2D runs
        cd(subj_scan_folders_names(iRun,:)); % go to run folder
        if strcmp(study_nm,'fMRI_pilots')
            if ismember(sub_nm,{'pilot_s1'})
                filenames = cellstr(spm_select('ExtFPList',pwd,'^LGCM_.*\.nii$'));
            elseif ismember(sub_nm,{'pilot_s2'})
                filenames = cellstr(spm_select('ExtFPList',pwd,'^run.*\.nii$'));
            elseif ismember(sub_nm,{'pilot_s3'})
                filenames = cellstr(spm_select('ExtFPList',pwd,'^ABNC.*\.img$'));
            else
                filenames = cellstr(spm_select('ExtFPList',pwd,'^CID.*\.nii$'));
            end
        else
            filenames = cellstr(spm_select('ExtFPList',pwd,'^CID.*\.nii$'));
            %             error('please check the format (nii/img) and the start of the name of each run because it has to be stabilized now...');
        end
        runFileNames.(['run_',num2str(iRun)]) = filenames;
        cd(subj_scans_folder);
    end
    
    %% realignement
    cd(subj_scans_folder);
    preproc_step = 1; % first step of preprocessing
    realign_step = nb_preprocessingSteps*(iS-1) + preproc_step; % number depends also on the number of subjects
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.data = cell(nb_runs, 1);
    for iRun = 1:nb_runs
        matlabbatch{realign_step}.spm.spatial.realign.estwrite.data{iRun,1} = runFileNames.(['run_',num2str(iRun)]);
    end
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{realign_step}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    %% coregistration
    preproc_step = 2;
    coreg_step = nb_preprocessingSteps*(iS-1) + preproc_step;
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    cd(newAnatFolder);
    if ismember(sub_nm,{'pilot_s1'})
        anat_file = ls('LGCM_*.nii');
    elseif ismember(sub_nm,{'pilot_s3'})
        anat_file = ls('ABNC_*.img');
    else
        %         anat_file = ls('mp2rage_*.nii');
        anat_file = ls('CID*.nii');
    end
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.source = {[newAnatFolder, anat_file]};
    cd(subj_scans_folder);
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{coreg_step}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    %% segmentation
    cd(subj_scans_folder);
    preproc_step = 3;
    segm_step = nb_preprocessingSteps*(iS-1) + preproc_step;
    matlabbatch{segm_step}.spm.spatial.preproc.channel.vols = {[newAnatFolder,anat_file]};
    matlabbatch{segm_step}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{segm_step}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{segm_step}.spm.spatial.preproc.channel.write = [1 1];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).tpm = {[spmTemplatePath,',1']};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(1).warped = [1 1]; % normalise grey matter and record modulated and unmodulated warped tissue
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).tpm = {[spmTemplatePath,',2']};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).tpm = {[spmTemplatePath,',3']};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).tpm = {[spmTemplatePath,',4']};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).tpm = {[spmTemplatePath,',5']};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).tpm = {[spmTemplatePath,',6']};
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{segm_step}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{segm_step}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{segm_step}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{segm_step}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{segm_step}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{segm_step}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{segm_step}.spm.spatial.preproc.warp.write = [1 1];
    %% normalization
    preproc_step = 4;
    normf_step = nb_preprocessingSteps*(iS-1) + preproc_step;
    matlabbatch{normf_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    for iRun = 1:nb_runs
        matlabbatch{normf_step}.spm.spatial.normalise.write.subj.resample(iRun) = cfg_dep(['Realign: Estimate & Reslice: Resliced Images (Sess ',num2str(iRun),')'], substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{iRun}, '.','rfiles'));
    end
    matlabbatch{normf_step}.spm.spatial.normalise.write.subj.resample(nb_runs+1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{realign_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{normf_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
    %% normalization anatomical scan
    preproc_step = 5;
    norma_step = nb_preprocessingSteps*(iS-1) + preproc_step;
    matlabbatch{norma_step}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{norma_step}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{segm_step}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{norma_step}.spm.spatial.normalise.write.woptions.prefix = 'w';
    %% smoothing
    preproc_step = 6;
    smooth_step = nb_preprocessingSteps*(iS-1) + preproc_step;
    matlabbatch{smooth_step}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{normf_step}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{smooth_step}.spm.spatial.smooth.fwhm = [4 4 4];
    matlabbatch{smooth_step}.spm.spatial.smooth.dtype = 0;
    matlabbatch{smooth_step}.spm.spatial.smooth.im = 0;
    matlabbatch{smooth_step}.spm.spatial.smooth.prefix = 's';

    cd(root);
end

% display spm batch before running it
spm_jobman('interactive',matlabbatch);
% spm_jobman('run',mtatlabbatch);