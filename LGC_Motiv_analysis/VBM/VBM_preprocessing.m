

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% checking = 0 if you want to run directly or = 1 if you want to look at the batch first
checking = 1;

%% working directory
computer_root = ['E:',filesep];
studyRoot = fullfile(computer_root,study_nm);

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

anaPath = cell(NS,1);
% load anatomical scans
for iS = 1:NS
    sub_anat_folder = [studyRoot,filesep,'CID',subject_id{iS},filesep,...
        'fMRI_analysis',filesep,'anatomical',filesep];
    wms_anat_name = ls([sub_anat_folder,'CID*_UNI-DEN.mat']);
    anaPath(iS) = {[sub_anat_folder, wms_anat_name]};
end

unil_TPM_priorsPath = fullfile('C:','Users','clairis','Desktop','GitHub',...
    'LGC_motiv','Matlab_DIY_functions','VBM','enhanced_TPM.nii');

%% initialize size of matlabbatch
nBatchs = 3;
matlabbatch = cell(1,nBatchs);
nativeOnly = [0 0];
nativeAndDartel = [1 1];

%% segmentation
jSeg = 1;
matlabbatch{jSeg}.spm.spatial.preproc.channel.vols = anaPath;
matlabbatch{jSeg}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{jSeg}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{jSeg}.spm.spatial.preproc.channel.write = [0 0];
for iTissue = 1:6
    matlabbatch{jSeg}.spm.spatial.preproc.tissue(iTissue).tpm = {[unil_TPM_priorsPath,',',num2str(iTissue)]};
    % gaussian intensity distribution differs depending on tissue type
    switch iTissue
        case {1,2}
            matlabbatch{jSeg}.spm.spatial.preproc.tissue(iTissue).ngaus = 1;
        case {3,6}
            matlabbatch{jSeg}.spm.spatial.preproc.tissue(iTissue).ngaus = 2;
        case 4
            matlabbatch{jSeg}.spm.spatial.preproc.tissue(iTissue).ngaus = 3;
        case 5
            matlabbatch{jSeg}.spm.spatial.preproc.tissue(iTissue).ngaus = 4;
    end
    % extract DARTEL only for 3 first tissues
    switch iTissue
        case {1,2,3}
            matlabbatch{jSeg}.spm.spatial.preproc.tissue(iTissue).native = nativeAndDartel;
        case {4,5,6}
            matlabbatch{jSeg}.spm.spatial.preproc.tissue(iTissue).native = nativeOnly;
    end
    matlabbatch{jSeg}.spm.spatial.preproc.tissue(iTissue).warped = [0 0];
end
matlabbatch{jSeg}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{jSeg}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{jSeg}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{jSeg}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{jSeg}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{jSeg}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{jSeg}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{jSeg}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{jSeg}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
    NaN NaN NaN];
%% geodesic shooting
jShoot = 2;
matlabbatch{jShoot}.spm.tools.shoot.warp.images{1}(1) = cfg_dep('Segment: rc1 Images', substruct('.','val', '{}',{jSeg},...
    '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','rc', '()',{':'}));
matlabbatch{jShoot}.spm.tools.shoot.warp.images{2}(1) = cfg_dep('Segment: rc2 Images', substruct('.','val', '{}',{jSeg},...
    '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','rc', '()',{':'}));
%% normalise in MNI space
jMNI = 3;
matlabbatch{jMNI}.spm.tools.shoot.norm.template(1) = cfg_dep('Run Shooting (create Templates): Template (4)', substruct('.','val', '{}',{2},...
    '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','template', '()',{5}));
matlabbatch{jMNI}.spm.tools.shoot.norm.data.subjs.deformations(1) = cfg_dep('Run Shooting (create Templates): Deformation Fields', substruct('.','val', '{}',{2},...
    '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','def', '()',{':'}));
matlabbatch{jMNI}.spm.tools.shoot.norm.data.subjs.images{1}(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{jSeg},...
    '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{jMNI}.spm.tools.shoot.norm.data.subjs.images{2}(1) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{jSeg},...
    '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
%     matlabbatch{jMNI}.spm.tools.shoot.norm.vox = [NaN NaN NaN];
matlabbatch{jMNI}.spm.tools.shoot.norm.bb = [NaN NaN NaN
    NaN NaN NaN];
matlabbatch{jMNI}.spm.tools.shoot.norm.preserve = 1;
matlabbatch{jMNI}.spm.tools.shoot.norm.fwhm = 6;
matlabbatch{jMNI}.spm.tools.shoot.norm.vox = [0.6 0.6 0.6];


%% display spm batch before running it or run it directly
if checking == 1
    spm_jobman('interactive',matlabbatch);
elseif checking == 0
    spm_jobman('run',matlabbatch);
end