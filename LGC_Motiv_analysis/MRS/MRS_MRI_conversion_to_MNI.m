function[] = MRS_MRI_conversion_to_MNI(condition, study_nm)
% MRS_MRI_conversion_to_MNI
%
% INPUTS
% condition: condition of subjects to include
%
% study_nm: study name

%% subject selection
% condition
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
% study name
if ~exist('study_nm','var')
    study_nm = 'study1';
end
% list of subjects
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
study_path = fullfile('E:',study_nm);
spmFolderPath = fullfile('C:','Users','clairis','Desktop');

%% launch conversion from anatomical file to MNI
%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');
spmTemplatePath = fullfile(spmFolderPath,'spm12','spm12','tpm','TPM.nii');

%% prepare batch
% count number of batch to run per subject (will depend on whether one or
% two MRI files to normalize)
n_steps = 0;
n_steps_perMRI = 2;
for iS1 = 1:NS
    sub_nm1 = subject_id{iS1};
    switch sub_nm1
        case {'021','056','088'} % one MRI for aINS and one for dmPFC
            n_steps = n_steps + n_steps_perMRI*2;
        otherwise
            n_steps = n_steps + n_steps_perMRI;
    end
end % subject loop
% initialize batch
matlabbatch = cell(n_steps,1);

%% loop through subjects
jStep = 0;
for iS2 = 1:NS
    sub_nm2 = subject_id{iS2};
    sub_fullNm = ['CID',sub_nm2];
    switch sub_nm2
        case {'021','056','088'} % one MRI for aINS and one for dmPFC
            n_MRI = 2;
            sub_MRI_folder_dmPFC = [fullfile(study_path, sub_fullNm,'MRS','MRI','dmPFC_MRI'),filesep];
            dmPFC_anat_file = ls([sub_MRI_folder_dmPFC,'CID*UNI-DEN.nii']);
            sub_MRI_folder_ai = [fullfile(study_path, sub_fullNm,'MRS','MRI','ai_MRI'),filesep];
            ai_anat_file = ls([sub_MRI_folder_ai,'CID*UNI-DEN.nii']);
        otherwise
            n_MRI = 1;
            sub_MRI_folder = [fullfile(study_path, sub_fullNm,'MRS','MRI'),filesep];
            anat_file = ls([sub_MRI_folder,'CID*UNI-DEN.nii']);
    end % subject filter (for those who have 2 MRI)
    
    for iMRI = 1:n_MRI
        %% segment MRI
        jStep = jStep + 1;
        switch sub_nm2
            case {'021','056','088'} % one MRI for aINS and one for dmPFC
                switch iMRI
                    case 1 % dmPFC
                        matlabbatch{jStep}.spm.spatial.preproc.channel.vols = {[sub_MRI_folder_dmPFC,dmPFC_anat_file]};
                    case 2 % anterior insula
                        matlabbatch{jStep}.spm.spatial.preproc.channel.vols = {[sub_MRI_folder_ai,ai_anat_file]};
                end
            otherwise
                matlabbatch{jStep}.spm.spatial.preproc.channel.vols = {[sub_MRI_folder,anat_file]};
        end
        matlabbatch{jStep}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{jStep}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{jStep}.spm.spatial.preproc.channel.write = [1 1];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(1).tpm = {[spmTemplatePath,',1']};
        matlabbatch{jStep}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{jStep}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(1).warped = [1 1]; % normalise grey matter and record modulated and unmodulated warped tissue
        matlabbatch{jStep}.spm.spatial.preproc.tissue(2).tpm = {[spmTemplatePath,',2']};
        matlabbatch{jStep}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{jStep}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(2).warped = [1 1]; % normalise white matter and record modulated and unmodulated warped tissue
        matlabbatch{jStep}.spm.spatial.preproc.tissue(3).tpm = {[spmTemplatePath,',3']};
        matlabbatch{jStep}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{jStep}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(4).tpm = {[spmTemplatePath,',4']};
        matlabbatch{jStep}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{jStep}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(5).tpm = {[spmTemplatePath,',5']};
        matlabbatch{jStep}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{jStep}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(6).tpm = {[spmTemplatePath,',6']};
        matlabbatch{jStep}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{jStep}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{jStep}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{jStep}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{jStep}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{jStep}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{jStep}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{jStep}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{jStep}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{jStep}.spm.spatial.preproc.warp.write = [1 1]; % record inverse and forward deformation fields
        
        %% normalise MRI to MNI space
        jStep = jStep + 1;
        matlabbatch{jStep}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations',...
            substruct('.','val', '{}',{jStep-1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','fordef', '()',{':'}));
        matlabbatch{jStep}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Segment: Bias Corrected (1)',...
            substruct('.','val', '{}',{jStep-1}, '.','val', '{}',{1}, '.','val', '{}',{1}),...
            substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
        matlabbatch{jStep}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{jStep}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
        matlabbatch{jStep}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{jStep}.spm.spatial.normalise.write.woptions.prefix = 'w';
    end % MRI loop
end % subject loop

%% launch batch
spm_jobman('run',matlabbatch);

end % function