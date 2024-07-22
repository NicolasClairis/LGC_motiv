function[] = MRS_voxel_rmv_NaNs()
% MRS_voxel_rmv_NaNs will change NaN voxels inside MRS voxels into zeros 
% to allow timeseries extraction (which is crashing as long as the MRS
% voxels contain some NaN values).

%% subject selection
% condition
condition = 'fMRI';
% study name
if ~exist('study_nm','var')
    study_nm = 'study1';
end
% list of subjects
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
study_path = fullfile('E:',study_nm);

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% main parameters
MRS_ROIs = {'dmPFC','ai'};
nROIs = length(MRS_ROIs);
% initialize batch
matlabbatch = cell(1,1);

%% subject loop
jBatch_idx = 0; % batch index
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_fullNm = ['CID',sub_nm];
    
    for iROI = 1:nROIs
        MRS_ROI_nm = MRS_ROIs{iROI};
        sub_ROI_folder = [fullfile(study_path, sub_fullNm,'MRS','MRS_voxels'),filesep];
        ROI_file_nm = [sub_ROI_folder, 'bw', MRS_ROI_nm];
        
        %% filter case when ROI not extracted (because signal too bad or else)
        if exist(ROI_file_nm,'file')
            jBatch_idx = jBatch_idx + 1;
            
            matlabbatch{jBatch_idx}.spm.util.imcalc.input = {[ROI_file_nm,',1']};
            matlabbatch{jBatch_idx}.spm.util.imcalc.output = ['bw', MRS_ROI_nm,'_noNaN'];
            matlabbatch{jBatch_idx}.spm.util.imcalc.outdir = {sub_ROI_folder};
            matlabbatch{jBatch_idx}.spm.util.imcalc.expression = 'i1 + 0.*isnan(i1)'; % take initial image + change NaNs into zeros
            matlabbatch{jBatch_idx}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{jBatch_idx}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{jBatch_idx}.spm.util.imcalc.options.mask = 0;
            matlabbatch{jBatch_idx}.spm.util.imcalc.options.interp = 1;
            matlabbatch{jBatch_idx}.spm.util.imcalc.options.dtype = 4;
            
        end % ROI existing filter
    end % ROI loop
end % subject loop

%% launch batch
spm_jobman('run',matlabbatch);

end % function