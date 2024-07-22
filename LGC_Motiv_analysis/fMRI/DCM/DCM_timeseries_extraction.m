function[] = DCM_timeseries_extraction(GLM, checking)
% DCM_timeseries_extraction(GLM, checking)
% DCM_timeseries_extraction serves to extract the timeseries in the ROI of
% interest in order to perform the DCM.
%
% INPUTS
% GLM: GLM number
%
% checking: display batch before launching it or not
% (0) launch batch directly
% (1) display batch before launching it

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% GLM
if ~exist('GLM','var') || isempty(GLM)
   GLM_info = inputdlg({'GLM number?'});
   GLM = str2double(GLM_info{1});
end

%% which DCM mode to use? will define how sessions are grouped
[DCM_mode] = which_DCM_mode_for_GLM;

%% value of the smoothing during preprocessing?
preproc_sm_kernel = 8;

%% index of F-contrast used to adjust data
% define the way to define the effects of interest
% (1) only ignores movement + run constant
% (2) ignores movement + run constant + fixation cross
% (3) ignores movement + run constant + fixation cross + dispChosen + feedback
Fcon_to_adjust = 2;

%% loop over subjects to load matlabbatch
VOI_names = {'dmPFCdACC','aIns'};
n_VOIs_per_sub = length(VOI_names);
%% initialize batch
matlabbatch = cell(1,1);
jBatch = 0;
%% loop over subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subFolder = [fullfile('E:',study_nm,['CID',sub_nm],...
        'fMRI_analysis','functional',...
        ['preproc_sm_',num2str(preproc_sm_kernel),'mm_DCM'],...
        ['GLM',num2str(GLM),'_no_satRun_DCM_mode',num2str(DCM_mode)]),filesep];
    MRS_ROI_folder = [fullfile('E:',...
        study_nm,['CID',sub_nm],...
        'MRS','MRS_voxels'),filesep];
    
    %% loop over volumes of interest
    for iVOI = 1:n_VOIs_per_sub
        VOI_nm = VOI_names{iVOI};
        % name of MNI-normalized and binarized VOI file
        switch VOI_nm
            case 'dmPFCdACC'
                VOI_file_nm = 'bwdmpfc';
            case 'aIns'
                VOI_file_nm = 'bwai';
        end
        ROI_file_nm = [MRS_ROI_folder, VOI_file_nm,'.nii'];
        
        if exist(ROI_file_nm,'file')
            jBatch = jBatch + 1;
            matlabbatch{jBatch}.spm.util.voi.spmmat = {[subFolder,'SPM.mat']};
            matlabbatch{jBatch}.spm.util.voi.adjust = Fcon_to_adjust; % index of F-contrast used to adjust data
            matlabbatch{jBatch}.spm.util.voi.session = 1; % all sessions pooled for DCM
            matlabbatch{jBatch}.spm.util.voi.name = [VOI_nm,'_timeseries'];
            % load binarized MRS voxels
            switch VOI_nm
                case {'dmPFCdACC','aIns'}
                    matlabbatch{jBatch}.spm.util.voi.roi{1}.mask.image = {[ROI_file_nm,',1']};
                otherwise
                    error(['ROI = ',VOI_nm,' not ready yet']);
            end
            matlabbatch{jBatch}.spm.util.voi.roi{1}.mask.threshold = 0.1;
            matlabbatch{jBatch}.spm.util.voi.roi{2}.mask.image = {[subFolder,'mask.nii,1']};
            matlabbatch{jBatch}.spm.util.voi.roi{2}.mask.threshold = 0.1;
            matlabbatch{jBatch}.spm.util.voi.expression = '(i1>0).*(i2>0)==1';
        end % filter if file name exists or not
    end % VOI loop
end % subject loop

%% display spm batch before running it
if exist('checking','var') && checking == 1
    spm_jobman('interactive',matlabbatch);
%     spm_jobman('run',matlabbatch);
elseif ~exist('checking','var') || isempty(checking) || checking == 0
    spm_jobman('run',matlabbatch);
end

end % function