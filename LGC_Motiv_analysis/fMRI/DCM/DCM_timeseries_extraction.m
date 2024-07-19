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
matlabbatch = cell(NS*n_VOIs_per_sub,1);
matlabbatch = cell(1,1);
for iS = 1%:NS
    sub_nm = subject_id{iS};
    subFolder = [fullfile('E:',study_nm,['CID',sub_nm],...
        'fMRI_analysis','functional',...
        ['preproc_sm_',num2str(preproc_sm_kernel),'mm_DCM'],...
        ['GLM',num2str(GLM),'_no_satRun_DCM_mode',num2str(DCM_mode)]),filesep];
    MRS_ROI_folder = fullfile('E:',...
        study_nm,['CID',sub_nm],...
        'MRS','MRS_voxels');
    
    %% loop over volumes of interest
    for iVOI = 1%:n_VOIs_per_sub
        VOI_nm = VOI_names{iVOI};
        batch_idx = iVOI + n_VOIs_per_sub.*(iS - 1);
        
        matlabbatch{batch_idx}.spm.util.voi.spmmat = {[subFolder,'SPM.mat']};
        matlabbatch{batch_idx}.spm.util.voi.adjust = Fcon_to_adjust; % index of F-contrast used to adjust data
        matlabbatch{batch_idx}.spm.util.voi.session = 1; % all sessions pooled for DCM
        matlabbatch{batch_idx}.spm.util.voi.name = VOI_nm;
        switch VOI_nm
            case 'dmPFCdACC'
                matlabbatch{batch_idx}.spm.util.voi.roi{1}.mask.image = {[fullfile(MRS_ROI_folder,...
                    'wdmpfc.nii'),',1']};
            case 'aIns'
                matlabbatch{batch_idx}.spm.util.voi.roi{1}.mask.image = {[fullfile(MRS_ROI_folder,...
                    'wai.nii'),',1']};
            otherwise
                error(['ROI = ',VOI_nm,' not ready yet']);
        end
        matlabbatch{batch_idx}.spm.util.voi.roi{1}.mask.threshold = 0.1;
        matlabbatch{batch_idx}.spm.util.voi.expression = 'i1';
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