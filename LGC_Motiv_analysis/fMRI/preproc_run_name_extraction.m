function[filenames] = preproc_run_name_extraction(study_nm, sub_nm, runPath, prefix)
%[filenames] = preproc_run_name_extraction(study_nm, sub_nm, runPath, prefix)
% preproc_run_name_extraction will extract the list of files corresponding
% to the selected study, subject and run, and based on the prefix defined
% in 'prefix'.
%
% INPUTS
% study_nm: 'study1'/'study2'
%
% sub_nm: string with subject identification number 'XXX'
%
% runPath: path (with file separator at the end) indicating where you
% should check
%
% prefix: prefix indicating whether to look at the raw data, or to the
% re-aligned data ('r') or any other filter you may want to apply
%
% OUTPUT
% filenames: list of files corresponding to the run selected
%
% See also preprocessing_batch.m

%% by default, prefix empty
if ~exist('prefix','var') || isempty(prefix)
    prefix = '';
end

%% extract files
if strcmp(study_nm,'fMRI_pilots')
    if ismember(sub_nm,{'pilot_s1'})
        filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'LGCM_.*\.(img|nii)$']));
    elseif ismember(sub_nm,{'pilot_s2'})
        filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'run.*\.(img|nii)$']));
    elseif ismember(sub_nm,{'pilot_s3'})
        filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'ABNC.*\.(img|nii)$']));
    else
        filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'CID.*\.(img|nii)$']));
    end
elseif strcmp(study_nm,'study2_pilots')
    if ismember(sub_nm,{'fMRI_pilot1_AC'})
        filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'AC.*\.(img|nii)$']));
    else
        error('pilots other than fMRI_pilot1_AC not ready yet');
    end
else
    filenames = cellstr(spm_select('ExtFPList',runPath,['^',prefix,'CID.*\.(img|nii)$']));
    %             error('please check the format (nii/img) and the start of the name of each run because it has to be stabilized now...');
end

end % function