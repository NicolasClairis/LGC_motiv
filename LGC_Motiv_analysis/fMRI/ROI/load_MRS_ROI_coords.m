function[sxyz_ROI] = load_MRS_ROI_coords(study_nm, sub_nm, ROI_nm)
% [sxyz_ROI] = load_MRS_ROI_coords(study_nm, sub_nm, ROI_nm)
% load_MRS_ROI_coords will extract the coordinates of the MRS voxel
% selected in input in order to perform the extraction of the data in this
% voxel through the corresponding ROI_extraction_group.m function.
%
% INPUTS
% study_nm: study name ('study1' by default if left empty)
%
% sub_nm: subject name as string made of 3 numbers 'XXX'
%
% ROI_nm:
% 'MRS_dmPFC': dorsomedial prefrontal cortex voxel used for spectroscopy
% 'MRS_aINS': anterior insula voxel used for spectroscopy
%
% OUTPUTS
% sxyz_ROI: coordinates of the voxel

%% study
if ~exist('study_nm','var') || ~strcmp(study_nm,'study1')
    study_nm = 'study1';
elseif exist('study_nm','var') && ~strcmp(study_nm,'study1')
    error('case not ready yet');
end

%% working directory
voxel_path = fullfile('E:',study_nm,['CID',sub_nm],'MRS','MRS_voxels');
%% file name
switch ROI_nm
    case 'MRS_dmPFC'
        roi_filenm = 'bwdmpfc.nii';
    case 'MRS_aINS'
        roi_filenm = 'bwai.nii';
end
fullfilename = fullfile(voxel_path, roi_filenm);

%% extract data
if exist(fullfilename,'file')
    % same as in ROI_selection.m to get the coordinates
    % extract ROI mask as a 3-D matrix
    ROI_mask = spm_data_read(fullfilename);
    % extract index of non-zero voxels (= voxels of the mask = coordinates in voxel-space of the mask)
    [sx, sy, sz] = ind2sub(size(ROI_mask), find(ROI_mask > 0));
    % combine the three dims to one variable + add a final column
    % of 1 (necessary for proper matrix multiplication)
    sxyz = [sx(:), sy(:), sz(:), ones(length(sx),1)];
    
    %% convert from voxel space to MNI coordinates
    ROI_vol = spm_vol(fullfilename); % extract info about voxel size and position of the center in 4x4 matrix
    sxyz_ROI = sxyz * ROI_vol.mat'; % convert mask from voxel space to MNI space
else
    sxyz_ROI = [];
end

end % function