function[VS_ROI_infos] = load_VS_ROI()


%% path
gitFolder = 
mask_img
%% get ROI coordinates
% extract ROI mask as a 3-D matrix
ROI_mask = spm_data_read(mask_img);

% extract index of non-zero voxels (= voxels of the mask = coordinates in voxel-space of the mask)
[sx, sy, sz] = ind2sub(size(ROI_mask), find(ROI_mask > 0));
% combine the three dims to one variable + add a final column
% of 1 (necessary for proper matrix multiplication)
sxyz = [sx(:), sy(:), sz(:), ones(length(sx),1)];

%% convert from voxel space to MNI coordinates
ROI_vol = spm_vol(mask_img); % extract info about voxel size and position of the center in 4x4 matrix
sxyz_ROI = sxyz * ROI_vol.mat'; % convert mask from voxel space to MNI space

%% extract in output
ROI_xyz.(['ROI_',num2str(iROI)])            = sxyz_ROI;
ROI_sphere_or_mask.(['ROI_',num2str(iROI)]) = sphere_mask;
ROI_nm.(['ROI_',num2str(iROI),'_shortName'])= maskName;
ROI_nm.(['ROI_',num2str(iROI)])             = maskName;
ROI_nm.fullpath.(['ROI_',num2str(iROI)])    = mask_img;

end % function