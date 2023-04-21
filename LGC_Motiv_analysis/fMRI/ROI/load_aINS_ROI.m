function[aINS_ROI_infos] = load_aINS_ROI()
% [aINS_ROI_infos] = load_aINS_ROI()
% load_aINS_ROI will extract informations about anterior insula mask
%
% OUTPUT
% aINS_ROI_infos: big structure with all the relevant information

%% path
gitFolder = fullfile('C:','Users','clairis','Desktop','GitHub',...
    'LGC_motiv','Matlab_DIY_functions',...
    'ROI','NicoC_masks','aIN');
maskName = 'aIN_BartraNeg3_AAL-IN';
mask_img = [gitFolder, filesep, maskName,'.nii'];
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
aINS_ROI_infos.ROI_xyz.ROI_1 = sxyz_ROI;
aINS_ROI_infos.ROI_sphere_or_mask.ROI_1 = 1;
aINS_ROI_infos.ROI_nm.ROI_1 = 'aINS';
aINS_ROI_infos.ROI_nm.ROI_1_shortName = 'aINS';
aINS_ROI_infos.ROI_nm.fullpath.ROI_1 = mask_img;
aINS_ROI_infos.n_ROIs = 1;

end % function