%% script to test if ROI extraction works well

%% identify location of brodmann template (one Brodmann area number within 
% Brodmann area => easy to control all works good)
brodmann_path = fullfile('C:','Users','clairis','Desktop',...
    'MRIcron_windows','MRIcron','Resources','templates','brodmann.nii',...
    'brodmann.nii');

%% define coordinates to check
ROI_coord_center = [9 36 34]; % area 32
% ROI_coord_center = [13 61 11]; % area 10
% ROI_coord_center = [13 41 -10]; % frontier between area 10 and 11 => should give value between 10 and 11

%% create sphere
% define radius of the sphere to use
r = 3;
% create sphere with a given radius
% first use ndgrid to create a cube from -radius to +radius around [0,0,0]
[sx, sy, sz] = ndgrid(-r:r, -r:r, -r:r);
% combine the three dims to one variable
sxyz = [sx(:), sy(:), sz(:)];
% compute the distance of each coordinate from [0,0,0]
srad = sqrt(sum(sxyz .* sxyz, 2));
% remove voxels/coordinates where the distance > radius
sxyz(srad > r, :) = [];

%% relocate the sphere at the desired coordinates (in mm == MNI space)
roimm = ROI_coord_center;
sxyz_ROI = sxyz + ones(size(sxyz, 1), 1)*roimm;
% sxyz is just a sphere of the diameter we input at first
% sxyz_ROI is the same sphere but centered on the specific ROI
% coordinates

% for matrix multiplication, we need a column of 1's at the end
sxyz_ROI(:, 4) = 1;

%% extract the data
betaVol     = spm_vol(brodmann_path);
betadata    = spm_read_vols(betaVol);
vxyz        = unique(floor((inv(betaVol.mat) * sxyz_ROI')'), 'rows');
vi          = sub2ind(betaVol.dim, vxyz(:, 1), vxyz(:, 2), vxyz(:, 3));
con_value   = mean(betadata(vi),'omitnan');