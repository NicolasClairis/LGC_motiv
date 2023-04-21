function[dmPFC_ROI_infos] = load_dmPFC_ROI()
% [VS_ROI_infos] = load_dmPFC_ROI()
% load_dmPFC_ROI will extract informations about dorsomedial prefrontal mask
%
% OUTPUT
% dmPFC_ROI_infos: big structure with all the relevant information

%% define dmPFC sphere
ROI_coord_center = [-3 18 45];
sphereRadius = 8;
% create sphere with a given radius
% first use ndgrid to create a cube from -radius to +radius around [0,0,0]
[sx, sy, sz] = ndgrid(-sphereRadius:sphereRadius,...
    -sphereRadius:sphereRadius,...
    -sphereRadius:sphereRadius);
% combine the three dims to one variable
sxyz = [sx(:), sy(:), sz(:)];
% compute the distance of each coordinate from [0,0,0]
srad = sqrt(sum(sxyz .* sxyz, 2));
% remove voxels/coordinates where the distance > radius
sxyz(srad > sphereRadius, :) = [];

%% relocate the sphere at the desired coordinates (in mm == MNI space)
roimm = ROI_coord_center;
sxyz_ROI = sxyz + ones(size(sxyz, 1), 1)*roimm;
% sxyz is just a sphere of the diameter we input at first
% sxyz_ROI is the same sphere but centered on the specific ROI
% coordinates

% for matrix multiplication, we need a column of 1's at the end
sxyz_ROI(:, 4) = 1;

%% extract in output
dmPFC_ROI_infos.ROI_xyz.ROI_1 = sxyz_ROI;
dmPFC_ROI_infos.ROI_sphere_or_mask.ROI_1 = 0;
dmPFC_ROI_infos.ROI_nm.ROI_1 = 'dmPFC';
dmPFC_ROI_infos.ROI_nm.ROI_1_shortName = 'dmPFC';
dmPFC_ROI_infos.ROI_nm.ROI_1_centerCoords = num2str(ROI_coord_center); % store coordinates of the center of the sphere
dmPFC_ROI_infos.n_ROIs = 1;

end % function