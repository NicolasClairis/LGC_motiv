function [ con_value ] = ROI_extraction( s_con, iCon, path, sxyz_ROI, beta_or_t_value )
%[ con_value ] = ROI_extraction( s_con, iCon, path, sxyz_ROI, beta_or_t_value )
% extract the mean beta in the voxels of sxyz_ROI for the contrast
% identified with the number iCon in the folder [path].
%
% Note that this function requires the SPM toolbox to be added in your path
% to be able to work. (version 12 works for sure, the previous ones have to 
% be tested)
%
% INPUTS
% s_con:
% 's': use contrasts smoothed at the first level (scon files)
% '': use unsmoothed contrasts (con files)
%
% iCon: number of the contrast you want to extract (format: number). Should
% be the number which is in the file name.
%
% path: folder where the contrast in which you want to extract your ROI is.
%
% sxyz_ROI: 2-D matrix containing the coordinates of the ROI to extract in
% MNI coordinates (=mm). (Format: nVoxels*(x, y ,z, ones) the fourth and
% last column should only be made of ones to be able to do the conversion
% afterwards in the voxel-space)
%
% beta_or_t_value: extract beta ('beta_value') or t.value ('t_value') for the ROI
%
% OUTPUTS
% con_value:
% 1 number corresponding to the average of the betas inside the ROI of
% interest for the contrast selected in the input.
%
% See also ROI_selection
%
% Written by N.Clairis - june 2018 (adapted from J.Geerts 2016)


%% path
root = pwd;
cd(path);

switch beta_or_t_value
    case 'beta_value'
        %% extract smoothed or unsmoothed contrasts?
        if strcmp(s_con,'s')
            s_con_nm = 'scon_';
        elseif strcmp(s_con,'')
            s_con_nm = 'con_';
        end
    case 't_value'
        s_con_nm = 'spmT_';
end

%% add relevant zeros number to contrast name
if iCon < 10
    conZeros = '000';
elseif iCon >= 10 && iCon < 100
    conZeros = '00';
elseif iCon >= 100 && iCon < 1000
    conZeros = '0';
else
    conZeros = '';
end

%% extract the full whole-brain contrast matrix
betaNum     = strcat([s_con_nm, conZeros, num2str(iCon), '.nii']);% extract name of the file
betaVol     = spm_vol(betaNum); % extracts header to read con file
betadata    = spm_read_vols(betaVol); % reads the con file

%% extract indices inside the con file based on the
% MNI coordinates entered previously for the
% sphere/mask
vxyz        = unique(floor((inv(betaVol.mat) * sxyz_ROI')'), 'rows');   % converts the ROI to extract from MNI (mm) to voxel-space of your own specific study
vi          = sub2ind(betaVol.dim, vxyz(:, 1), vxyz(:, 2), vxyz(:, 3)); % extracts coordinates of the ROI in the voxel-space of your study
con_value   = nanmean(betadata(vi));                                    % extracts mean beta for the selected ROI sphere/mask


%% script should be adapted for the case where ROI comes from another study with different voxel size and gap between slices.
% Indeed, if so, there will be gaps between voxels that spm_dilate.m should
% be able to fill.
% Then you can use spm_erode.m to limit the inflated number of voxels
% Otherwise you can use it that way but some voxels might be missed


%% go back to root
cd(root);

end