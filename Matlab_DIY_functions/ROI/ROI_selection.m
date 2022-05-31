function [ ROI_xyz, ROI_sphere_or_mask, ROI_nm, nb_ROIs, ROI_vol, ROI_mask ] = ROI_selection(gitFolder)
%[ ROI_xyz, ROI_sphere_or_mask, ROI_nm, nb_ROIs, ROI_vol, ROI_mask ] = ROI_selection(gitFolder)
% proposes different ROI and asks you which one(s) you want to use.
%
% Note that this script requires the SPM toolbox (version 12 works,
% previous versions should be tested)
%
% INPUTS
% gitFolder: path where git ROI scripts are stored
%
% OUTPUTS
% ROI_xyz: structure containing one subfield per ROI named as
% 'ROI_1/2/3/...'
% Each subfield entails a 3-D matrix containing the (xyz) coordinates of
% the ROI in MNI and one column of ones to be used for the ROI conversion
% in the voxel space by ROI_extraction.m
% Format is (nb_ROIs rows * ([x, y, z, 1] columns)
%
% ROI_sphere_or_mask: structure indicating, for each ROI, if it's a sphere
% (=0) or a mask (=1). Organized similarly to ROI_xyz with one subfield per
% ROI named as 'ROI_1/2/3/...'
%
% ROI_nm: structure indicating, for each ROI, the name of the mask if it's
% a mask or the coordinates of the center of the sphere if it's a sphere.
% Organized similarly to ROI_xyz with one subfield per
% ROI named as 'ROI_1/2/3/...'
% ROI_nm also contains a subfield 'fullpath' which contains the full path
% for each ROI
%
% nb_ROIs: number indicating the number of ROI selected by the subject
% inside this script.
%
% ROI_vol: structure containing, in the mat subfield, a 4x4 matrix with
% informations about the ROI voxel-space (3 first columns correspond to the
% x,y,z info; 4th column corresponds to the rescaling needed for rescaling)
%
% ROI_mask: matrix the size of the whole brain with 1 for the voxels of the
% ROI
%
% See also ROI_extraction
%
% Written by N.Clairis 14/06/18

%% number of ROIS to test
nb_ROIs = spm_input('How many ROI?',1,'e',[]);

for iROI = 1:nb_ROIs
    
    %% sphere or mask?
    % use sphere (manually defining sphere size centered on coordinates
    % selected in the list provided or manually defined) or a mask?
    sphere_mask = spm_input('Use a sphere or a mask?',1,'m','sphere| mask',[0 1], 0);
    
    %% define ROI
    if sphere_mask == 0
        ROI_center = spm_input('What ROI do you want to analyze?',1,'m',...
            ['manual input |'...
            'vmPFC value ((-10,48,-12) Clairis 2022) |'...
            'vmPFC ((0,40,-12) Bartra 2013) |'...
            'vmPFC value ((-10,44,-8) Lebreton 2009) |'...
            'vmPFC value ((-2,52,-2) Lebreton 2015 |'...
            'vmPFC/pgACC ((2,40,10) Lopez 2017) |'...
            'left NAcc ((-16,4,-4) Bartra 2013) |'...
            'right Nacc ((18,6,-6), Bartra 2013) |'...
            'left hippocampus ((-18,-36,-10) Lebreton 2009) |'...
            'left hippocampus ((-28,-18,-16) Neurosynth) |'...
            'right hippocampus ((28,-18,-16) Neurosynth) |'...
            'mPFC Confidence ((-8,52,18) Clairis 2022) |'...
            'dmPFC DT/Effort ((10,12,48) Clairis 2022) |'...
            'left dmFC ((-2,16,48) Engstrom 2015) |'...
            'left dmFC ((-8,20,46) Lopez 2017) |'...
            'right dmFC ((4,24,40) Engstrom 2015) |'...
            'right dmFC ((8,18,46) Lopez 2017) |'...
            'left anterior insula ((-32,26,0) Bartra 2013) |'...
            'left anterior insula ((-30,26,2) Lopez 2017) |'...
            'right anterior insula ((32,20,-6) Bartra 2013) |'...
            'right anterior insula ((30,28,0) Lopez 2017) |'...
            'left SMA ((-9, -7, 58) Klein-Flugge 2016) |'...
            'left dACC ((-3, 11, 34) Klein-Flugge 2016) |'...
            'right dACC ((10,26,34) Lopez 2017)'], ...
            1:24, 0);
        
        switch ROI_center
            case 1 % manual input
                ROI_coord_center = spm_input('Enter x y z coordinates of your ROI',1,'c');
            case 2 % vmPFC value ((-10,48,-12) Clairis 2022)
                ROI_coord_center = [-10, 48, -12];
            case 3 % vmPFC ((0,40,-12) Bartra 2013)
                ROI_coord_center = [0, 40, -12];
            case 4 % vmPFC ((-10,44,-8) Lebreton 2009)
                ROI_coord_center = [-10, 44, -8];
            case 5 % vmPFC ((-2,52,-2) Lebreton 2015)
                ROI_coord_center = [-2, 52, -2];
            case 6 % vmPFC/pgACC ((2,40,10) Lopez 2017)
                ROI_coord_center = [2, 40, 10];
            case 7 % left NAcc ((-16,4,-4) Bartra 2013)
                ROI_coord_center = [-16, 4,-4];
            case 8 % right Nacc ((18,6,-6), Bartra 2013)
                ROI_coord_center = [18, 6, -6];
            case 9 % left hippocampus ((-18,-36,-10) Lebreton 2009)
                ROI_coord_center = [-18, -36, -10];
            case 10 % left hippocampus ((-28, -18, -16) Neurosynth)
                ROI_coord_center = [-28, -18, -16];
            case 11 % right hippocampus ((28, -18, -16) Neurosynth)
                ROI_coord_center = [28, -18, -16];
            case 12 % mPFC confidence ((-8,52,18) Clairis 2022)
                ROI_coord_center = [-8, 52, 18];
            case 13 % dmPFC DT/Effort ((10,12,48) Clairis 2022)
                ROI_coord_center = [10, 12, 48];
            case 14 % left dmFC ((-2,16,48) Engstrom 2015)
                ROI_coord_center = [-2, 16, 48];
            case 15 % left dmFC ((-8,20,46) Lopez 2017)
                ROI_coord_center = [-8,20,46];
            case 16 % right dmFC ((4,24,40) Engstrom 2015)
                ROI_coord_center = [4, 24, 40];
            case 17 % right dmFC ((8,18,46) Lopez 2017)
                ROI_coord_center = [8,18,46];
            case 18 % left anterior insula ((-32,26,0) Bartra 2013)
                ROI_coord_center = [-32, 26, 0];
            case 19 % left anterior insula ((-30,26,2) Lopez 2017)
                ROI_coord_center = [-30,26,2];
            case 20 % right anterior insula ((32,20,-6) Bartra 2013)
                ROI_coord_center = [32, 20, -6];
            case 21 % right anterior insula ((30,28,0) Lopez 2017)
                ROI_coord_center = [30,28,0];
            case 22 % left SMA ((-9, -7, 58) Klein-Flugge 2016)
                ROI_coord_center = [-9, -7, 58];
            case 23 % left dACC ((-3, 11, 34) Klein-Flugge 2016)
                ROI_coord_center = [-3, 11, 34];
            case 24 % right dACC ((10,26,34) Lopez 2017)
                ROI_coord_center = [10,26,34];
        end
        
        %% name for the sphere
        sphere_nm = spm_input('Short name for the sphere?',1,'s');
        
        %% sphere size
        r = 0;
        while r <= 0
            r = spm_input('Sphere radius (in mm)?',1,'e',8);
        end
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
        
        %% extract in output
        ROI_xyz.(['ROI_',num2str(iROI)])                = sxyz_ROI;
        ROI_sphere_or_mask.(['ROI_',num2str(iROI)])     = sphere_mask;
        ROI_nm.(['ROI_',num2str(iROI)])                 = [sphere_nm,'_',num2str(r),'mm_sphere']; % should not start with a number or next scripts will bug
        ROI_nm.(['ROI_',num2str(iROI),'_centerCoords']) = num2str(ROI_coord_center); % store coordinates of the center of the sphere
        
        [ROI_vol, ROI_mask] = deal([]);
        
        %% using a mask
    elseif sphere_mask == 1
        
        maskFolder = [gitFolder, filesep, 'NicoC_masks', filesep];
        
        % choice of the area(s) on which you want to focus
        possibleMaskAreas = {'vmPFC','PCC',...
            'Striatum','aIN',...
            'dACC','mPFC',...
            'dmPFC','dlPFC',...
            'GLM-based_masks'};
        nb_possible_areas = length(possibleMaskAreas);
        
        ROI_area = spm_input('What area(s) do you want to analyze?',1,'m',...
            ['vmPFC |',...
            'PCC |',...
            'striatum |',...
            'anterior insula |',...
            'dACC |',...
            'mPFC |',...
            'dmPFC |',...
            'dlPFC |',...
            'GLM-based masks'],...
            1:nb_possible_areas, 0);
        
        % define mask to be used
        mask_img = spm_select(1,'.nii','Which mask do you want to use?',{},...
            [maskFolder, possibleMaskAreas{ROI_area}]);
        
        % store name for each area
        switch ROI_area
            case 1
                maskName = 'vmPFC';
            case 2
                maskName = 'PCC';
            case 3
                maskName = 'striatum';
            case 4
                maskName = 'aIN';
            case 5
                maskName = 'dACC';
            case 6
                maskName = 'mPFC';
            case 7
                maskName = 'dmPFC';
            case 8
                maskName = 'dlPFC';
            case 9
                % extract index of the last filesep
                idx_lastFilesep = find(mask_img == '\',1,'last');
                maskName = mask_img(idx_lastFilesep+1:end-4);
        end
        
        
        %% extract ROI coordinates in MNI space
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
        ROI_nm.(['ROI_',num2str(iROI)])             = maskName;
        ROI_nm.fullpath.(['ROI_',num2str(iROI)])    = mask_img;
        
    end % sphere or mask
    
end % ROI loop

end % function

