function [  ] = extract_cluster_from_map(  )
% [  ] = extract_cluster_from_map(  )
% to select one or multiple clusters from a .nii map in order to create a
% new ROI containing only the selected clusters/voxels.

%% working directory
root = fullfile('B:','resultats','ROI_definition');
[mask_filenm, select_ok] = spm_select(1,'.nii','Select images you want to pool',{},root);
if select_ok == 0
    warning('problem in file selection');
    return;
end
mask_img = spm_data_read(mask_filenm);
mask_vol = spm_vol(mask_filenm);
full_brain = ones( size(mask_img));

%% extract all the clusters
[x_full,y_full,z_full] = ind2sub( size( mask_img ), find(mask_img ~= 0 | mask_img == 0));
[x, y, z] = ind2sub( size( mask_img ), find(mask_img ~= 0));

clusters = spm_clusters([x,y,z]');
n_clusters = length( unique(clusters) );
%% select the mask you want to extract
fig();

cluster_to_include = [];

for iC = 1:n_clusters
    
    scatter3(x_full,y_full,z_full);
    hold on;
    voxels_idx_tmp = clusters == iC;
    scatter3(x(voxels_idx_tmp), y(voxels_idx_tmp), z(voxels_idx_tmp),'r' );
    
    answer = questdlg('Inclure ce cluster?');
    if strcmp(answer,'Yes')
        clusters = [cluster_to_include, iC];
    end
    hold off;
end % cluster loop

%% extract the index for all the clusters included
x_voxels_idx = x(clusters(cluster_to_include));
y_voxels_idx = y(clusters(cluster_to_include));
z_voxels_idx = z(clusters(cluster_to_include));

%% save an ROI with the voxels of the ROI you want to save
new_mask = zeros( size(mask_img) );
new_mask( [x_voxels_idx, y_voxels_idx, z_voxels_idx] ) = 1;

ROI_nm = inputdlg('Name for the new ROI?');

% spm_data_write(, new_mask)

end % function