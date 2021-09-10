function [ new_mask_data ] = isolate_nii_cluster(  )
%[ new_mask_data ] = isolate_nii_cluster(  )
% isolate_nii_cluster goal is to help you to isolate one cluster from a .nii map 
% and save a new .nii file based on the selected clusters.
%
% OUTPUTS
% new_mask_data: structure with info for the output file after selecting
% the clusters of interest

%% working directory
root = pwd; % by default start where you are currently

%% select baseline .nii file
original_mask_img_nm = spm_select(1,'.nii','Which mask do you want to use?',{},root);

%% open file
original_mask_data = spm_data_read(original_mask_img_nm);
hdr = spm_data_hdr_read(original_mask_img_nm);

%% get info about all non-zero voxels
% extract index of non-zero voxels (= voxels of the mask = coordinates in voxel-space of the mask)
[sx, sy, sz] = ind2sub(size(original_mask_data), find(original_mask_data > 0));

%% extract clusters index
clusters_idx = spm_clusters([sx,sy,sz]');
n_clusters = nanmax(clusters_idx); 

%% select the cluster of interest
fig();
hold on;
keep_cluster = false(1,n_clusters);
for iC = 1:n_clusters
    hold off;
    
    % plot all clusters
    scatter3(sx, sy, sz, 'filled');
    hold on;
    
    % extract coordinates of the current cluster
    sx_tmp = sx(clusters_idx == iC);
    sy_tmp = sy(clusters_idx == iC);
    sz_tmp = sz(clusters_idx == iC);
    
    scatter3(sx_tmp, sy_tmp, sz_tmp, 'filled',...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor',[0 0.25 0.25]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    keep_yn = questdlg('Keep this cluster?');
    if strcmp(keep_yn,'Yes')
        keep_cluster(iC) = true;
    end
end % cluster loop

%% create and save file with only the cluster(s) of interest
new_mask_data_coord = zeros(size(original_mask_data));
% clusters to keep
n_clusters_to_keep = sum(keep_cluster);
clusters_idx_to_keep = find(keep_cluster);
for iC_to_keep = 1:n_clusters_to_keep
    sx_to_keep = sx(clusters_idx == clusters_idx_to_keep(iC_to_keep));
    sy_to_keep = sy(clusters_idx == clusters_idx_to_keep(iC_to_keep));
    sz_to_keep = sz(clusters_idx == clusters_idx_to_keep(iC_to_keep));
    n_voxels = length(sx_to_keep);
    for iVoxel = 1:n_voxels
        new_mask_data_coord(sx_to_keep(iVoxel),...
            sy_to_keep(iVoxel),...
            sz_to_keep(iVoxel)) = 1;
    end
end

new_filenm = inputdlg('Name of the new file?');
new_hdr = hdr;
new_hdr.fname = [root,filesep,new_filenm{1},'.nii'];
new_mask_data = spm_data_write(new_hdr,new_mask_data_coord);

end % function