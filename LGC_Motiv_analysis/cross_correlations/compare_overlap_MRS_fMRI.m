function[overlap] = compare_overlap_MRS_fMRI
% [overlap] = compare_overlap_MRS_fMRI
% compare_overlap_MRS_fMRI will compare how much the 1H-MRS dmPFC/dACC
% voxel and the fMRI dmPFC/dACC cluster overlap with each other using
% various levels of threshold on the 1H-MRS density map (each number
% corresponding to a given proportion of subjects who have to be included
% in a given voxel)
%
% OUTPUTS
% overlap: structure with the result for the fMRI and the MRS voxel
% depending on the threshold used for the density map

%% extract number of subjects
[study_nm, ~, ~, ~, NS] = sub_id;
NS_str = num2str(NS);

%% working directory
MRS_folder = fullfile('E:',study_nm,'MRS_mask');
gitFolder = fullfile('C:','Users','clairis','Desktop','GitHub');
fMRI_folder = fullfile(gitFolder,'LGC_motiv','Matlab_DIY_functions',...
    'ROI','NicoC_masks','GLM-based_masks',['LGC_',study_nm]);

%% load MRS mask
dmPFC_MRS_mask_vol = spm_vol([MRS_folder,filesep,...
    'density_map_dmPFC_MRS_voxel__fMRI_noSatTaskSub_noSatRun_',NS_str,'_subjects.nii']);
dmPFC_MRS_mask_data = spm_read_vols(dmPFC_MRS_mask_vol);
%% load fMRI cluster
GLM_str = 'GLM235';
dmPFC_fMRI_cluster_vol = spm_vol([fMRI_folder,filesep,...
    GLM_str,'_EpEm_dmPFC_Ech_pFWE005voxelCorr_',NS_str,'subs.nii']);
dmPFC_fMRI_cluster_data = spm_read_vols(dmPFC_fMRI_cluster_vol);

%% extract total number of voxels in each mask to compute percentage at the end
% extract number of voxels depending on density map threshold
n_dmPFC_MRS.density0 = sum(sum(sum(dmPFC_MRS_mask_data ~= 0)));
n_dmPFC_MRS.density20 = sum(sum(sum(dmPFC_MRS_mask_data > 0.2)));
n_dmPFC_MRS.density40 = sum(sum(sum(dmPFC_MRS_mask_data > 0.4)));
n_dmPFC_MRS.density50 = sum(sum(sum(dmPFC_MRS_mask_data > 0.5)));
n_dmPFC_MRS.density60 = sum(sum(sum(dmPFC_MRS_mask_data > 0.6)));
n_dmPFC_MRS.density80 = sum(sum(sum(dmPFC_MRS_mask_data > 0.8)));
n_dmPFC_MRS.density90 = sum(sum(sum(dmPFC_MRS_mask_data > 0.9)));
% extract number of voxels in fMRI cluster
n_dmPFC_fMRI_cluster = sum(sum(sum(dmPFC_fMRI_cluster_data ~= 0)));

%% compute overlap
[n_overlap.density0, n_overlap.density20,...
    n_overlap.density40, n_overlap.density50, n_overlap.density60,...
    n_overlap.density80, n_overlap.density90] = deal(0);
for iX = 1:size(dmPFC_MRS_mask_data,1)
    for iY = 1:size(dmPFC_MRS_mask_data, 2)
        for iZ = 1:size(dmPFC_MRS_mask_data, 3)
            
            % filter if voxel present in fMRI cluster
            if (dmPFC_fMRI_cluster_data(iX, iY, iZ)) ~= 0
                
                % filter according to proportion of subjects included in
                % MRS voxel
                if dmPFC_MRS_mask_data(iX, iY, iZ) > 0
                    n_overlap.density0 = n_overlap.density0 + 1;
                    if dmPFC_MRS_mask_data(iX, iY, iZ) > 0.2
                        n_overlap.density20 = n_overlap.density20 + 1;
                        if dmPFC_MRS_mask_data(iX, iY, iZ) > 0.4
                            n_overlap.density40 = n_overlap.density40 + 1;
                            if dmPFC_MRS_mask_data(iX, iY, iZ) > 0.5
                                n_overlap.density50 = n_overlap.density50 + 1;
                                if dmPFC_MRS_mask_data(iX, iY, iZ) > 0.6
                                    n_overlap.density60 = n_overlap.density60 + 1;
                                    if dmPFC_MRS_mask_data(iX, iY, iZ) > 0.8
                                        n_overlap.density80 = n_overlap.density80 + 1;
                                        if dmPFC_MRS_mask_data(iX, iY, iZ) > 0.9
                                            n_overlap.density90 = n_overlap.density90 + 1;
                                        end % 90% MRS voxel
                                    end % 80% MRS voxel
                                end % 60% MRS voxel
                            end % 50% MRS voxel
                        end % 40% MRS voxel
                    end % 20% MRS voxel
                end % filter according to if some subjects in MRS voxel
                
            end % fMRI filter
        end % loop over z dimension
    end % loop over y dimension
end % loop over x dimension

%% assess percentage
% percentage overlap for the MRS voxel
overlap.density0.MRS_voxel = n_overlap.density0./n_dmPFC_MRS.density0;
overlap.density20.MRS_voxel = n_overlap.density20./n_dmPFC_MRS.density20;
overlap.density40.MRS_voxel = n_overlap.density40./n_dmPFC_MRS.density40;
overlap.density50.MRS_voxel = n_overlap.density50./n_dmPFC_MRS.density50;
overlap.density60.MRS_voxel = n_overlap.density60./n_dmPFC_MRS.density60;
overlap.density80.MRS_voxel = n_overlap.density80./n_dmPFC_MRS.density80;
overlap.density90.MRS_voxel = n_overlap.density90./n_dmPFC_MRS.density90;
% percentage overlap for the fMRI cluster
overlap.density0.fMRI_cluster = n_overlap.density0./n_dmPFC_fMRI_cluster;
overlap.density20.fMRI_cluster = n_overlap.density20./n_dmPFC_fMRI_cluster;
overlap.density40.fMRI_cluster = n_overlap.density40./n_dmPFC_fMRI_cluster;
overlap.density50.fMRI_cluster = n_overlap.density50./n_dmPFC_fMRI_cluster;
overlap.density60.fMRI_cluster = n_overlap.density60./n_dmPFC_fMRI_cluster;
overlap.density80.fMRI_cluster = n_overlap.density80./n_dmPFC_fMRI_cluster;
overlap.density90.fMRI_cluster = n_overlap.density90./n_dmPFC_fMRI_cluster;

end % function