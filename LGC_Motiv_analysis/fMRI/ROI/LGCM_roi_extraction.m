%% extract ROI fMRI pilots


%% extract working directories
[computerRoot, spmFolderPath] = LGCM_root_paths();

study_nm = 'fMRI_pilots';
n_GLM = 4;
GLM_str = num2str(n_GLM);
dataRoot = [computerRoot,filesep,study_nm,filesep];
[subject_id, NS] = LGCM_subject_selection(study_nm);

gitFolder = fullfile('C:','Users','Loco','Documents','GitHub','LGC_motiv','Matlab_DIY_functions','ROI');

pSize = 40;

%% select the ROI to use
[ ROI_xyz, ROI_sphere_or_mask, ROI_nm, nb_ROIs, ROI_vol, ROI_mask ] = ROI_selection(gitFolder);

%% select contrasts of interest
beta_or_t_value = 'beta_value';

%% extract the ROI for each contrast of the selected subjects
for iROI = 1:nb_ROIs
    sxyz_ROI = ROI_xyz.(['ROI_',num2str(iROI)]);
    for iS = 1:NS
        sub_nm = subject_id{iS};
        sub_fMRI_path = [fullfile(dataRoot, ['CID',sub_nm], 'fMRI_analysis','functional',['GLM',GLM_str]),filesep];
        
        switch study_nm
            case 'fMRI_pilots'
                switch sub_nm
                    case 'pilot_s1' % 1 run Ep + 1 run Em
                        % define index
                        switch beta_or_t_value
                            case 'b_value'
                                R_Ep = 3;
                                E_Ep = 4;
                                R_Em = 14;
                                E_Em = 15;
                            case 't_value'
                                R_Ep = 5;
                                E_Ep = 7;
                                R_Em = 15;
                                E_Em = 17;
                        end
                        % extract the data
                        [ con_value_R_Ep ] = ROI_extraction( '', R_Ep, sub_fMRI_path, sxyz_ROI, beta_or_t_value );
                        [ con_value_E_Ep ] = ROI_extraction( '', E_Ep, sub_fMRI_path, sxyz_ROI, beta_or_t_value );
                        [ con_value_R_Em ] = ROI_extraction( '', R_Em, sub_fMRI_path, sxyz_ROI, beta_or_t_value );
                        [ con_value_E_Em ] = ROI_extraction( '', E_Em, sub_fMRI_path, sxyz_ROI, beta_or_t_value );
                        fig;
                        bar([con_value_R_Ep, con_value_E_Ep, con_value_R_Em, con_value_E_Em]);
                        legend_size(pSize);
                        xticks(1:4);
                        xticklabels({'R','E','R','E'});
                        switch beta_or_t_value
                            case 'b_value'
                                ylabel('regression estimate');
                            case 't_value'
                                ylabel('t.value');
                        end
                    case 'pilot_s2' % only 1 run Em
                        % define index
                        switch beta_or_t_value
                            case 'b_value'
                                R_Em = 3;
                                E_Em = 5;
                            case 't_value'
                                R_Em = 5;
                                E_Em = 7;
                        end
                        % extract the data
                        [ con_value_R_Em ] = ROI_extraction( '', R_Em, sub_fMRI_path, sxyz_ROI, beta_or_t_value );
                        [ con_value_E_Em ] = ROI_extraction( '', E_Em, sub_fMRI_path, sxyz_ROI, beta_or_t_value );
                        
                        fig;
                        bar([con_value_R_Em, con_value_E_Em]);
                        legend_size(pSize);
                        xticks(1:2);
                        xticklabels({'R','E'});
                        switch beta_or_t_value
                            case 'b_value'
                                ylabel('regression estimate');
                            case 't_value'
                                ylabel('t.value');
                        end
                end
        end
        
    end % subject loop
end % roi loop