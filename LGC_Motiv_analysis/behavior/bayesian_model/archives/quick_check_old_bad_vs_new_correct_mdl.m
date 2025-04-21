
%% run Arthur's (corrected) model to extract behavioral parameters
[CID_nb, mod] = ModelSelectionAndFeatureExtraction_NC_modified;

% load previous version
old_version_folder = 'C:\Users\clairis\Desktop\GitHub\LGC_motiv\LGC_Motiv_results\study1\bayesian_modeling\archives_wrong_values_due_to_pb_in_prevEfficacy\';
old_prm = getfield(getfield(load([old_version_folder,'behavioral_prm.mat']),'bayesian_mdl'),'mdl_3');
old_subject_id = old_prm.subject_id;
NS = length(old_subject_id);

%% extract corresponding subjects
full_subject_id = convert_sub_id_from_num_to_cell(CID_nb);
% initialize parameters of interest with correct N
[prm.kR, prm.kP, prm.kEp, prm.kEm,...
    prm.kFp, prm.kLm, prm.kBiasM] = deal(NaN(1,NS));
prm_names = fieldnames(prm);
nPrm = length(prm_names);

for iS = 1:NS
   sub_nm = old_subject_id{iS};
   sub_idx = strcmp(sub_nm, full_subject_id);
   for iP = 1:nPrm
       prm_nm = prm_names{iP};
       switch prm_nm
           case 'kBiasM'
               prm.(prm_nm)(iS) = mod.biais(sub_idx);
           otherwise
               prm.(prm_nm)(iS) = mod.(prm_nm)(sub_idx);
       end
   end % parameter loop
end % subject loop

%% display scatter plots
[pSize, lWidth, col] = general_fig_prm;
fig;
for iP = 1:nPrm
    prm_nm = prm_names{iP};
    subplot(2,4,iP);
    scat_hdl = scatter(old_prm.(prm_nm), prm.(prm_nm));
    scat_hdl.MarkerEdgeColor = 'k';
    scat_hdl.MarkerFaceColor = col.grey;
    xlabel(['Old ',prm_nm]);
    ylabel(['Fixed ',prm_nm]);
    
    %% add correlation
    [coef.(prm_nm), pval.(prm_nm)] = corr(old_prm.(prm_nm)', prm.(prm_nm)');
    place_r_and_pval(coef.(prm_nm), pval.(prm_nm));
end % parameter loop