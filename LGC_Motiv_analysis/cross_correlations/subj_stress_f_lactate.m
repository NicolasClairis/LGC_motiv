function[betas, pval, r_corr] = subj_stress_f_lactate(fig_disp, rmv_outliers_yn)
% [betas, pval, rho] = subj_stress_f_lactate(fig_disp, rmv_outliers_yn)
% subj_stress_f_lactate will compare plasma and brain levels of lactate to
% subjective ratings of stress across the experiment.
%
% INPUTS
% fig_disp: display figure (1) or not (0)? By default will display it
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% betas: regression estimates for correlations
%
% pval: structure with p.value for correlations
%
% rho: correlation coefficient for correlations

%% default inputs
if ~exist('fig_disp','var') || isempty(fig_disp) || ~ismember(fig_disp,[0,1])
    fig_disp = 1;
end
% remove outliers by default
if ~exist('rmv_outliers_yn','var') || isempty(rmv_outliers_yn) || ~ismember(rmv_outliers_yn,[0,1])
    rmv_outliers_yn = 1;
end

%% subject selection
[study_nm, ~, ~, subject_id, NS] = sub_id;

%% initialize variable of interest
[plasma_Lac, dmPFC_Lac, aIns_Lac] = deal(NaN(1,NS));
stress_rating_names = {'preMRS','postMRS','prefMRI','postfMRI',...
    'postfMRI_min_preMRS','postfMRI_min_prefMRI'};
n_stress_rtgs = length(stress_rating_names);
subj_stress = NaN(n_stress_rtgs, NS); % 1 measure pre-MRS, 1 post-MRS, 1 pre-fMRI and 1 post-fMRI + delta pre-MRS/post-fMRI + delta pre-fMRI/post-fMRI

%% load subjective stress
[~,...
    subj_stress(1,:), subj_stress(2,:),...
    subj_stress(3,:), subj_stress(4,:),...
    subj_stress(5,:), subj_stress(6,:)] = extract_subjective_stress_ratings(study_nm, subject_id, NS);

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac(:) = metabolites.dmPFC.Lac;
aIns_Lac(:) = metabolites.aIns.Lac;

%% load plasma lactate
[plasma_Lac_struct] = load_plasma_Lac(subject_id);
plasma_Lac(:) = plasma_Lac_struct.Lac;

%% remove outliers
if rmv_outliers_yn == 1
    [~,~,subj_stress(1,:)] = rmv_outliers_3sd(subj_stress(1,:));
    [~,~,subj_stress(2,:)] = rmv_outliers_3sd(subj_stress(2,:));
    [~,~,subj_stress(3,:)] = rmv_outliers_3sd(subj_stress(3,:));
    [~,~,subj_stress(4,:)] = rmv_outliers_3sd(subj_stress(4,:));
    [~,~,subj_stress(5,:)] = rmv_outliers_3sd(subj_stress(5,:));
    [~,~,subj_stress(6,:)] = rmv_outliers_3sd(subj_stress(6,:));
    [~,~,dmPFC_Lac] = rmv_outliers_3sd(dmPFC_Lac);
    [~,~,aIns_Lac] = rmv_outliers_3sd(aIns_Lac);
    [~,~,plasma_Lac] = rmv_outliers_3sd(plasma_Lac);
end % remove outliers

%% test correlations
for iStress = 1:n_stress_rtgs
    stress_nm = stress_rating_names{iStress};
    
    stress_plasma_nm = [stress_nm,'_stress_f_plasma_Lac'];
    [r_corr.(stress_plasma_nm), betas.(stress_plasma_nm),...
        pval.(stress_plasma_nm),...
        ~, plasma_Lac_sorted,...
        subj_stress_fit.(stress_plasma_nm)] = glm_package(plasma_Lac, subj_stress(iStress,:), 'normal');
    stress_dmPFC_nm = [stress_nm,'_stress_f_dmPFC_Lac'];
    [r_corr.(stress_dmPFC_nm), betas.(stress_dmPFC_nm),...
        pval.(stress_dmPFC_nm),...
        ~, dmPFC_Lac_sorted,...
        subj_stress_fit.(stress_dmPFC_nm)] = glm_package(dmPFC_Lac, subj_stress(iStress,:), 'normal');
    stress_aIns_nm = [stress_nm,'_stress_f_aIns_Lac'];
    [r_corr.(stress_aIns_nm), betas.(stress_aIns_nm),...
        pval.(stress_aIns_nm),...
        ~, aIns_Lac_sorted,...
        subj_stress_fit.(stress_aIns_nm)] = glm_package(aIns_Lac, subj_stress(iStress,:), 'normal');
end % loop over stress periods

%% figure
if fig_disp == 1
    %% aggregate data in one big matrix for the heatmap
    dmPFC_r_vector = [r_corr.preMRS_stress_f_dmPFC_Lac;...
        r_corr.postMRS_stress_f_dmPFC_Lac;...
        r_corr.prefMRI_stress_f_dmPFC_Lac;...
        r_corr.postfMRI_stress_f_dmPFC_Lac;...
        r_corr.postfMRI_min_preMRS_stress_f_dmPFC_Lac;...
        r_corr.postfMRI_min_prefMRI_stress_f_dmPFC_Lac];
    aIns_r_vector = [r_corr.preMRS_stress_f_aIns_Lac;...
        r_corr.postMRS_stress_f_aIns_Lac;...
        r_corr.prefMRI_stress_f_aIns_Lac;...
        r_corr.postfMRI_stress_f_aIns_Lac;...
        r_corr.postfMRI_min_preMRS_stress_f_aIns_Lac;...
        r_corr.postfMRI_min_prefMRI_stress_f_aIns_Lac];
    plasma_r_vector = [r_corr.preMRS_stress_f_plasma_Lac;...
        r_corr.postMRS_stress_f_plasma_Lac;...
        r_corr.prefMRI_stress_f_plasma_Lac;...
        r_corr.postfMRI_stress_f_plasma_Lac;...
        r_corr.postfMRI_min_preMRS_stress_f_plasma_Lac;...
        r_corr.postfMRI_min_prefMRI_stress_f_plasma_Lac];
    corr_mtrx = [dmPFC_r_vector, aIns_r_vector, plasma_r_vector];
    
    % same with p.value
    dmPFC_pval_vector = [pval.preMRS_stress_f_dmPFC_Lac(2);...
        pval.postMRS_stress_f_dmPFC_Lac(2);...
        pval.prefMRI_stress_f_dmPFC_Lac(2);...
        pval.postfMRI_stress_f_dmPFC_Lac(2);...
        pval.postfMRI_min_preMRS_stress_f_dmPFC_Lac(2);...
        pval.postfMRI_min_prefMRI_stress_f_dmPFC_Lac(2)];
    aIns_pval_vector = [pval.preMRS_stress_f_aIns_Lac(2);...
        pval.postMRS_stress_f_aIns_Lac(2);...
        pval.prefMRI_stress_f_aIns_Lac(2);...
        pval.postfMRI_stress_f_aIns_Lac(2);...
        pval.postfMRI_min_preMRS_stress_f_aIns_Lac(2);...
        pval.postfMRI_min_prefMRI_stress_f_aIns_Lac(2)];
    plasma_pval_vector = [pval.preMRS_stress_f_plasma_Lac(2);...
        pval.postMRS_stress_f_plasma_Lac(2);...
        pval.prefMRI_stress_f_plasma_Lac(2);...
        pval.postfMRI_stress_f_plasma_Lac(2);...
        pval.postfMRI_min_preMRS_stress_f_plasma_Lac(2);...
        pval.postfMRI_min_prefMRI_stress_f_plasma_Lac(2)];
    pval_mtrx = [dmPFC_pval_vector, aIns_pval_vector, plasma_pval_vector];
    
    %% general figure parameters
    [pSize, ~, col] = general_fig_prm;
    % define which colormap you want to use (see full list here if you are not
    % happy with the selection:
    % https://ch.mathworks.com/help/matlab/ref/colormap.html)
    % color_range_choices = 'hot';
    % color_range_choices = 'turbo';
    % color_range_choices = 'jet';
    color_range_choices = redblue(45);
    
    % correlation range
    corr_range = [-0.65 0.65];
    
    % x/y-axis ratio
    ax_ratio = [1.5 1 1];
    
    %% heatmap with correlation coefficients
    fig_heatmap = fig;
    subP = subplot(1,2,1);
    imagesc(corr_mtrx, corr_range);
    colormap(subP, color_range_choices);
    cbar = colorbar;
    cbar.Label.String = 'r';
    xticks(1:3);
    xticklabels({'dmPFC/dACC','aIns','plasma'});
    yticks(1:size(corr_mtrx,1));
    yticklabels({'S1','S2','S3','S4','S4-S1','S4-S3'});
    % add stars in the graph if some correlations are significant
    for iLac_measure = 1:size(corr_mtrx,2)
        for iStressRtg = 1:size(corr_mtrx,1)
            if pval_mtrx(iStressRtg, iLac_measure) <= 0.05
                if pval_mtrx(iStressRtg, iLac_measure) > 0.01 && pval_mtrx(iStressRtg, iLac_measure) <= 0.05
                    pval_hdl = text(iLac_measure, iStressRtg, '*');
                elseif pval_mtrx(iStressRtg, iLac_measure) > 0.001 && pval_mtrx(iStressRtg, iLac_measure) <= 0.01
                    pval_hdl = text(iLac_measure, iStressRtg, '**');
                elseif pval_mtrx(iStressRtg, iLac_measure) <= 0.001
                    pval_hdl = text(iLac_measure, iStressRtg, '***');
                end % p.value
                % adjust p.value parameters
                pval_hdl.Color = col.white;
                pval_hdl.FontSize = 70;
                pval_hdl.FontWeight = 'bold';
                pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
                pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
            end % when p.value is significant
        end % loop over Y variables
    end % loop over X variables
end % figure display
end % function