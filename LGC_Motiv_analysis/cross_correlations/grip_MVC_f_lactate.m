function[betas, pval, r_corr] = grip_MVC_f_lactate(fig_disp, rmv_outliers_yn)
% [betas, pval, rho] = grip_MVC_f_lactate(fig_disp, rmv_outliers_yn)
% grip_MVC_f_lactate will compare plasma and brain levels of lactate to
% maximum voluntary contraction (MVC) force in the physical effort task.
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
[~, ~, ~, subject_id, NS] = sub_id;

%% initialize variable of interest
[plasma_Lac, dmPFC_Lac, aIns_Lac,...
    MVC, theoretical_MVC] = deal(NaN(1,NS));

%% load MVC
[~, ~, MVC_subject_id, ~, ~,...
    MVC_all, PCSA_all] = grip_MVC_vs_PCSA(0);
% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    MVC_sub_idx = find(strcmp(sub_nm,MVC_subject_id));
    MVC(iS) = MVC_all(MVC_sub_idx);
    theoretical_MVC(iS) = PCSA_all(MVC_sub_idx);
end % subject loop

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac(:) = metabolites.dmPFC.Lac;
aIns_Lac(:) = metabolites.aIns.Lac;

%% load plasma lactate
[plasma_Lac_struct] = load_plasma_Lac(subject_id);
plasma_Lac(:) = plasma_Lac_struct.Lac./1000;

%% remove outliers
if rmv_outliers_yn == 1
    [~,~,MVC] = rmv_outliers_3sd(MVC);
    [~,~,theoretical_MVC] = rmv_outliers_3sd(theoretical_MVC);
    [~,~,dmPFC_Lac] = rmv_outliers_3sd(dmPFC_Lac);
    [~,~,aIns_Lac] = rmv_outliers_3sd(aIns_Lac);
    [~,~,plasma_Lac] = rmv_outliers_3sd(plasma_Lac);
end % remove outliers

%% test correlations
% tests with MVC
[r_corr.MVC_f_plasma_Lac, betas.MVC_f_plasma_Lac,...
    pval.MVC_f_plasma_Lac,...
    ~, plasma_Lac_sorted,...
    MVC_fit_f_plasma_Lac] = glm_package(plasma_Lac, MVC, 'normal');
[r_corr.MVC_f_dmPFC_Lac, betas.MVC_f_dmPFC_Lac,...
    pval.MVC_f_dmPFC_Lac,...
    ~, dmPFC_Lac_sorted,...
    MVC_fit_f_dmPFC_Lac] = glm_package(dmPFC_Lac, MVC, 'normal');
[r_corr.MVC_f_aIns_Lac, betas.MVC_f_aIns_Lac,...
    pval.MVC_f_aIns_Lac,...
    ~, aIns_Lac_sorted,...
    MVC_fit_f_aIns_Lac] = glm_package(aIns_Lac, MVC, 'normal');

% same but with theoretical MVC
[r_corr.theoretical_MVC_f_plasma_Lac, betas.theoretical_MVC_f_plasma_Lac,...
    pval.theoretical_MVC_f_plasma_Lac,...
    ~, plasma_Lac_sorted,...
    theoretical_MVC_fit_f_plasma_Lac] = glm_package(plasma_Lac, theoretical_MVC, 'normal');
[r_corr.theoretical_MVC_f_dmPFC_Lac, betas.theoretical_MVC_f_dmPFC_Lac,...
    pval.theoretical_MVC_f_dmPFC_Lac,...
    ~, ~,...
    theoretical_MVC_fit_f_dmPFC_Lac] = glm_package(dmPFC_Lac, theoretical_MVC, 'normal');
[r_corr.theoretical_MVC_f_aIns_Lac, betas.theoretical_MVC_f_aIns_Lac,...
    pval.theoretical_MVC_f_aIns_Lac,...
    ~, ~,...
    theoretical_MVC_fit_f_aIns_Lac] = glm_package(aIns_Lac, theoretical_MVC, 'normal');

%% figure
if fig_disp == 1
    % general parameters
    [pSize, lWidth, col, mSize] = general_fig_prm;
    mSize = 100;
    lWidthScatter = 1.5;

    %% MVC = f(plasma Lac)
    fig;
    hold on;
    scatter(plasma_Lac, MVC,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',col.black,'MarkerFaceColor',col.grey);
    plot(plasma_Lac_sorted, MVC_fit_f_plasma_Lac,...
        'LineWidth',lWidth,'LineStyle','-','Color',col.black);
    ylabel('MVC (N)');
    xlabel('Plasma lactate (mM)');
    place_r_and_pval(r_corr.MVC_f_plasma_Lac, pval.MVC_f_plasma_Lac(2));
    legend_size(pSize);
    
    %% MVC = f(dmPFC/dACC Lac)
    fig;
    hold on;
    scatter(dmPFC_Lac, MVC,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',col.black,'MarkerFaceColor',col.grey);
    plot(dmPFC_Lac_sorted, MVC_fit_f_dmPFC_Lac,...
        'LineWidth',lWidth,'LineStyle','-','Color',col.black);
    ylabel('MVC (N)');
    xlabel('dmPFC/dACC lactate (mM)');
    place_r_and_pval(r_corr.MVC_f_dmPFC_Lac, pval.MVC_f_dmPFC_Lac(2));
    legend_size(pSize);
    
    %% MVC = f(aIns Lac)
    fig;
    hold on;
    scatter(aIns_Lac, MVC,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',col.black,'MarkerFaceColor',col.grey);
    plot(aIns_Lac_sorted, MVC_fit_f_aIns_Lac,...
        'LineWidth',lWidth,'LineStyle','-','Color',col.black);
    ylabel('MVC (N)');
    xlabel('aIns lactate (mM)');
    place_r_and_pval(r_corr.MVC_f_aIns_Lac, pval.MVC_f_aIns_Lac(2));
    legend_size(pSize);
    
    %% theoretical MVC = f(plasma Lac)
    fig;
    hold on;
    scatter(plasma_Lac, theoretical_MVC,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',col.black,'MarkerFaceColor',col.grey);
    plot(plasma_Lac_sorted, theoretical_MVC_fit_f_plasma_Lac,...
        'LineWidth',lWidth,'LineStyle','-','Color',col.black);
    ylabel('theoretical MVC (N)');
    xlabel('Plasma lactate (mM)');
    place_r_and_pval(r_corr.theoretical_MVC_f_plasma_Lac, pval.theoretical_MVC_f_plasma_Lac(2));
    legend_size(pSize);
    
    %% theoretical MVC = f(dmPFC/dACC Lac)
    fig;
    hold on;
    scatter(dmPFC_Lac, theoretical_MVC,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',col.black,'MarkerFaceColor',col.grey);
    plot(dmPFC_Lac_sorted, theoretical_MVC_fit_f_dmPFC_Lac,...
        'LineWidth',lWidth,'LineStyle','-','Color',col.black);
    ylabel('theoretical MVC (N)');
    xlabel('dmPFC/dACC lactate (mM)');
    place_r_and_pval(r_corr.theoretical_MVC_f_dmPFC_Lac, pval.theoretical_MVC_f_dmPFC_Lac(2));
    legend_size(pSize);
    
    %% theoretical MVC = f(aIns Lac)
    fig;
    hold on;
    scatter(aIns_Lac, theoretical_MVC,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',col.black,'MarkerFaceColor',col.grey);
    plot(aIns_Lac_sorted, theoretical_MVC_fit_f_aIns_Lac,...
        'LineWidth',lWidth,'LineStyle','-','Color',col.black);
    ylabel('theoretical MVC (N)');
    xlabel('aIns lactate (mM)');
    place_r_and_pval(r_corr.theoretical_MVC_f_aIns_Lac, pval.theoretical_MVC_f_aIns_Lac(2));
    legend_size(pSize);
end % figure
end % function