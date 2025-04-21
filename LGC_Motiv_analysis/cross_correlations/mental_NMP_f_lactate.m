function[betas, pval, r_corr] = mental_NMP_f_lactate(fig_disp, rmv_outliers_yn)
% [betas, pval, rho] = mental_NMP_f_lactate(fig_disp, rmv_outliers_yn)
% mental_NMP_f_lactate will compare plasma and brain levels of lactate to
% number of maximal performance (NMP) in the mental effort task.
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
[plasma_Lac, dmPFC_Lac, aIns_Lac, NMP] = deal(NaN(1,NS));

%% load NMP
[NMP(:)] = extract_Nback_NMP('study1', subject_id, NS);

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac(:) = metabolites.dmPFC.Lac;
aIns_Lac(:) = metabolites.aIns.Lac;

%% load plasma lactate
[plasma_Lac_struct] = load_plasma_Lac(subject_id);
plasma_Lac(:) = plasma_Lac_struct.Lac./1000;

%% remove outliers
if rmv_outliers_yn == 1
    [~,~,NMP] = rmv_outliers_3sd(NMP);
    [~,~,dmPFC_Lac] = rmv_outliers_3sd(dmPFC_Lac);
    [~,~,aIns_Lac] = rmv_outliers_3sd(aIns_Lac);
    [~,~,plasma_Lac] = rmv_outliers_3sd(plasma_Lac);
end % remove outliers

%% test correlations
% tests with NMP
[r_corr.NMP_f_plasma_Lac, betas.NMP_f_plasma_Lac,...
    pval.NMP_f_plasma_Lac,...
    ~, plasma_Lac_sorted,...
    NMP_fit_f_plasma_Lac] = glm_package(plasma_Lac, NMP, 'normal');
[r_corr.NMP_f_dmPFC_Lac, betas.NMP_f_dmPFC_Lac,...
    pval.NMP_f_dmPFC_Lac,...
    ~, dmPFC_Lac_sorted,...
    NMP_fit_f_dmPFC_Lac] = glm_package(dmPFC_Lac, NMP, 'normal');
[r_corr.NMP_f_aIns_Lac, betas.NMP_f_aIns_Lac,...
    pval.NMP_f_aIns_Lac,...
    ~, aIns_Lac_sorted,...
    NMP_fit_f_aIns_Lac] = glm_package(aIns_Lac, NMP, 'normal');

%% figure
if fig_disp == 1
    % general parameters
    [pSize, lWidth, col, mSize] = general_fig_prm;
    mSize = 100;
    lWidthScatter = 1.5;

    %% NMP = f(plasma Lac)
    fig;
    hold on;
    scatter(plasma_Lac, NMP,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',col.black,'MarkerFaceColor',col.grey);
    plot(plasma_Lac_sorted, NMP_fit_f_plasma_Lac,...
        'LineWidth',lWidth,'LineStyle','-','Color',col.black);
    ylabel('NMP');
    xlabel('Plasma lactate (mM)');
    place_r_and_pval(r_corr.NMP_f_plasma_Lac, pval.NMP_f_plasma_Lac(2));
    legend_size(pSize);
    
    %% NMP = f(dmPFC/dACC Lac)
    fig;
    hold on;
    scatter(dmPFC_Lac, NMP,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',col.black,'MarkerFaceColor',col.grey);
    plot(dmPFC_Lac_sorted, NMP_fit_f_dmPFC_Lac,...
        'LineWidth',lWidth,'LineStyle','-','Color',col.black);
    ylabel('NMP');
    xlabel('dmPFC/dACC lactate (mM)');
    place_r_and_pval(r_corr.NMP_f_dmPFC_Lac, pval.NMP_f_dmPFC_Lac(2));
    legend_size(pSize);
    
    %% NMP = f(aIns Lac)
    fig;
    hold on;
    scatter(aIns_Lac, NMP,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',col.black,'MarkerFaceColor',col.grey);
    plot(aIns_Lac_sorted, NMP_fit_f_aIns_Lac,...
        'LineWidth',lWidth,'LineStyle','-','Color',col.black);
    ylabel('NMP');
    xlabel('aIns lactate (mM)');
    place_r_and_pval(r_corr.NMP_f_aIns_Lac, pval.NMP_f_aIns_Lac(2));
    legend_size(pSize);
    
end % figure
end % function