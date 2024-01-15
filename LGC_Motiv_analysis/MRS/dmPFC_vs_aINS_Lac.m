
%% subjects selection
[study_nm, condition, gender, subject_id, NS] = sub_id;

%% load metabolites
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac = metabolites.dmPFC.Lac;
aIns_Lac = metabolites.aIns.Lac;

%% extract and correlate dmPFC/dACC Lac vs aINS Lac
[r_corr, betas, pval, y_fit,...
    aIns_Lac_sorted, dmPFC_Lac_fit_aInsSorted] = glm_package(aIns_Lac, dmPFC_Lac, 'normal', 'on');

% figure
[pSize, lWidth, col, mSize] = general_fig_prm;
fig;
scat_hdl = scatter(aIns_Lac, dmPFC_Lac);
fit_hdl = plot(aIns_Lac_sorted, dmPFC_Lac_fit_aInsSorted);
scat_hdl_upgrade(scat_hdl);
fit_hdl_upgrade(fit_hdl);
xlabel('aINS Lac (mM)');
ylabel('dmPFC/dACC Lac (mM)');
place_r_and_pval(r_corr, pval(2));
legend_size(pSize);