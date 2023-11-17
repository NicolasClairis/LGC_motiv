% correlation between blood Lactate and brain metabolites

%% subject identification
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load plasma Lac
[Lac_struct] = load_plasma_Lac();
Lac = NaN(1,NS);
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = strcmp(Lac_struct.CID, ['CID',sub_nm]);
    Lac(iS) = Lac_struct.Lac(sub_idx);
end % subject loop

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac = metabolites.dmPFC.Lac';
aIns_Lac = metabolites.aIns.Lac';

%% test significance
[r_corr.dmPFC, betas.dmPFC, pval.dmPFC,...
    ~, Lac_sorted.dmPFC, Lac_fit_xSorted.dmPFC] = glm_package(Lac', dmPFC_Lac, 'normal', 'on');
[r_corr.aIns, betas.aIns, pval.aIns,...
    ~, Lac_sorted.aIns, Lac_fit_xSorted.aIns] = glm_package(Lac', aIns_Lac, 'normal', 'on');

%% figure
[pSize, lWidth, col, mSize] = general_fig_prm;

fig;
% dmPFC
scat_dmPFC_hdl = scatter(Lac, dmPFC_Lac);
scat_dmPFC_hdl.MarkerEdgeColor = col.green;
scat_dmPFC_hdl.LineWidth = lWidth;
plot_dmPFC_hdl = plot(Lac_sorted.dmPFC, Lac_fit_xSorted.dmPFC);
plot_dmPFC_hdl.Color = col.green;
plot_dmPFC_hdl.LineWidth = lWidth;
plot_dmPFC_hdl.LineStyle = '-';
% anterior insula
scat_AI_hdl = scatter(Lac, aIns_Lac);
scat_AI_hdl.MarkerEdgeColor = col.orange;
scat_AI_hdl.LineWidth = lWidth;
plot_aIns_hdl = plot(Lac_sorted.aIns, Lac_fit_xSorted.aIns);
plot_aIns_hdl.Color = col.orange;
plot_aIns_hdl.LineWidth = lWidth;
plot_aIns_hdl.LineStyle = '-';
xlabel('plasma Lac (μM)');
ylabel('brain Lac');
legend([plot_dmPFC_hdl, plot_aIns_hdl],{'dmPFC','aINS'})
legend('boxoff');
legend_size(pSize);