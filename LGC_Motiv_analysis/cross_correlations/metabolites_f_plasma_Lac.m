% correlation between blood Lactate and brain metabolites

%% subject identification
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load plasma Lac
[Lac_struct] = load_plasma_Lac();
plasma_Lac = NaN(1,NS);
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = strcmp(Lac_struct.CID, ['CID',sub_nm]);
    plasma_Lac(iS) = Lac_struct.Lac(sub_idx);
end % subject loop

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac = metabolites.dmPFC.Lac';
aIns_Lac = metabolites.aIns.Lac';

%% check number of subjects
NS_goodS.dmPFC = sum(~isnan(plasma_Lac.*dmPFC_Lac'));
NS_goodS.aINS = sum(~isnan(plasma_Lac.*aIns_Lac'));

%% test significance
[r_corr.dmPFC, betas.dmPFC, pval.dmPFC,...
    ~, Lac_sorted.dmPFC, Lac_fit_xSorted.dmPFC] = glm_package(plasma_Lac', dmPFC_Lac, 'normal', 'on');
[r_corr.aIns, betas.aIns, pval.aIns,...
    ~, Lac_sorted.aIns, Lac_fit_xSorted.aIns] = glm_package(plasma_Lac', aIns_Lac, 'normal', 'on');
[r_corr.dmPFC_vs_aIns, betas.dmPFC_vs_aIns, pval.dmPFC_vs_aIns,...
    ~, Lac_sorted.dmPFC_vs_aIns, Lac_fit_xSorted.dmPFC_vs_aIns] = glm_package(plasma_Lac', dmPFC_Lac-aIns_Lac, 'normal', 'on');
%% figure
[pSize, lWidth, col, mSize] = general_fig_prm;

fig;
% dmPFC
scat_dmPFC_hdl = scatter(plasma_Lac./1000, dmPFC_Lac);
scat_dmPFC_hdl.MarkerEdgeColor = col.black;
scat_dmPFC_hdl.LineWidth = lWidth;
scat_dmPFC_hdl.SizeData = 60;
plot_dmPFC_hdl = plot(Lac_sorted.dmPFC./1000, Lac_fit_xSorted.dmPFC);
plot_dmPFC_hdl.Color = col.black;
plot_dmPFC_hdl.LineWidth = lWidth;
plot_dmPFC_hdl.LineStyle = '-';
% anterior insula
scat_AI_hdl = scatter(plasma_Lac./1000, aIns_Lac);
scat_AI_hdl.MarkerEdgeColor = col.grey;
scat_AI_hdl.LineWidth = lWidth;
scat_AI_hdl.SizeData = 60;
plot_aIns_hdl = plot(Lac_sorted.aIns./1000, Lac_fit_xSorted.aIns);
plot_aIns_hdl.Color = col.grey;
plot_aIns_hdl.LineWidth = lWidth;
plot_aIns_hdl.LineStyle = '-';
xlabel('plasma Lac (mM)');
ylabel('brain Lac (mM)');
legend([plot_dmPFC_hdl, plot_aIns_hdl],{'dmPFC','aINS'})
legend('boxoff');
legend_size(35);