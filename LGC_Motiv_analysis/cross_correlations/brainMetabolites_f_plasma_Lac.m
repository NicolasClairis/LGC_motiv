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
goodS.dmPFC = ~isnan(plasma_Lac.*dmPFC_Lac');
NS_goodS.dmPFC = sum(goodS.dmPFC);
goodS.aINS = ~isnan(plasma_Lac.*aIns_Lac');
NS_goodS.aINS = sum(goodS.aINS);
goodS.all = ~isnan(plasma_Lac.*aIns_Lac'.*dmPFC_Lac');
NS_goodS.all = sum(goodS.all);

%% test significance
[r_corr.dmPFC, betas.dmPFC, pval.dmPFC,...
    ~, Lac_sorted.dmPFC, Lac_fit_xSorted.dmPFC] = glm_package(plasma_Lac', dmPFC_Lac, 'normal', 'on');
[r_corr.aIns, betas.aIns, pval.aIns,...
    ~, Lac_sorted.aIns, Lac_fit_xSorted.aIns] = glm_package(plasma_Lac', aIns_Lac, 'normal', 'on');
[r_corr.dmPFC_vs_aIns, betas.dmPFC_vs_aIns, pval.dmPFC_vs_aIns,...
    ~, Lac_sorted.dmPFC_vs_aIns, Lac_fit_xSorted.dmPFC_vs_aIns] = glm_package(plasma_Lac', dmPFC_Lac-aIns_Lac, 'normal', 'on');

%% redo tests but restricting the analysis to the common subjects
% and extract relevant material for corcor comparison of correlation
% coefficients
[r_corr.overlap.dmPFC] = corr(plasma_Lac(goodS.all)', dmPFC_Lac(goodS.all));
[r_corr.overlap.aIns] = corr(plasma_Lac(goodS.all)', aIns_Lac(goodS.all));
[r_corr.overlap.dmPFC_vs_aIns] = corr(dmPFC_Lac(goodS.all), aIns_Lac(goodS.all));
% display relevant information for corcor:
disp('relevant information for corcor R package:');
disp(['r(plasma-Lac/dmPFC-dACC-Lac = ',num2str(r_corr.overlap.dmPFC)]);
disp(['r(plasma-Lac/aIns-Lac = ',num2str(r_corr.overlap.aIns)]);
disp(['r(dmPFC/dACC-Lac/aIns-Lac = ',num2str(r_corr.overlap.dmPFC_vs_aIns)]);
disp(['number of subjects ok for both dmPFC/dACC, aIns and plasma lactate = ',num2str(NS_goodS.all)]);

%% figure
[pSize, lWidth, col, mSize] = general_fig_prm;
dmPFC_col = col.red;
aINS_col = col.blue_light;

%% both dmPFC/dACC and aIns in same graph
fig;
% dmPFC/dACC
scat_dmPFC_hdl = scatter(plasma_Lac./1000, dmPFC_Lac);
scat_dmPFC_hdl.MarkerEdgeColor = dmPFC_col;
scat_dmPFC_hdl.LineWidth = lWidth;
scat_dmPFC_hdl.SizeData = 60;
plot_dmPFC_hdl = plot(Lac_sorted.dmPFC./1000, Lac_fit_xSorted.dmPFC);
plot_dmPFC_hdl.Color = dmPFC_col;
plot_dmPFC_hdl.LineWidth = lWidth;
plot_dmPFC_hdl.LineStyle = '-';
% anterior insula
scat_AI_hdl = scatter(plasma_Lac./1000, aIns_Lac);
scat_AI_hdl.MarkerEdgeColor = aINS_col;
scat_AI_hdl.LineWidth = lWidth;
scat_AI_hdl.SizeData = 60;
plot_aIns_hdl = plot(Lac_sorted.aIns./1000, Lac_fit_xSorted.aIns);
plot_aIns_hdl.Color = aINS_col;
plot_aIns_hdl.LineWidth = lWidth;
plot_aIns_hdl.LineStyle = '-';
xlabel('plasma Lac (mM)');
ylabel('brain Lac (mM)');
legend([plot_dmPFC_hdl, plot_aIns_hdl],{'dmPFC','aINS'})
legend('boxoff');
legend_size(35);

%% split dmPFC/dACC and aIns in 2 separate graphs
% dmPFC/dACC
fig;
scat_dmPFC_hdl = scatter(plasma_Lac./1000, dmPFC_Lac);
fit_dmPFC_hdl = plot(Lac_sorted.dmPFC./1000, Lac_fit_xSorted.dmPFC);
scat_hdl_upgrade(scat_dmPFC_hdl);
fit_hdl_upgrade(fit_dmPFC_hdl);
xlabel('plasma lactate (mM)');
ylabel('dmPFC/dACC lactate (mM)');
place_r_and_pval(r_corr.dmPFC, pval.dmPFC(2));
legend_size(pSize);

% anterior insula
fig;
scat_AI_hdl = scatter(plasma_Lac./1000, aIns_Lac);
fit_aIns_hdl = plot(Lac_sorted.aIns./1000, Lac_fit_xSorted.aIns);
scat_hdl_upgrade(scat_AI_hdl);
fit_hdl_upgrade(fit_aIns_hdl);
xlabel('plasma lactate (mM)');
ylabel('aIns lactate (mM)');
place_r_and_pval(r_corr.aIns, pval.aIns(2));
legend_size(pSize);