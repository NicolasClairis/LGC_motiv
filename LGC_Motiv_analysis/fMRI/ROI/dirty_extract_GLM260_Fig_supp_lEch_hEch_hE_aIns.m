% load data, remove outliers, and recompute statistics

%% load data
folderPath =  'P:\boulot\postdoc_CarmenSandi\papers\Clairis_mediation_Lac\figures\figS6_dmPFC_aIns_intraindiv';
dmPFC_loadStruct = load([folderPath,filesep,'GLM260_aIns_60subs.mat']);
con_vec_all1 = dmPFC_loadStruct.con_vec_all;
con_avg1 = dmPFC_loadStruct.con_avg;
con_sem1 = dmPFC_loadStruct.con_sem;
selectedCon = dmPFC_loadStruct.selectedCon;
ttest_pval1 = dmPFC_loadStruct.ttest_ROI.p_value;
con_names = dmPFC_loadStruct.con_names;
n_ROIs = dmPFC_loadStruct.n_ROIs;

con_nms = {'lEch','hEch'};
% % display graph before outlier removal
% [roi_fig] = roi_graph(selectedCon, n_ROIs,...
%     con_vec_all1, con_avg1, con_sem1,...
%     con_nms, ttest_pval1);
pval.raw.lEch = ttest_pval1(selectedCon(1));
pval.raw.hEch = ttest_pval1(selectedCon(2));
% compare the two:
[~,pval.raw.hEch_vs_lEch] = ttest(con_vec_all1(selectedCon(1),:), con_vec_all1(selectedCon(2),:));

%% remove outliers and recompute stats
lEch_data = con_vec_all1(selectedCon(1),:);
hEch_data = con_vec_all1(selectedCon(2),:);
% remove outlier
% [~, idx_badSubs_lEch] = rmv_outliers_3sd(lEch_data);
idx_badSubs_lEch = lEch_data < -30; % to remove 2 weird outliers
[~, idx_badSubs_hEch] = rmv_outliers_3sd(hEch_data);
con_vec_all2 = con_vec_all1;
con_vec_all2(:,(idx_badSubs_lEch | idx_badSubs_hEch)) = NaN;
[con_avg2, con_sem2] = mean_sem_sd(con_vec_all2,2);
% perform stats
[~,pval.filtered.lEch] = ttest(con_vec_all2(selectedCon(1),:));
[~,pval.filtered.hEch] = ttest(con_vec_all2(selectedCon(2),:));
% compare the two:
[~,pval.filtered.hEch_vs_lEch] = ttest(con_vec_all2(selectedCon(1),:), con_vec_all2(selectedCon(2),:));

% do all tests
nCons = size(con_vec_all2,1);
ttest_pval2 = NaN(nCons,1);
for iC = 1:nCons
   [~,ttest_pval2(iC)] = ttest(con_vec_all2(iC,:));
end

% display graph before outlier removal
[roi_fig] = roi_graph(selectedCon, n_ROIs,...
    con_vec_all2, con_avg2, con_sem2,...
    con_nms, ttest_pval2);
ylabel('hE regression estimate');