% load data, remove outliers, and recompute statistics

load('GLM260_dmPFCdACC_60subs.mat')
% [roi_fig] = roi_graph(selectedCon, n_ROIs,...
%     con_vec_all, con_avg, con_sem,...
%     {'lEch','hEch'}, ttest_ROI.p_value);

% remove outliers
rmv_outliers_3sd()