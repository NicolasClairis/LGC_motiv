% This function will check the correlation between brain metabolites
% measured with 1H-MRS and the fMRI contrast selected

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% load brain metabolites
[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);

%% load fMRI
fMRI_fig_disp = 0;
biasFieldCorr = 0;
[con_vec_all,...
    con_avg, con_sem, con_sd,...
    con_names,...
    ROI_coords, ttest_ROI, subject_id, fMRI_GLM,...
    ~, n_ROIs, ~,...
    con_vec_all_no_outliers,...
    con_avg_no_outliers, con_sem_no_outliers, con_sd_no_outliers,...
    ttest_ROI_no_outliers] = ROI_extraction_group(study_nm, [],...
    subject_id, condition, fMRI_fig_disp, biasFieldCorr);
% select the contrast of interest
[selectedContrast, selectedCon_nm] = fMRI_contrast_selection(con_names);
con_perSub = con_vec_all(selectedContrast, :);
con_short_nm = inputdlg('Can you give a short name for the contrast?');

%% test the correlation
goodSubs = ~isnan(con_perSub.*metabolite_allSubs);
[r_corr, pval] = corr(metabolite_allSubs(goodSubs)', con_perSub(goodSubs)');
% add glm for fit
[~, betas, pval2, con_fitted,...
    mb_sorted, con_fitted_sorted_by_mb] = glm_package(metabolite_allSubs(goodSubs)', con_perSub(goodSubs)', 'normal', 'on');

%% add figure
fig;
scat_hdl = scatter(metabolite_allSubs(goodSubs)', con_perSub(goodSubs)');
scat_hdl_upgrade(scat_hdl);
fit_hdl = plot(mb_sorted, con_fitted_sorted_by_mb);
fit_hdl_upgrade(fit_hdl);
xlabel([MRS_ROI_nm, ' - ',metabolite_nm,' (mM)']);
ylabel(con_short_nm{1});
place_r_and_pval(r_corr, pval);