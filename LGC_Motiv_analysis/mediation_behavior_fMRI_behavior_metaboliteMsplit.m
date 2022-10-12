% script to perform a median split based on some metabolite level and look
% whether the mediation from behavior to fMRI to behavior varies according
% to the levels of the metabolite.


%% first step: select the metabolite and perform the median split on it
% define study and subjects of interest
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);
% define metabolite and ROI you want to focus on
% + extract corresponding metabolites across individuals
[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);
med_metab = median(metabolite_allSubs,2,'omitnan');
metab_idx.low = metabolite_allSubs <= med_metab;
metab_idx.high = metabolite_allSubs > med_metab;

%% extract the mediation for each group
[pval,...
    a_path_allSubs, b_path_allSubs, c_path_allSubs, c_prime_path_allSubs,...
    subject_id,...
    input_f_input_bin, ROI_f_input_bin, output_f_input_bin,...
    ROI_f_ROI_bin, output_f_ROI_bin,...
    input_prm_nm, ROI_short_nm, output_prm_nm] = mediation_behavior_fMRI_behavior(study_nm, condition, subject_id);

mSplit_grps = {'low','high'};
n_grps = length(mSplit_grps);
for iGrp = 1:n_grps
    grp_nm = mSplit_grps{iGrp};
    a_path.(grp_nm) = a_path_allSubs(metab_idx.(grp_nm));
    b_path.(grp_nm) = b_path_allSubs(metab_idx.(grp_nm));
    c_path.(grp_nm) = c_path_allSubs(metab_idx.(grp_nm));
    c_prime_path.(grp_nm) = c_prime_path_allSubs(metab_idx.(grp_nm));
    
    % look at how significant is each path for each group
    [~,pval.(grp_nm).a] = ttest(a_path.(grp_nm));
    [~,pval.(grp_nm).b] = ttest(b_path.(grp_nm));
    [~,pval.(grp_nm).c] = ttest(c_path.(grp_nm));
    [~,pval.(grp_nm).c_prime] = ttest(c_prime_path.(grp_nm));
end

%% compare the paths between the two groups with unpaired t.test
[~,pval.low_vs_high.a] = ttest2(a_path.low, a_path.high);
[~,pval.low_vs_high.b] = ttest2(b_path.low, b_path.high);
[~,pval.low_vs_high.c] = ttest2(c_path.low, c_path.high);
[~,pval.low_vs_high.c_prime] = ttest2(c_prime_path.low, c_prime_path.high);

%% look at average within each group
for iGrp = 1:n_grps
    grp_nm = mSplit_grps{iGrp};
    a_path.mean.(grp_nm) = mean(a_path.(grp_nm),2,'omitnan');
    b_path.mean.(grp_nm) = mean(b_path.(grp_nm),2,'omitnan');
    c_path.mean.(grp_nm) = mean(c_path.(grp_nm),2,'omitnan');
    c_prime_path.mean.(grp_nm) = mean(c_prime_path.(grp_nm),2,'omitnan');
end % loop

%% average data per group

%% figure
