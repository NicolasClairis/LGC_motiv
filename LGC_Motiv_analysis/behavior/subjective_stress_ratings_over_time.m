% This script will look how subjective stress ratings evolve over time.

%% subject selection
[study_nm, condition, gender, subject_id, NS] = sub_id;

%% load stress ratings
[~,...
    preMRS_stress, postMRS_stress,...
    prefMRI_stress, postfMRI_stress] = extract_subjective_stress_ratings(study_nm, subject_id, NS);
stress = NaN(4, NS);
stress(1,:) = preMRS_stress;
stress(2,:) = postMRS_stress;
stress(3,:) = prefMRI_stress;
stress(4,:) = postfMRI_stress;

beta.global = NaN(1,NS);
goodS = false(1,NS);
for iS = 1:NS
    beta_tmp = glmfit(1:4, stress(:,iS), 'normal');
    beta.global(iS) = beta_tmp(2);
    goodS(iS) = sum(~isnan(stress(:,iS))) == 4;
end % subject loop

% global stats
[beta.mean.global, beta.sem.global] = mean_sem_sd(beta.global,2);
[stress_global.mean,...
    stress_global.sem, stress_global.sd] = mean_sem_sd(stress, 2);
[~,pval.global] = ttest(beta.global);

%% stats
[pval.global, tbl, stats] = anova1(stress');
[comparison,means,h,gnames] = multcompare(stats,'alpha',0.05,'ctype','bonferroni');

% pre-MRS/post-MRS
preMRS_vs_postMRS_idx = comparison(:,1) == 1 & comparison(:,2) == 2;
pval.preMRS_vs_postMRS = comparison(preMRS_vs_postMRS_idx,6);
% pre-MRS/pre-fMRI
preMRS_vs_prefMRI_idx = comparison(:,1) == 1 & comparison(:,2) == 3;
pval.preMRS_vs_prefMRI = comparison(preMRS_vs_prefMRI_idx,6);
% pre-MRS/post-fMRI
preMRS_vs_postfMRI_idx = comparison(:,1) == 1 & comparison(:,2) == 4;
pval.preMRS_vs_postfMRI = comparison(preMRS_vs_postfMRI_idx,6);
% post-MRS/pre-fMRI
postMRS_vs_prefMRI_idx = comparison(:,1) == 2 & comparison(:,2) == 3;
pval.postMRS_vs_prefMRI = comparison(postMRS_vs_prefMRI_idx,6);
% post-MRS/post-fMRI
postMRS_vs_postfMRI_idx = comparison(:,1) == 2 & comparison(:,2) == 4;
pval.postMRS_vs_postfMRI = comparison(postMRS_vs_postfMRI_idx,6);
% pre-fMRI/post-fMRI
prefMRI_vs_postfMRI_idx = comparison(:,1) == 3 & comparison(:,2) == 4;
pval.prefMRI_vs_postfMRI = comparison(prefMRI_vs_postfMRI_idx,6);

%% display figure
pSize = general_fig_prm;
fig;
bar_hdl = bar(1:4,stress_global.mean);
erbar_hdl = errorbar(1:4,stress_global.mean, stress_global.sd);
erbar_hdl.LineStyle = 'none';
erbar_hdl.LineWidth = 3;
erbar_hdl.Color = [0 0 0];
ylim([0 5]);
ylabel('Stress');
xticks(1:4);
xticklabels({'pre-MRS','post-MRS','pre-fMRI','post-fMRI'})
legend_size(pSize);