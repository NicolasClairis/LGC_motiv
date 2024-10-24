% This script will look how subjective fatigue ratings evolve over time.

%% subject selection
[study_nm, condition, gender, subject_id, NS] = sub_id;

%% load fatigue ratings
[~,...
    preMRS_fatigue, postMRS_fatigue,...
    prefMRI_fatigue, postfMRI_fatigue] = extract_subjective_fatigue_ratings(study_nm, subject_id, NS);
fatigue = NaN(4, NS);
fatigue(1,:) = preMRS_fatigue;
fatigue(2,:) = postMRS_fatigue;
fatigue(3,:) = prefMRI_fatigue;
fatigue(4,:) = postfMRI_fatigue;

beta.global = NaN(1,NS);
goodS = false(1,NS);
for iS = 1:NS
    beta_tmp = glmfit(1:4, fatigue(:,iS), 'normal');
    beta.global(iS) = beta_tmp(2);
    goodS(iS) = sum(~isnan(fatigue(:,iS))) == 4;
end % subject loop

% global stats
[beta.mean.global, beta.sem.global] = mean_sem_sd(beta.global,2);
[fatigue_global.mean,...
    fatigue_global.sem, fatigue_global.sd] = mean_sem_sd(fatigue, 2);
[~,pval.global] = ttest(beta.global);

%% stats
[pval.global, tbl, stats] = anova1(fatigue');
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
bar_hdl = bar(1:4,fatigue_global.mean);
erbar_hdl = errorbar(1:4,fatigue_global.mean, fatigue_global.sd);
erbar_hdl.LineStyle = 'none';
erbar_hdl.LineWidth = 3;
erbar_hdl.Color = [0 0 0];
ylim([0 8]);
ylabel('Fatigue');
xticks(1:4);
xticklabels({'pre-MRS','post-MRS','pre-fMRI','post-fMRI'})
legend_size(pSize);