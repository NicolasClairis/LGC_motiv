function[r_corr, pval, mkR, mkP] = compare_kR_vs_kP()
% [r_corr, pval, mkR, mkP] = compare_kR_vs_kP()
% compare_kR_vs_kP will compare kR and kP to see if they are significantly
% correlated and if they are significantly different from each other.
%
% OUTPUTS
% r_corr: structure with correlation coefficient between kR and kP with and
% without removing outliers of each variable
%
% pval: structure with p.values for correlation between kR and kP and for
% paired t.test comparing kR and kP with and without removing outliers of 
% each variable
%
% mkR: structure with reward parameter across subjects, before/after
% removing outliers and with mean values
%
% mkP: punishment parameter across subjects, before/after
% removing outliers and with mean values

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;
%% initialize variables of interest
[kR, kP] = deal(NaN(NS,1));
    
%% load parameters
prm = prm_extraction(study_nm, subject_id,'bayesian','5');
kR(:) = prm.kR;
kP(:) = prm.kP;

%% identify outliers and remove them
[~,~,~,idx_kR_ok] = rmv_outliers_3sd(kR);
[~,~,~,idx_kP_ok] = rmv_outliers_3sd(kP);
idx_ok_subs = (idx_kR_ok.*idx_kP_ok) == 1;

%% perform correlations with and without outliers
[r_corr.allSubs, pval.r_corr.allSubs] = corr(kR, kP);
[r_corr.noOutliers, pval.r_corr.noOutliers] = corr(kR(idx_ok_subs), kP(idx_ok_subs));

%% perform a t.test comparing kR and kP with and without outliers
[~,pval.ttest.allSubs] = ttest(kR, kP);
[~,pval.ttest.noOutliers2] = ttest(kR(idx_ok_subs), kP(idx_ok_subs));

%% extract info regarding kR and kP
% all subjects
mkR.allSubs = kR;
mkP.allSubs = kP;
[mkR.mean_allSubs, mkR.sem_allSubs] = mean_sem_sd(kR, 1);
[mkP.mean_allSubs, mkP.sem_allSubs] = mean_sem_sd(kP, 1);
% subjects after removing outliers of each individual variable
mkR.noOutliers = kR(idx_kR_ok);
mkP.noOutliers = kP(idx_kP_ok);
[mkR.mean_noOutliers, mkR.sem_noOutliers] = mean_sem_sd(kR(idx_kR_ok), 1);
[mkP.mean_noOutliers, mkP.sem_noOutliers] = mean_sem_sd(kP(idx_kP_ok), 1);
% subjects after removing outliers of both kR and kP
mkR.noOutliers2 = kR(idx_ok_subs);
mkP.noOutliers2 = kP(idx_ok_subs);
[mkR.mean_noOutliers2, mkR.sem_noOutliers2] = mean_sem_sd(kR(idx_ok_subs), 1);
[mkP.mean_noOutliers2, mkP.sem_noOutliers2] = mean_sem_sd(kP(idx_ok_subs), 1);

end % function