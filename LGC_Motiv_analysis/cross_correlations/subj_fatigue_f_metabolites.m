% script to check whether subjective fatigue levels depends on the level of metabolites

%% define all subjects
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% extract fatigue levels
[preMRS_fatigue, postMRS_fatigue,...
    prefMRI_fatigue, postfMRI_fatigue,...
    deltaFatiguePrePostExp] = deal(NaN(1,NS));
% load fatigue
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    fatigue_sub_idx = strcmp(sub_nm, excelReadGeneralFile.CID);
    
    preMRS_fatigue(iS)          = excelReadGeneralFile.FatiguePr__MRS(fatigue_sub_idx);
    postMRS_fatigue(iS)          = excelReadGeneralFile.FatiguePost_MRS(fatigue_sub_idx);
    prefMRI_fatigue(iS)          = excelReadGeneralFile.FatiguePr__fMRI(fatigue_sub_idx);
    postfMRI_fatigue(iS)          = excelReadGeneralFile.FatiguePost_fMRI(fatigue_sub_idx);
    deltaFatiguePrePostExp(iS)  = postfMRI_fatigue(iS) - preMRS_fatigue(iS);
end
%% define metabolite and ROI you want to focus on and extract subjects accordingly
[low_met_subs, high_met_subs, mSplit_metabolite_nm] = medSplit_metabolites(study_nm,subject_id);
[metabolite_allSubs, MRS_ROI_nm, linear_metabolite_nm] = metabolite_extraction(study_nm, subject_id);
% remove '_' from name for labels
[mSplit_metabolite_nm_bis] = metab_div_rnm(mSplit_metabolite_nm);
[linear_metabolite_nm_bis] = metab_div_rnm(linear_metabolite_nm);
%% extract the data
preMRS_fatigue_metSplit.low = preMRS_fatigue(low_met_subs);
preMRS_fatigue_metSplit.high = preMRS_fatigue(high_met_subs);
postMRS_fatigue_metSplit.low = postMRS_fatigue(low_met_subs);
postMRS_fatigue_metSplit.high = postMRS_fatigue(high_met_subs);
prefMRI_fatigue_metSplit.low = prefMRI_fatigue(low_met_subs);
prefMRI_fatigue_metSplit.high = prefMRI_fatigue(high_met_subs);
postfMRI_fatigue_metSplit.low = postfMRI_fatigue(low_met_subs);
postfMRI_fatigue_metSplit.high = postfMRI_fatigue(high_met_subs);
deltaFatiguePrePostExp_metSplit.low = deltaFatiguePrePostExp(low_met_subs);
deltaFatiguePrePostExp_metSplit.high = deltaFatiguePrePostExp(high_met_subs);

% test difference between two groups
[~,pval.baseline] = ttest2(preMRS_fatigue_metSplit.low, preMRS_fatigue_metSplit.high);
[~,pval.pre_vs_post_exp] = ttest2(deltaFatiguePrePostExp_metSplit.low, deltaFatiguePrePostExp_metSplit.high);

% test linear correlation with each level of fatigue
metabolite_ascOrder = sort(metabolite_allSubs(~isnan(metabolite_allSubs)));
% pre-MRS
goodSubs_preMRS = ~isnan(metabolite_allSubs.*preMRS_fatigue);
[b_preMRS,~,stats_preMRS] = glmfit(metabolite_allSubs(goodSubs_preMRS), preMRS_fatigue(goodSubs_preMRS),'normal');
preMRS_fatigue_fit = glmval(b_preMRS,metabolite_ascOrder,'identity');
pval.linearCorrel.preMRS = stats_preMRS.p(2);
% post-MRS
goodSubs_postMRS = ~isnan(metabolite_allSubs.*postMRS_fatigue);
[b_postMRS,~,stats_postMRS] = glmfit(metabolite_allSubs(goodSubs_postMRS), postMRS_fatigue(goodSubs_postMRS),'normal');
postMRS_fatigue_fit = glmval(b_postMRS,metabolite_ascOrder,'identity');
pval.linearCorrel.postMRS = stats_postMRS.p(2);
% pre-fMRI
goodSubs_prefMRI = ~isnan(metabolite_allSubs.*prefMRI_fatigue);
[b_prefMRI,~,stats_prefMRI] = glmfit(metabolite_allSubs(goodSubs_prefMRI), prefMRI_fatigue(goodSubs_prefMRI),'normal');
prefMRI_fatigue_fit = glmval(b_prefMRI,metabolite_ascOrder,'identity');
pval.linearCorrel.prefMRI = stats_prefMRI.p(2);
% post-fMRI
goodSubs_postfMRI = ~isnan(metabolite_allSubs.*postfMRI_fatigue);
[b_postfMRI,~,stats_postfMRI] = glmfit(metabolite_allSubs(goodSubs_postfMRI), postfMRI_fatigue(goodSubs_postfMRI),'normal');
postfMRI_fatigue_fit = glmval(b_postfMRI,metabolite_ascOrder,'identity');
pval.linearCorrel.postfMRI = stats_postfMRI.p(2);
% delta fatigue end of experiment - baseline
goodSubs_deltaPrePostExp = ~isnan(metabolite_allSubs.*deltaFatiguePrePostExp);
[b_deltaPrePostExp,~,stats_deltaPrePostExp] = glmfit(metabolite_allSubs(goodSubs_deltaPrePostExp), deltaFatiguePrePostExp(goodSubs_deltaPrePostExp),'normal');
deltaPrePostExp_fatigue_fit = glmval(b_deltaPrePostExp,metabolite_ascOrder,'identity');
pval.linearCorrel.deltaPrePostExp = stats_deltaPrePostExp.p(2);

% extract average and sem
[preMRS_fatigue_metSplit.avg.low,...
    preMRS_fatigue_metSplit.sem.low] = mean_sem_sd(preMRS_fatigue_metSplit.low,2);
[preMRS_fatigue_metSplit.avg.high,...
    preMRS_fatigue_metSplit.sem.high] = mean_sem_sd(preMRS_fatigue_metSplit.high,2);
[postMRS_fatigue_metSplit.avg.low,...
    postMRS_fatigue_metSplit.sem.low] = mean_sem_sd(postMRS_fatigue_metSplit.low,2);
[postMRS_fatigue_metSplit.avg.high,...
    postMRS_fatigue_metSplit.sem.high] = mean_sem_sd(postMRS_fatigue_metSplit.high,2);
[prefMRI_fatigue_metSplit.avg.low,...
    prefMRI_fatigue_metSplit.sem.low] = mean_sem_sd(prefMRI_fatigue_metSplit.low,2);
[prefMRI_fatigue_metSplit.avg.high,...
    prefMRI_fatigue_metSplit.sem.high] = mean_sem_sd(prefMRI_fatigue_metSplit.high,2);
[postfMRI_fatigue_metSplit.avg.low,...
    postfMRI_fatigue_metSplit.sem.low] = mean_sem_sd(postfMRI_fatigue_metSplit.low,2);
[postfMRI_fatigue_metSplit.avg.high,...
    postfMRI_fatigue_metSplit.sem.high] = mean_sem_sd(postfMRI_fatigue_metSplit.high,2);
[deltaFatiguePrePostExp_metSplit.avg.low,...
    deltaFatiguePrePostExp_metSplit.sem.low] = mean_sem_sd(deltaFatiguePrePostExp_metSplit.low,2);
[deltaFatiguePrePostExp_metSplit.avg.high,...
    deltaFatiguePrePostExp_metSplit.sem.high] = mean_sem_sd(deltaFatiguePrePostExp_metSplit.high,2);

%% display results
% figure parameters
lWidth = 3;
pSize = 25;
bWidth = 0.4;
bDist = 0.2;
grey = [143 143 143]./255;

%% errorbar graph
fig;
hold on;
% low Ep
bar([1-bDist, 2-bDist, 3-bDist, 4-bDist],...
    [preMRS_fatigue_metSplit.avg.low, postMRS_fatigue_metSplit.avg.low,...
    prefMRI_fatigue_metSplit.avg.low, postfMRI_fatigue_metSplit.avg.low],...
    'FaceColor','b','BarWidth',bWidth);
errorbar([1-bDist, 2-bDist, 3-bDist, 4-bDist],...
    [preMRS_fatigue_metSplit.avg.low, postMRS_fatigue_metSplit.avg.low,...
    prefMRI_fatigue_metSplit.avg.low, postfMRI_fatigue_metSplit.avg.low],...
    [preMRS_fatigue_metSplit.sem.low, postMRS_fatigue_metSplit.sem.low,...
    prefMRI_fatigue_metSplit.sem.low, postfMRI_fatigue_metSplit.sem.low],...
    'k','LineStyle','none','LineWidth',3);
% high Ep
bar([1+bDist, 2+bDist, 3+bDist, 4+bDist],...
    [preMRS_fatigue_metSplit.avg.high, postMRS_fatigue_metSplit.avg.high,...
    prefMRI_fatigue_metSplit.avg.high, postfMRI_fatigue_metSplit.avg.high],...
    'FaceColor','g','BarWidth',bWidth);
errorbar([1+bDist, 2+bDist, 3+bDist, 4+bDist],...
    [preMRS_fatigue_metSplit.avg.high, postMRS_fatigue_metSplit.avg.high,...
    prefMRI_fatigue_metSplit.avg.high, postfMRI_fatigue_metSplit.avg.high],...
    [preMRS_fatigue_metSplit.sem.high, postMRS_fatigue_metSplit.sem.high,...
    prefMRI_fatigue_metSplit.sem.high, postfMRI_fatigue_metSplit.sem.high],...
    'k','LineStyle','none','LineWidth',3);

% add infos
xticks([1-bDist, 1+bDist, 2-bDist, 2+bDist, 3-bDist, 3+bDist, 4-bDist, 4+bDist]);
xticklabels({['l',mSplit_metabolite_nm_bis],['h',mSplit_metabolite_nm_bis],...
    ['l',mSplit_metabolite_nm_bis],['h',mSplit_metabolite_nm_bis],...
    ['l',mSplit_metabolite_nm_bis],['h',mSplit_metabolite_nm_bis],...
    ['l',mSplit_metabolite_nm_bis],['h',mSplit_metabolite_nm_bis]});
ylabel('subjective fatigue level');
legend_size(pSize);

%% delta between start and end of experiment
fig;
hold on;
% low Ep
bar(1-bDist, deltaFatiguePrePostExp_metSplit.avg.low,...
    'FaceColor','b','BarWidth',bWidth);
errorbar(1-bDist,...
    deltaFatiguePrePostExp_metSplit.avg.low,...
    deltaFatiguePrePostExp_metSplit.sem.low,...
    'k');
% high Ep
bar(1+bDist, deltaFatiguePrePostExp_metSplit.avg.high,...
    'FaceColor','g','BarWidth',bWidth);
errorbar(1+bDist,...
    deltaFatiguePrePostExp_metSplit.avg.high,...
    deltaFatiguePrePostExp_metSplit.sem.high,...
    'k');

% add infos
xticks([1-bDist, 1+bDist]);
xticklabels({['l',mSplit_metabolite_nm_bis],['h',mSplit_metabolite_nm_bis]});
ylabel('subjective fatigue end - start');
legend_size(pSize);

fig;
hold on;
% low Ep
violinplot([deltaFatiguePrePostExp_metSplit.low', deltaFatiguePrePostExp_metSplit.high'])

% add infos
xticks([1, 2]);
xticklabels({['l',mSplit_metabolite_nm_bis],['h',mSplit_metabolite_nm_bis]});
ylabel('subjective fatigue end - start');
legend_size(pSize);

%% figure for linear correlation
fig;

% pre-MRS
subplot(2,2,1); hold on;
scat_hdl = scatter(metabolite_allSubs(goodSubs_preMRS),preMRS_fatigue(goodSubs_preMRS));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(metabolite_ascOrder, preMRS_fatigue_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel([MRS_ROI_nm,' - ',linear_metabolite_nm_bis]);
ylabel('pre-MRS fatigue');
legend_size(pSize);

% post-MRS
subplot(2,2,2); hold on;
scat_hdl = scatter(metabolite_allSubs(goodSubs_postMRS),postMRS_fatigue(goodSubs_postMRS));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(metabolite_ascOrder, postMRS_fatigue_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel([MRS_ROI_nm,' - ',linear_metabolite_nm_bis]);
ylabel('post-MRS fatigue');
legend_size(pSize);

% pre-fMRI
subplot(2,2,3); hold on;
scat_hdl = scatter(metabolite_allSubs(goodSubs_prefMRI),prefMRI_fatigue(goodSubs_prefMRI));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(metabolite_ascOrder, prefMRI_fatigue_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel([MRS_ROI_nm,' - ',linear_metabolite_nm_bis]);
ylabel('pre-fMRI fatigue');
legend_size(pSize);

% post-fMRI
subplot(2,2,4); hold on;
scat_hdl = scatter(metabolite_allSubs,postfMRI_fatigue);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(metabolite_ascOrder, postfMRI_fatigue_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel([MRS_ROI_nm,' - ',linear_metabolite_nm_bis]);
ylabel('post-fMRI fatigue');
legend_size(pSize);

%% figure for delta as well
fig;
scat_hdl = scatter(metabolite_allSubs(goodSubs_deltaPrePostExp),deltaFatiguePrePostExp(goodSubs_deltaPrePostExp));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(metabolite_ascOrder, deltaPrePostExp_fatigue_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel([MRS_ROI_nm,' - ',linear_metabolite_nm_bis]);
ylabel('subjective fatigue end - start');
legend_size(pSize);