%% script to check whether subjective stress levels depends on the level of metabolites

%% define all subjects
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% extract stress levels
[preMRS_stress, postMRS_stress,...
    prefMRI_stress, postfMRI_stress,...
    deltaStressPrePostExp] = deal(NaN(1,NS));
% load stress
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    stress_sub_idx = strcmp(sub_nm, excelReadGeneralFile.CID);
    
    preMRS_stress(iS)          = excelReadGeneralFile.StressPr__MRS(stress_sub_idx);
    postMRS_stress(iS)          = excelReadGeneralFile.StressPost_MRS(stress_sub_idx);
    prefMRI_stress(iS)          = excelReadGeneralFile.StressPr__IRMf(stress_sub_idx);
    postfMRI_stress(iS)          = excelReadGeneralFile.StressPost_IRMf(stress_sub_idx);
    deltaStressPrePostExp(iS)  = postfMRI_stress(iS) - preMRS_stress(iS);
end
%% define metabolite and ROI you want to focus on and extract subjects accordingly
[low_met_subs, high_met_subs, mSplit_metabolite_nm] = medSplit_metabolites(study_nm, subject_id);
[metabolite_allSubs, MRS_ROI_nm, linear_metabolite_nm] = metabolite_extraction(study_nm, subject_id);
% remove '_' from name for labels
[mSplit_metabolite_nm_bis] = metab_div_rnm(mSplit_metabolite_nm);
[linear_metabolite_nm_bis] = metab_div_rnm(linear_metabolite_nm);
%% extract the data
preMRS_stress_metSplit.low = preMRS_stress(low_met_subs);
preMRS_stress_metSplit.high = preMRS_stress(high_met_subs);
postMRS_stress_metSplit.low = postMRS_stress(low_met_subs);
postMRS_stress_metSplit.high = postMRS_stress(high_met_subs);
prefMRI_stress_metSplit.low = prefMRI_stress(low_met_subs);
prefMRI_stress_metSplit.high = prefMRI_stress(high_met_subs);
postfMRI_stress_metSplit.low = postfMRI_stress(low_met_subs);
postfMRI_stress_metSplit.high = postfMRI_stress(high_met_subs);
deltaStressPrePostExp_metSplit.low = deltaStressPrePostExp(low_met_subs);
deltaStressPrePostExp_metSplit.high = deltaStressPrePostExp(high_met_subs);

% test difference between two groups
[~,pval.baseline] = ttest2(preMRS_stress_metSplit.low, preMRS_stress_metSplit.high);
[~,pval.pre_vs_post_exp] = ttest2(deltaStressPrePostExp_metSplit.low, deltaStressPrePostExp_metSplit.high);

% test linear correlation with each level of stress
metabolite_ascOrder = sort(metabolite_allSubs(~isnan(metabolite_allSubs)));
% pre-MRS
goodSubs_preMRS = ~isnan(metabolite_allSubs.*preMRS_stress);
[b_preMRS,~,stats_preMRS] = glmfit(metabolite_allSubs(goodSubs_preMRS), preMRS_stress(goodSubs_preMRS),'normal');
preMRS_stress_fit = glmval(b_preMRS,metabolite_ascOrder,'identity');
pval.linearCorrel.preMRS = stats_preMRS.p(2);
% post-MRS
goodSubs_postMRS = ~isnan(metabolite_allSubs.*postMRS_stress);
[b_postMRS,~,stats_postMRS] = glmfit(metabolite_allSubs(goodSubs_postMRS), postMRS_stress(goodSubs_postMRS),'normal');
postMRS_stress_fit = glmval(b_postMRS,metabolite_ascOrder,'identity');
pval.linearCorrel.postMRS = stats_postMRS.p(2);
% pre-fMRI
goodSubs_prefMRI = ~isnan(metabolite_allSubs.*prefMRI_stress);
[b_prefMRI,~,stats_prefMRI] = glmfit(metabolite_allSubs(goodSubs_prefMRI), prefMRI_stress(goodSubs_prefMRI),'normal');
prefMRI_stress_fit = glmval(b_prefMRI,metabolite_ascOrder,'identity');
pval.linearCorrel.prefMRI = stats_prefMRI.p(2);
% post-fMRI
goodSubs_postfMRI = ~isnan(metabolite_allSubs.*postfMRI_stress);
[b_postfMRI,~,stats_postfMRI] = glmfit(metabolite_allSubs(goodSubs_postfMRI), postfMRI_stress(goodSubs_postfMRI),'normal');
postfMRI_stress_fit = glmval(b_postfMRI,metabolite_ascOrder,'identity');
pval.linearCorrel.postfMRI = stats_postfMRI.p(2);
% delta stress end of experiment - baseline
goodSubs_deltaPrePostExp = ~isnan(metabolite_allSubs.*deltaStressPrePostExp);
[b_deltaPrePostExp,~,stats_deltaPrePostExp] = glmfit(metabolite_allSubs(goodSubs_deltaPrePostExp), deltaStressPrePostExp(goodSubs_deltaPrePostExp),'normal');
deltaPrePostExp_stress_fit = glmval(b_deltaPrePostExp,metabolite_ascOrder,'identity');
pval.linearCorrel.deltaPrePostExp = stats_deltaPrePostExp.p(2);

% extract average and sem
[preMRS_stress_metSplit.avg.low,...
    preMRS_stress_metSplit.sem.low] = mean_sem_sd(preMRS_stress_metSplit.low,2);
[preMRS_stress_metSplit.avg.high,...
    preMRS_stress_metSplit.sem.high] = mean_sem_sd(preMRS_stress_metSplit.high,2);
[postMRS_stress_metSplit.avg.low,...
    postMRS_stress_metSplit.sem.low] = mean_sem_sd(postMRS_stress_metSplit.low,2);
[postMRS_stress_metSplit.avg.high,...
    postMRS_stress_metSplit.sem.high] = mean_sem_sd(postMRS_stress_metSplit.high,2);
[prefMRI_stress_metSplit.avg.low,...
    prefMRI_stress_metSplit.sem.low] = mean_sem_sd(prefMRI_stress_metSplit.low,2);
[prefMRI_stress_metSplit.avg.high,...
    prefMRI_stress_metSplit.sem.high] = mean_sem_sd(prefMRI_stress_metSplit.high,2);
[postfMRI_stress_metSplit.avg.low,...
    postfMRI_stress_metSplit.sem.low] = mean_sem_sd(postfMRI_stress_metSplit.low,2);
[postfMRI_stress_metSplit.avg.high,...
    postfMRI_stress_metSplit.sem.high] = mean_sem_sd(postfMRI_stress_metSplit.high,2);
[deltaStressPrePostExp_metSplit.avg.low,...
    deltaStressPrePostExp_metSplit.sem.low] = mean_sem_sd(deltaStressPrePostExp_metSplit.low,2);
[deltaStressPrePostExp_metSplit.avg.high,...
    deltaStressPrePostExp_metSplit.sem.high] = mean_sem_sd(deltaStressPrePostExp_metSplit.high,2);

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
    [preMRS_stress_metSplit.avg.low, postMRS_stress_metSplit.avg.low,...
    prefMRI_stress_metSplit.avg.low, postfMRI_stress_metSplit.avg.low],...
    'FaceColor','b','BarWidth',bWidth);
errorbar([1-bDist, 2-bDist, 3-bDist, 4-bDist],...
    [preMRS_stress_metSplit.avg.low, postMRS_stress_metSplit.avg.low,...
    prefMRI_stress_metSplit.avg.low, postfMRI_stress_metSplit.avg.low],...
    [preMRS_stress_metSplit.sem.low, postMRS_stress_metSplit.sem.low,...
    prefMRI_stress_metSplit.sem.low, postfMRI_stress_metSplit.sem.low],...
    'k','LineStyle','none','LineWidth',3);
% high Ep
bar([1+bDist, 2+bDist, 3+bDist, 4+bDist],...
    [preMRS_stress_metSplit.avg.high, postMRS_stress_metSplit.avg.high,...
    prefMRI_stress_metSplit.avg.high, postfMRI_stress_metSplit.avg.high],...
    'FaceColor','g','BarWidth',bWidth);
errorbar([1+bDist, 2+bDist, 3+bDist, 4+bDist],...
    [preMRS_stress_metSplit.avg.high, postMRS_stress_metSplit.avg.high,...
    prefMRI_stress_metSplit.avg.high, postfMRI_stress_metSplit.avg.high],...
    [preMRS_stress_metSplit.sem.high, postMRS_stress_metSplit.sem.high,...
    prefMRI_stress_metSplit.sem.high, postfMRI_stress_metSplit.sem.high],...
    'k','LineStyle','none','LineWidth',3);

% add infos
xticks([1-bDist, 1+bDist, 2-bDist, 2+bDist, 3-bDist, 3+bDist, 4-bDist, 4+bDist]);
xticklabels({['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm]});
ylabel('subjective stress level');
legend_size(pSize);

%% delta between start and end of experiment
fig;
hold on;
% low Ep
bar(1-bDist, deltaStressPrePostExp_metSplit.avg.low,...
    'FaceColor','b','BarWidth',bWidth);
errorbar(1-bDist,...
    deltaStressPrePostExp_metSplit.avg.low,...
    deltaStressPrePostExp_metSplit.sem.low,...
    'k');
% high Ep
bar(1+bDist, deltaStressPrePostExp_metSplit.avg.high,...
    'FaceColor','g','BarWidth',bWidth);
errorbar(1+bDist,...
    deltaStressPrePostExp_metSplit.avg.high,...
    deltaStressPrePostExp_metSplit.sem.high,...
    'k');

% add infos
xticks([1-bDist, 1+bDist]);
xticklabels({['l',metabolite_nm],['h',metabolite_nm]});
ylabel('subjective stress end - start');
legend_size(pSize);

fig;
hold on;
% low Ep
if sum(low_met_subs) == sum(high_met_subs) % even number
    violinplot([deltaStressPrePostExp_metSplit.low',...
        deltaStressPrePostExp_metSplit.high']);
elseif sum(low_met_subs) == sum(high_met_subs) + 1
    violinplot([deltaStressPrePostExp_metSplit.low',...
        [deltaStressPrePostExp_metSplit.high,NaN]']);
end

% add infos
xticks([1, 2]);
xticklabels({['l',metabolite_nm],['h',metabolite_nm]});
ylabel('subjective stress end - start');
legend_size(pSize);
line(xlim(),[0 0],'LineWidth',1,'Color','k');

%% figure for linear correlation
fig;

% pre-MRS
subplot(2,2,1); hold on;
scat_hdl = scatter(metabolite_allSubs(goodSubs_preMRS),preMRS_stress(goodSubs_preMRS));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(metabolite_ascOrder, preMRS_stress_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel([MRS_ROI_nm,' - ',linear_metabolite_nm_bis]);
ylabel('pre-MRS stress');
legend_size(pSize);

% post-MRS
subplot(2,2,2); hold on;
scat_hdl = scatter(metabolite_allSubs(goodSubs_postMRS),postMRS_stress(goodSubs_postMRS));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(metabolite_ascOrder, postMRS_stress_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel([MRS_ROI_nm,' - ',linear_metabolite_nm_bis]);
ylabel('post-MRS stress');
legend_size(pSize);

% pre-fMRI
subplot(2,2,3); hold on;
scat_hdl = scatter(metabolite_allSubs(goodSubs_prefMRI),prefMRI_stress(goodSubs_prefMRI));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(metabolite_ascOrder, prefMRI_stress_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel([MRS_ROI_nm,' - ',linear_metabolite_nm_bis]);
ylabel('pre-fMRI stress');
legend_size(pSize);

% post-fMRI
subplot(2,2,4); hold on;
scat_hdl = scatter(metabolite_allSubs,postfMRI_stress);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(metabolite_ascOrder, postfMRI_stress_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel([MRS_ROI_nm,' - ',linear_metabolite_nm_bis]);
ylabel('post-fMRI stress');
legend_size(pSize);

%% try to do mixed-effects GLM to check for repeated measures ANOVA (based on chatgpt)
dataTable = table([subject_id';subject_id';subject_id';subject_id'],...
    [ones(NS,1); ones(NS,1)*2; ones(NS,1)*3; ones(NS,1)*4],...
    [preMRS_stress'; postMRS_stress'; prefMRI_stress'; postfMRI_stress'],...
    [metabolite_allSubs';metabolite_allSubs';metabolite_allSubs';metabolite_allSubs';],...
    'VariableNames', {'SubjectID', 'TimePoint', 'SubjectiveStress', linear_metabolite_nm});

% Define the mixed-effects model formula
formula = ['SubjectiveStress ~ 1 + ',linear_metabolite_nm,' + TimePoint + (1|SubjectID)'];

% Fit the mixed-effects model using fitlme
mdl = fitlme(dataTable, formula);

% Display the results
disp(mdl);

% You can also inspect fixed effects and random effects
disp(fixedEffects(mdl));
disp(randomEffects(mdl));
