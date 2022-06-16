%% script to check whether subjective fatigue levels depends on the level of metabolites

%% working directory
[computerRoot] = LGCM_root_paths();
study_nm = 'study1';
studyFolder = fullfile(computerRoot, study_nm);

%% define all subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection('study1', condition);

%% extract fatigue levels
[preMRS_fatigue, postMRS_fatigue,...
    preIRMf_fatigue, postIRMf_fatigue,...
    deltaFatiguePrePostExp] = deal(NaN(1,NS));
% load fatigue
fatiguePath = fullfile('C:','Users','Loco','Documents','GitHub',...
    'LGC_motiv','LGC_Motiv_analysis','behavior');
fatigue = getfield(load([fatiguePath, filesep, 'fatigue_tmp.mat'],'fatigue'),'fatigue');
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = strcmp(fatigue.subject_id, sub_nm);
    preMRS_fatigue(iS)          = fatigue.preMRS(sub_idx);
    postMRS_fatigue(iS)          = fatigue.postMRS(sub_idx);
    preIRMf_fatigue(iS)          = fatigue.preIRMf(sub_idx);
    postIRMf_fatigue(iS)          = fatigue.postIRMf(sub_idx);
    deltaFatiguePrePostExp(iS)  = fatigue.end_vs_start_exp(sub_idx);
end
%% define metabolite and ROI you want to focus on and extract subjects accordingly
[low_met_subs, high_met_subs, metabolite_nm] = medSplit_metabolites(subject_id);

%% extract the data
preMRS_fatigue_metSplit.low = preMRS_fatigue(low_met_subs);
preMRS_fatigue_metSplit.high = preMRS_fatigue(high_met_subs);
postMRS_fatigue_metSplit.low = postMRS_fatigue(low_met_subs);
postMRS_fatigue_metSplit.high = postMRS_fatigue(high_met_subs);
preIRMf_fatigue_metSplit.low = preIRMf_fatigue(low_met_subs);
preIRMf_fatigue_metSplit.high = preIRMf_fatigue(high_met_subs);
postIRMf_fatigue_metSplit.low = postIRMf_fatigue(low_met_subs);
postIRMf_fatigue_metSplit.high = postIRMf_fatigue(high_met_subs);
deltaFatiguePrePostExp_metSplit.low = deltaFatiguePrePostExp(low_met_subs);
deltaFatiguePrePostExp_metSplit.high = deltaFatiguePrePostExp(high_met_subs);

% test difference between two groups
[~,pval.baseline] = ttest2(preMRS_fatigue_metSplit.low, preMRS_fatigue_metSplit.high);
[~,pval.pre_vs_post_exp] = ttest2(deltaFatiguePrePostExp_metSplit.low, deltaFatiguePrePostExp_metSplit.high);

% extract average and sem
[preMRS_fatigue_metSplit.avg.low,...
    preMRS_fatigue_metSplit.sem.low] = mean_sem_sd(preMRS_fatigue_metSplit.low,2);
[preMRS_fatigue_metSplit.avg.high,...
    preMRS_fatigue_metSplit.sem.high] = mean_sem_sd(preMRS_fatigue_metSplit.high,2);
[postMRS_fatigue_metSplit.avg.low,...
    postMRS_fatigue_metSplit.sem.low] = mean_sem_sd(postMRS_fatigue_metSplit.low,2);
[postMRS_fatigue_metSplit.avg.high,...
    postMRS_fatigue_metSplit.sem.high] = mean_sem_sd(postMRS_fatigue_metSplit.high,2);
[preIRMf_fatigue_metSplit.avg.low,...
    preIRMf_fatigue_metSplit.sem.low] = mean_sem_sd(preIRMf_fatigue_metSplit.low,2);
[preIRMf_fatigue_metSplit.avg.high,...
    preIRMf_fatigue_metSplit.sem.high] = mean_sem_sd(preIRMf_fatigue_metSplit.high,2);
[postIRMf_fatigue_metSplit.avg.low,...
    postIRMf_fatigue_metSplit.sem.low] = mean_sem_sd(postIRMf_fatigue_metSplit.low,2);
[postIRMf_fatigue_metSplit.avg.high,...
    postIRMf_fatigue_metSplit.sem.high] = mean_sem_sd(postIRMf_fatigue_metSplit.high,2);
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

% errorbar graph
fig;
hold on;
% low Ep
bar([1-bDist, 2-bDist, 3-bDist, 4-bDist],...
    [preMRS_fatigue_metSplit.avg.low, postMRS_fatigue_metSplit.avg.low,...
    preIRMf_fatigue_metSplit.avg.low, postIRMf_fatigue_metSplit.avg.low],...
    'FaceColor','b','BarWidth',bWidth);
errorbar([1-bDist, 2-bDist, 3-bDist, 4-bDist],...
    [preMRS_fatigue_metSplit.avg.low, postMRS_fatigue_metSplit.avg.low,...
    preIRMf_fatigue_metSplit.avg.low, postIRMf_fatigue_metSplit.avg.low],...
    [preMRS_fatigue_metSplit.sem.low, postMRS_fatigue_metSplit.sem.low,...
    preIRMf_fatigue_metSplit.sem.low, postIRMf_fatigue_metSplit.sem.low],...
    'k','LineStyle','none','LineWidth',3);
% high Ep
bar([1+bDist, 2+bDist, 3+bDist, 4+bDist],...
    [preMRS_fatigue_metSplit.avg.high, postMRS_fatigue_metSplit.avg.high,...
    preIRMf_fatigue_metSplit.avg.high, postIRMf_fatigue_metSplit.avg.high],...
    'FaceColor','g','BarWidth',bWidth);
errorbar([1+bDist, 2+bDist, 3+bDist, 4+bDist],...
    [preMRS_fatigue_metSplit.avg.high, postMRS_fatigue_metSplit.avg.high,...
    preIRMf_fatigue_metSplit.avg.high, postIRMf_fatigue_metSplit.avg.high],...
    [preMRS_fatigue_metSplit.sem.high, postMRS_fatigue_metSplit.sem.high,...
    preIRMf_fatigue_metSplit.sem.high, postIRMf_fatigue_metSplit.sem.high],...
    'k','LineStyle','none','LineWidth',3);

% add infos
xticks([1-bDist, 1+bDist, 2-bDist, 2+bDist, 3-bDist, 3+bDist, 4-bDist, 4+bDist]);
xticklabels({['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm]});
ylabel('subjective fatigue level');
legend_size(pSize);

% delta between start and end of experiment
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
xticklabels({['l',metabolite_nm],['h',metabolite_nm]});
ylabel('subjective fatigue end - start');
legend_size(pSize);

fig;
hold on;
% low Ep
violinplot([deltaFatiguePrePostExp_metSplit.low', deltaFatiguePrePostExp_metSplit.high'])

% add infos
xticks([1, 2]);
xticklabels({['l',metabolite_nm],['h',metabolite_nm]});
ylabel('subjective fatigue end - start');
legend_size(pSize);
