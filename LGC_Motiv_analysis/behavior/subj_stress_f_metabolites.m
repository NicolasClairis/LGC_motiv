%% script to check whether subjective stress levels depends on the level of metabolites
error('needs to be updated to load summary_participants_infos.xlsx');
%% working directory
[computerRoot] = LGCM_root_paths();
study_nm = 'study1';
studyFolder = fullfile(computerRoot, study_nm);
switch computerRoot
    case 'E:'
        gitPath = fullfile('C:','Users','clairis','Desktop');
    otherwise
        gitPath = fullfile('C:','Users','Loco','Documents');
end
%% define all subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection('study1', condition);

%% extract stress levels
[preMRS_stress, postMRS_stress,...
    preIRMf_stress, postIRMf_stress,...
    deltaStressPrePostExp] = deal(NaN(1,NS));
% load stress
stressPath = fullfile(gitPath,'GitHub',...
    'LGC_motiv','LGC_Motiv_analysis','behavior');
stress = getfield(load([stressPath, filesep, 'stress_tmp.mat'],'stress'),'stress');
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = strcmp(stress.subject_id, sub_nm);
    preMRS_stress(iS)          = stress.preMRS(sub_idx);
    postMRS_stress(iS)          = stress.postMRS(sub_idx);
    preIRMf_stress(iS)          = stress.preIRMf(sub_idx);
    postIRMf_stress(iS)          = stress.postIRMf(sub_idx);
    deltaStressPrePostExp(iS)  = stress.end_vs_start_exp(sub_idx);
end
%% define metabolite and ROI you want to focus on and extract subjects accordingly
[low_met_subs, high_met_subs, metabolite_nm] = medSplit_metabolites(subject_id);

%% extract the data
preMRS_stress_metSplit.low = preMRS_stress(low_met_subs);
preMRS_stress_metSplit.high = preMRS_stress(high_met_subs);
postMRS_stress_metSplit.low = postMRS_stress(low_met_subs);
postMRS_stress_metSplit.high = postMRS_stress(high_met_subs);
preIRMf_stress_metSplit.low = preIRMf_stress(low_met_subs);
preIRMf_stress_metSplit.high = preIRMf_stress(high_met_subs);
postIRMf_stress_metSplit.low = postIRMf_stress(low_met_subs);
postIRMf_stress_metSplit.high = postIRMf_stress(high_met_subs);
deltaStressPrePostExp_metSplit.low = deltaStressPrePostExp(low_met_subs);
deltaStressPrePostExp_metSplit.high = deltaStressPrePostExp(high_met_subs);

% test difference between two groups
[~,pval.baseline] = ttest2(preMRS_stress_metSplit.low, preMRS_stress_metSplit.high);
[~,pval.pre_vs_post_exp] = ttest2(deltaStressPrePostExp_metSplit.low, deltaStressPrePostExp_metSplit.high);

% extract average and sem
[preMRS_stress_metSplit.avg.low,...
    preMRS_stress_metSplit.sem.low] = mean_sem_sd(preMRS_stress_metSplit.low,2);
[preMRS_stress_metSplit.avg.high,...
    preMRS_stress_metSplit.sem.high] = mean_sem_sd(preMRS_stress_metSplit.high,2);
[postMRS_stress_metSplit.avg.low,...
    postMRS_stress_metSplit.sem.low] = mean_sem_sd(postMRS_stress_metSplit.low,2);
[postMRS_stress_metSplit.avg.high,...
    postMRS_stress_metSplit.sem.high] = mean_sem_sd(postMRS_stress_metSplit.high,2);
[preIRMf_stress_metSplit.avg.low,...
    preIRMf_stress_metSplit.sem.low] = mean_sem_sd(preIRMf_stress_metSplit.low,2);
[preIRMf_stress_metSplit.avg.high,...
    preIRMf_stress_metSplit.sem.high] = mean_sem_sd(preIRMf_stress_metSplit.high,2);
[postIRMf_stress_metSplit.avg.low,...
    postIRMf_stress_metSplit.sem.low] = mean_sem_sd(postIRMf_stress_metSplit.low,2);
[postIRMf_stress_metSplit.avg.high,...
    postIRMf_stress_metSplit.sem.high] = mean_sem_sd(postIRMf_stress_metSplit.high,2);
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

% errorbar graph
fig;
hold on;
% low Ep
bar([1-bDist, 2-bDist, 3-bDist, 4-bDist],...
    [preMRS_stress_metSplit.avg.low, postMRS_stress_metSplit.avg.low,...
    preIRMf_stress_metSplit.avg.low, postIRMf_stress_metSplit.avg.low],...
    'FaceColor','b','BarWidth',bWidth);
errorbar([1-bDist, 2-bDist, 3-bDist, 4-bDist],...
    [preMRS_stress_metSplit.avg.low, postMRS_stress_metSplit.avg.low,...
    preIRMf_stress_metSplit.avg.low, postIRMf_stress_metSplit.avg.low],...
    [preMRS_stress_metSplit.sem.low, postMRS_stress_metSplit.sem.low,...
    preIRMf_stress_metSplit.sem.low, postIRMf_stress_metSplit.sem.low],...
    'k','LineStyle','none','LineWidth',3);
% high Ep
bar([1+bDist, 2+bDist, 3+bDist, 4+bDist],...
    [preMRS_stress_metSplit.avg.high, postMRS_stress_metSplit.avg.high,...
    preIRMf_stress_metSplit.avg.high, postIRMf_stress_metSplit.avg.high],...
    'FaceColor','g','BarWidth',bWidth);
errorbar([1+bDist, 2+bDist, 3+bDist, 4+bDist],...
    [preMRS_stress_metSplit.avg.high, postMRS_stress_metSplit.avg.high,...
    preIRMf_stress_metSplit.avg.high, postIRMf_stress_metSplit.avg.high],...
    [preMRS_stress_metSplit.sem.high, postMRS_stress_metSplit.sem.high,...
    preIRMf_stress_metSplit.sem.high, postIRMf_stress_metSplit.sem.high],...
    'k','LineStyle','none','LineWidth',3);

% add infos
xticks([1-bDist, 1+bDist, 2-bDist, 2+bDist, 3-bDist, 3+bDist, 4-bDist, 4+bDist]);
xticklabels({['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm]});
ylabel('subjective stress level');
legend_size(pSize);

% delta between start and end of experiment
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
violinplot([deltaStressPrePostExp_metSplit.low', deltaStressPrePostExp_metSplit.high'])

% add infos
xticks([1, 2]);
xticklabels({['l',metabolite_nm],['h',metabolite_nm]});
ylabel('subjective stress end - start');
legend_size(pSize);
