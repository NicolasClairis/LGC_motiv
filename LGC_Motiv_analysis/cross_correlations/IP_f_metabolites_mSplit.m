%% script to check whether indifference point depends on the level of metabolites

%% working directory
[computerRoot] = LGCM_root_paths();
study_nm = 'study1';
studyFolder = fullfile(computerRoot, study_nm);

%% define all subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection('study1', condition);

%% extract indifference points
[deltaIP.Ep, deltaIP.Em] = deal(NaN(1,NS));
for iS = 1:NS
   sub_nm = subject_id{iS};
   subFolder = [studyFolder,filesep,'CID',sub_nm,filesep,'behavior',filesep];
       IPdata_tmp = getfield(load([subFolder,'delta_IP_CID',sub_nm,'.mat'],'IP_variables'),'IP_variables');
       deltaIP.Ep(iS) = IPdata_tmp.physicalDeltaIP;
       deltaIP.Em(iS) = IPdata_tmp.mentalDeltaIP;
end
%% define metabolite and ROI you want to focus on and extract subjects accordingly
[low_met_subs, high_met_subs, metabolite_nm] = medSplit_metabolites(subject_id);

%% extract the data
tasks = {'Ep','Em'};
for iTask = 1:length(tasks)
    task_nm = tasks{iTask};
    deltaIP_met.(task_nm).all.low = deltaIP.(task_nm)(low_met_subs);
    deltaIP_met.(task_nm).all.high = deltaIP.(task_nm)(high_met_subs);
    
    % test difference between two groups
    [~,pval.(task_nm)] = ttest2(deltaIP_met.(task_nm).all.low, deltaIP_met.(task_nm).all.high);
    
    % extract average and sem
    [deltaIP_met.(task_nm).avg.low,...
        deltaIP_met.(task_nm).sem.low] = mean_sem_sd(deltaIP_met.(task_nm).all.low,2);
    [deltaIP_met.(task_nm).avg.high,...
        deltaIP_met.(task_nm).sem.high] = mean_sem_sd(deltaIP_met.(task_nm).all.high,2);
end

%% display results
% figure parameters
lWidth = 3;
pSize = 30;
bWidth = 0.4;
bDist = 0.2;

% errorbar graph
fig;
hold on;
% low Ep
bar(1-bDist, deltaIP_met.Ep.avg.low,...
    'FaceColor','b','BarWidth',bWidth);
errorbar(1-bDist,...
    deltaIP_met.Ep.avg.low,...
    deltaIP_met.Ep.sem.low,...
    'k');
% high Ep
bar(1+bDist, deltaIP_met.Ep.avg.high,...
    'FaceColor','g','BarWidth',bWidth);
errorbar(1+bDist,...
    deltaIP_met.Ep.avg.high,...
    deltaIP_met.Ep.sem.high,...
    'k');
% low Em
bar(2-bDist, deltaIP_met.Em.avg.low,...
    'FaceColor','b','BarWidth',bWidth);
errorbar(2-bDist,...
    deltaIP_met.Em.avg.low,...
    deltaIP_met.Em.sem.low,...
    'k');
% high Em
bar(2+bDist, deltaIP_met.Em.avg.high,...
    'FaceColor','g','BarWidth',bWidth);
errorbar(2+bDist,...
    deltaIP_met.Em.avg.high,...
    deltaIP_met.Em.sem.high,...
    'k');
% add infos
xticks([1-bDist, 1+bDist, 2-bDist, 2+bDist]);
xticklabels({['l',metabolite_nm],['h',metabolite_nm],...
    ['l',metabolite_nm],['h',metabolite_nm]});
ylabel('delta with Indifference Point (CHF)');
legend_size(pSize);
