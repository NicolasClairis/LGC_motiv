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
[metabolite_allSubs, ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);
% filter subjects where no data could be obtained
good_subs = ~isnan(metabolite_allSubs);

%% extract the data
tasks = {'Ep','Em'};
nTasks = length(tasks);
for iTask = 1:nTasks
    task_nm = tasks{iTask};
    
    % perform the correlation
    [betas_tmp, ~, stats_tmp] = glmfit(metabolite_allSubs(good_subs),...
        deltaIP.(task_nm)(good_subs), 'normal');
    betas.(task_nm).(['IP_f_',ROI_nm,'_',metabolite_nm]) = betas_tmp;
    pval.(task_nm).(['IP_f_',ROI_nm,'_',metabolite_nm]) = stats_tmp.p;
    fitted_IP.(task_nm) = glmval(betas_tmp, metabolite_allSubs(good_subs), 'identity');
end

%% display results
% figure parameters
lWidth = 3;
pSize = 30;

% errorbar graph
fig;
% loop through tasks
for iTask = 1:nTasks
    task_nm = tasks{iTask};
    subplot(1,2,iTask);
    hold on;
    scatter(metabolite_allSubs(good_subs), deltaIP.(task_nm)(good_subs),...
        'LineWidth',lWidth);
    % re-order fit according to metabolite levels
    [reordered_metabolites, reorder_idx] = sort(metabolite_allSubs(good_subs));
    plot(reordered_metabolites, fitted_IP.(task_nm)(reorder_idx),...
        'LineStyle','--','LineWidth',lWidth,'Color','k');
    xlabel([ROI_nm,'-',metabolite_nm]);
    ylabel('delta with Indifference Point (CHF)');
    legend_size(pSize);
end % task loop
