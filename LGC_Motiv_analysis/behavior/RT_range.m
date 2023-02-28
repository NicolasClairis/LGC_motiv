function[RT_summary] = RT_range()
% [RT_summary] = RT_range()
% RT_range will give you the mean, SEM and SD of the RT (for choices) along
% with the minimum and maximum RT for each task separately. RT_range will
% also extract the mean, SEM and SD of the RT during the performance and
% will compare it with the RT during choices.
%
% OUTPUTS
% RT_summary: structure containing the information detailed above, split
% per task, along with the correlations results

%% if root not defined => ask for it
computerRoot = 'E:';
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% define main variables
nTrialsPerRun = 54;
tasks = {'Ep','Em'};
nTasks = length(tasks);
periods = {'choice','perf'};
nPeriods = length(periods);
for iTask = 1:nTasks
    task_nm = tasks{iTask};
    for iPeriod = 1:nPeriods
        period_nm = periods{iPeriod};
        [RT_summary.(task_nm).mean_RT.(period_nm),...
            RT_summary.(task_nm).median_RT.(period_nm),...
            RT_summary.(task_nm).sd_RT.(period_nm),...
            RT_summary.(task_nm).min_RT.(period_nm),...
            RT_summary.(task_nm).max_RT.(period_nm)] = deal(NaN(1,NS));
    end
end

%% extract the data
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];
    
    %% extract average choice RT across all tasks
    % extract information of runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    nTotalRuns = length(runsStruct.tasks);
    for iTask = 1:nTasks
        task_id = tasks{iTask};
        switch task_id
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end
        runs.(task_id) = strcmp(runsStruct.tasks,task_id);
        nRuns.(task_id) = sum(runs.(task_id));
        nSubjectTotalTrials.(task_id) = nTrialsPerRun*nRuns.(task_id);
        % prepare data to be extracted across runs
        [RT_choice_perTrial_tmp.(task_id),...
            RT_perf_perTrial_tmp.(task_id)] = deal(NaN(nSubjectTotalTrials.(task_id), 1));
        
        jRun = 0;
        for iRun = 1:nTotalRuns
            runToInclude = 0;
            if runs.(task_id)(iRun) == 1
                jRun = jRun + 1;
                runToInclude = 1;
            end
            run_nm = num2str(iRun);
            
            if runToInclude == 1
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);
                
                %% load RT during choice
                [RT_tmp] = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                RT_choice_perTrial_tmp.(task_id)(runTrials_idx) = RT_tmp;
                %% load RT during perf
                switch task_id
                    case 'Ep'
                        latency_tmp = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
                        RT_perf_perTrial_tmp.(task_id)(runTrials_idx) = latency_tmp.allTrials;
                    case 'Em'
                        [~, ~, ~, RT_avg] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
                        RT_perf_perTrial_tmp.(task_id)(runTrials_idx) = RT_avg.allTrials;
                end
            end % run to include
        end % run loop
        
        % average the RT across trials and runs within each subject (for
        % each task separately)
        % choice period
        RT_summary.(task_id).mean_RT.allSubs.choice(iS) = mean(RT_choice_perTrial_tmp.(task_id),1,'omitnan');
        RT_summary.(task_id).median_RT.allSubs.choice(iS) = median(RT_choice_perTrial_tmp.(task_id),1,'omitnan');
        RT_summary.(task_id).sd_RT.allSubs.choice(iS) = std(RT_choice_perTrial_tmp.(task_id),[],1,'omitnan');
        RT_summary.(task_id).min_RT.allSubs.choice(iS) = min(RT_choice_perTrial_tmp.(task_id),[],1,'omitnan');
        RT_summary.(task_id).max_RT.allSubs.choice(iS) = max(RT_choice_perTrial_tmp.(task_id),[],1,'omitnan');
        % perf RT
        RT_summary.(task_id).mean_RT.allSubs.perf(iS) = mean(RT_perf_perTrial_tmp.(task_id),1,'omitnan');
        RT_summary.(task_id).median_RT.allSubs.perf(iS) = median(RT_perf_perTrial_tmp.(task_id),1,'omitnan');
        RT_summary.(task_id).sd_RT.allSubs.perf(iS) = std(RT_perf_perTrial_tmp.(task_id),[],1,'omitnan');
        RT_summary.(task_id).min_RT.allSubs.perf(iS) = min(RT_perf_perTrial_tmp.(task_id),[],1,'omitnan');
        RT_summary.(task_id).max_RT.allSubs.perf(iS) = max(RT_perf_perTrial_tmp.(task_id),[],1,'omitnan');
    end % task loop
end % subject loop

%% average the RT across trials and runs within each subject (for
% each task separately)
fig;
pSize = 30;
for iT = 1:nTasks
     task_nm = tasks{iT};
    % choice period
    RT_summary.(task_nm).mean_RT.meanSubs.choice = mean(RT_summary.(task_nm).mean_RT.allSubs.choice,2,'omitnan');
    RT_summary.(task_nm).median_RT.meanSubs.choice = median(RT_summary.(task_nm).median_RT.allSubs.choice,2,'omitnan');
    RT_summary.(task_nm).sd_RT.meanSubs.choice = std(RT_summary.(task_nm).sd_RT.allSubs.choice,[],2,'omitnan');
    RT_summary.(task_nm).min_RT.meanSubs.choice = min(RT_summary.(task_nm).min_RT.allSubs.choice,[],2,'omitnan');
    RT_summary.(task_nm).max_RT.meanSubs.choice = max(RT_summary.(task_nm).max_RT.allSubs.choice,[],2,'omitnan');
    % perf RT
    RT_summary.(task_nm).mean_RT.meanSubs.perf = mean(RT_summary.(task_nm).mean_RT.allSubs.perf,2,'omitnan');
    RT_summary.(task_nm).median_RT.meanSubs.perf = median(RT_summary.(task_nm).median_RT.allSubs.perf,2,'omitnan');
    RT_summary.(task_nm).sd_RT.meanSubs.perf = std(RT_summary.(task_nm).sd_RT.allSubs.perf,[],2,'omitnan');
    RT_summary.(task_nm).min_RT.meanSubs.perf = min(RT_summary.(task_nm).min_RT.allSubs.perf,[],2,'omitnan');
    RT_summary.(task_nm).max_RT.meanSubs.perf = max(RT_summary.(task_nm).max_RT.allSubs.perf,[],2,'omitnan');
    
    %% test correlation between choice and perf RT
    [RT_summary.(task_nm).rho, RT_summary.(task_nm).pval] = corr(RT_summary.(task_nm).mean_RT.allSubs.perf', RT_summary.(task_nm).mean_RT.allSubs.choice');
    
    subplot(1,nTasks,iT);
    hold on;
    scatter(RT_summary.(task_nm).mean_RT.allSubs.perf',...
        RT_summary.(task_nm).mean_RT.allSubs.choice');
    
    legend_size(pSize);
end
end % function