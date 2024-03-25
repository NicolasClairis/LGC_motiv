function[RT_summary] = RT_range(subject_id, condition, fig_disp)
% [RT_summary] = RT_range(subject_id, condition, fig_disp)
% RT_range will give you the mean, SEM and SD of the RT (for choices) along
% with the minimum and maximum RT for each task separately. RT_range will
% also extract the mean, SEM and SD of the RT during the performance and
% will compare it with the RT during choices.
%
%
% INPUTS
% subject_id: list of subjects to include (asked if left empty)
%
% condition: condition to use to know which runs to include (asked if left
% empty)
%
% fig_disp: display final figure comparing perf and choice RT (1)
% or not (0)?
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
% condition selection
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
% subject selection
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm,condition);
else
    NS = length(subject_id);
end

% store data in output
RT_summary.subject_id = subject_id;
RT_summary.condition = condition;

%% figure display by default
if ~exist('fig_disp','var') || isempty(fig_disp) || ~ismember(fig_disp,[0,1])
    fig_disp = 1;
end

%% define main variables
nTrialsPerRun = 54;
tasks = {'EpEm','Ep','Em'};
nTasks = length(tasks);
tasks_bis = {'Ep','Em'};
nTasks_bis = length(tasks_bis);
periods = {'choice','perf'};
nPeriods = length(periods);
for iTask = 1:nTasks
    task_nm = tasks{iTask};
    for iPeriod = 1:nPeriods
        period_nm = periods{iPeriod};
        [RT_summary.(task_nm).mean_RT.allSubs.(period_nm),...
            RT_summary.(task_nm).median_RT.allSubs.(period_nm),...
            RT_summary.(task_nm).sd_RT.allSubs.(period_nm),...
            RT_summary.(task_nm).min_RT.allSubs.(period_nm),...
            RT_summary.(task_nm).max_RT.allSubs.(period_nm)] = deal(NaN(1,NS));
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
    for iTask = 1:nTasks_bis
        task_id = tasks_bis{iTask};
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
    
    % average RT across the two tasks
    % choice period
    RT_summary.EpEm.mean_RT.allSubs.choice(iS) = mean([RT_choice_perTrial_tmp.Ep; RT_choice_perTrial_tmp.Em],1,'omitnan');
    RT_summary.EpEm.median_RT.allSubs.choice(iS) = median([RT_choice_perTrial_tmp.Ep; RT_choice_perTrial_tmp.Em],1,'omitnan');
    RT_summary.EpEm.sd_RT.allSubs.choice(iS) = std([RT_choice_perTrial_tmp.Ep; RT_choice_perTrial_tmp.Em],[],1,'omitnan');
    RT_summary.EpEm.min_RT.allSubs.choice(iS) = min([RT_choice_perTrial_tmp.Ep; RT_choice_perTrial_tmp.Em],[],1,'omitnan');
    RT_summary.EpEm.max_RT.allSubs.choice(iS) = max([RT_choice_perTrial_tmp.Ep; RT_choice_perTrial_tmp.Em],[],1,'omitnan');
    % perf RT
    RT_summary.EpEm.mean_RT.allSubs.perf(iS) = mean([RT_perf_perTrial_tmp.Ep; RT_perf_perTrial_tmp.Em],1,'omitnan');
    RT_summary.EpEm.median_RT.allSubs.perf(iS) = median([RT_perf_perTrial_tmp.Ep; RT_perf_perTrial_tmp.Em],1,'omitnan');
    RT_summary.EpEm.sd_RT.allSubs.perf(iS) = std([RT_perf_perTrial_tmp.Ep; RT_perf_perTrial_tmp.Em],[],1,'omitnan');
    RT_summary.EpEm.min_RT.allSubs.perf(iS) = min([RT_perf_perTrial_tmp.Ep; RT_perf_perTrial_tmp.Em],[],1,'omitnan');
    RT_summary.EpEm.max_RT.allSubs.perf(iS) = max([RT_perf_perTrial_tmp.Ep; RT_perf_perTrial_tmp.Em],[],1,'omitnan');
end % subject loop

%% average the RT across trials and runs within each subject (for
% each task separately)
if fig_disp == 1
    fig;
    pSize = 30;
end
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
    [RT_summary.(task_nm).rho,...
        RT_summary.(task_nm).pval.corr] = corr(RT_summary.(task_nm).mean_RT.allSubs.perf', RT_summary.(task_nm).mean_RT.allSubs.choice');
    [~, RT_summary.(task_nm).betas, RT_summary.(task_nm).pval.glm, ~, RT_perf_tmp, RT_choice_tmp] = glm_package(RT_summary.(task_nm).mean_RT.allSubs.perf',...
        RT_summary.(task_nm).mean_RT.allSubs.choice', 'normal', 'on');
    
    if fig_disp == 1
        subplot(1,nTasks,iT);
        switch task_nm
            case 'EpEm'
                title('General');
            case 'Ep'
                title('Physical');
            case 'Em'
                title('Mental');
        end
        hold on;
        scat_hdl = scatter(RT_summary.(task_nm).mean_RT.allSubs.perf',...
            RT_summary.(task_nm).mean_RT.allSubs.choice');
        scat_hdl_upgrade(scat_hdl);
        fit_hdl = plot(RT_perf_tmp, RT_choice_tmp);
        fit_hdl_upgrade(fit_hdl);
        
        legend_size(pSize);
        xlabel('RT performance (s)');
        ylabel('RT choice (s)');
        place_r_and_pval(RT_summary.(task_nm).rho,...
            RT_summary.(task_nm).pval.corr);
    end % display figure
end % loop over physical/mental tasks
end % function