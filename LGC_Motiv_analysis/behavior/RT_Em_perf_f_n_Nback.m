function[RT_avg, RT_avg_sd]=RT_Em_perf_f_n_Nback()
% [RT_avg, RT_avg_sd]=RT_Em_perf_f_n_Nback()
% RT_Em_perf_f_n_Nback will extract the average RT per subject and across
% subjects depending on where people are in the N-back task. The idea is to
% have an idea towards the within-trial variability of RT to see whether
% averaging per trial is meaningful.
%
% OUTPUTS
% RT_avg: reaction times average per subject and across subjects per
% response number
%
% RT_avg_sd: average standard deviation of the reaction times per trial per effort
% level per subject.

%% subject identification
[study_nm, ~, ~, subject_id, NS] = sub_id;
%% working directories
% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% main parameters
task_id = 'Em';
task_fullName = 'mental';
nTrialsPerRun = 54;
n_Em_runs = 2;
nTotalEmTrials = nTrialsPerRun*n_Em_runs;
Ech_levels = 0:3;
n_E_ch = length(Ech_levels);
[RT_avg.perSub.perEch.with2first, RT_sd.perSub.perEch.with2first,...
    RT_avg.perSub.perEch.pureNback, RT_sd.perSub.perEch.pureNback] = deal(NaN(n_E_ch, NS));
n_potential_responses = 35; % high number on purpose
[RT_avg.perSub.per_n_Nback, RT_sd.perSub.per_n_Nback] = deal(NaN(n_potential_responses, NS));
[RT_avg.perSub.perEch_per_n_Nback,...
    RT_sd.perSub.perEch_per_n_Nback] = deal(NaN(n_E_ch, n_potential_responses, NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];
    %% extract average choice RT across all tasks
    % extract information of runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    nTotalRuns = length(runsStruct.tasks);
    runs = strcmp(runsStruct.tasks, task_id);
%     nEmRuns = sum(runs);
%     nSubjectTotalTrials = nTrialsPerRun*nEmRuns;
    [RT_avg_run_with2first_tmp,...
        RT_sd_run_with2first_tmp,...
        RT_avg_run_pureNback_tmp,...
        RT_sd_run_pureNback_tmp,...
        Ech_run_tmp] = deal(NaN(nTotalEmTrials,1));
    RT_per_n_Nback_tmp = NaN(nTotalEmTrials, n_potential_responses);
    jRun = 0;
    for iRun = 1:nTotalRuns
        runToInclude = 0;
        if runs(iRun) == 1
            jRun = jRun + 1;
            runToInclude = 1;
        end

        if runToInclude == 1
            runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);

            %% load data
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                '_task.mat']);
            mentalE_perf = behaviorStruct_tmp.mentalE_perf;
            
            for iTrial = 1:nTrialsPerRun
                jTrial = iTrial + nTrialsPerRun*(jRun - 1);
                rt_tmp = mentalE_perf.perfSummary{1,iTrial}.rt;
                % average RT including 2 first answers (which are just
                % for initialization)
                RT_avg_run_with2first_tmp(jTrial) = mean(rt_tmp, 2, 'omitnan');
                RT_sd_run_with2first_tmp(jTrial) = std(rt_tmp, 0, 2, 'omitnan');
                % average RT for N-back answers only
                RT_avg_run_pureNback_tmp(jTrial) = mean(rt_tmp(3:end), 2, 'omitnan');
                RT_sd_run_pureNback_tmp(jTrial) = std(rt_tmp(3:end), 0, 2, 'omitnan');
                % effort chosen for current trial
                Ech_run_tmp(jTrial) = mentalE_perf.E_chosen(iTrial);
                % RT per answer
                RT_per_n_Nback_tmp(jTrial, 1:length(rt_tmp)) = rt_tmp;
                if length(rt_tmp) > n_potential_responses
                    error(['need to increase n_potential_responses to ',...
                        'at least n = ',num2str(length(rt_tmp))]);
                end
            end % trial loop
        end % run to include?
    end % run loop
    %% average data per subject
    for iEch = Ech_levels
        jEch = iEch + 1;
        Ech_idx_tmp = Ech_run_tmp == iEch;
        % average and SD per effort level including 2 first answers
        RT_avg.perSub.perEch.with2first(jEch, iS) = mean(RT_avg_run_with2first_tmp(Ech_idx_tmp), 1, 'omitnan');
        RT_sd.perSub.perEch.with2first(jEch, iS) = std(RT_avg_run_with2first_tmp(Ech_idx_tmp), 0, 1, 'omitnan');
        % average and SD per effort level only for Nback replies
        RT_avg.perSub.perEch.pureNback(jEch, iS) = mean(RT_avg_run_pureNback_tmp(Ech_idx_tmp), 1, 'omitnan');
        RT_sd.perSub.perEch.pureNback(jEch, iS) = std(RT_avg_run_pureNback_tmp(Ech_idx_tmp), 0, 1, 'omitnan');
        % average and SD per effort level and answer
        RT_avg.perSub.perEch_per_n_Nback(jEch,:,iS) = mean(RT_per_n_Nback_tmp(Ech_idx_tmp,:), 1, 'omitnan');
        RT_sd.perSub.perEch_per_n_Nback(jEch,:,iS) = std(RT_per_n_Nback_tmp(Ech_idx_tmp,:), 0, 1, 'omitnan');
    end
    % average and SD per answer
    RT_avg.perSub.per_n_Nback(:,iS) = mean(RT_per_n_Nback_tmp, 1, 'omitnan');
    RT_sd.perSub.per_n_Nback(:,iS) = std(RT_per_n_Nback_tmp, 0, 1, 'omitnan');
end % subject loop

% average and within-subject SD per effort level including 2 first answers
[RT_avg.aSubs.perEch.with2first,...
    RT_sem.aSubs.perEch.with2first] = mean_sem_sd(RT_avg.perSub.perEch.with2first, 2);
[RT_avg_sd.aSubs.perEch.with2first,...
    RT_sem_sd.aSubs.perEch.with2first]  = mean_sem_sd(RT_sd.perSub.perEch.with2first, 2);
% average and within-subject SD per effort level only for Nback replies
[RT_avg.aSubs.perEch.pureNback,...
    RT_sem.aSubs.perEch.pureNback] = mean_sem_sd(RT_avg.perSub.perEch.pureNback, 2);
[RT_avg_sd.aSubs.perEch.pureNback,...
    RT_sem_sd.aSubs.perEch.pureNback]  = mean_sem_sd(RT_sd.perSub.perEch.pureNback, 2);
% average and SD per effort level and answer
[RT_avg.aSubs.perEch_per_n_Nback,...
    RT_sem.aSubs.perEch_per_n_Nback]  = mean_sem_sd(RT_avg.perSub.perEch_per_n_Nback, 3);
[RT_avg_sd.aSubs.perEch_per_n_Nback,...
    RT_sem_sd.aSubs.perEch_per_n_Nback]  = mean_sem_sd(RT_sd.perSub.perEch_per_n_Nback, 3);
% average and SD per answer
[RT_avg.aSubs.per_n_Nback,...
    RT_sem.aSubs.per_n_Nback]  = mean_sem_sd(RT_avg.perSub.per_n_Nback, 2);
[RT_avg_sd.aSubs.per_n_Nback,...
    RT_sem_sd.aSubs.per_n_Nback]  = mean_sem_sd(RT_sd.perSub.per_n_Nback, 2);
n_max_resp = find(~isnan(RT_avg.aSubs.per_n_Nback),1,'last');
%% figures
lWidth = 3;
pSize = 40;

%% average RT and avg SD per effort level including 2 first answers
fig;

subplot(1,2,1);
hold on;
avgRT_hdl = errorbar(Ech_levels, RT_avg.aSubs.perEch.with2first,...
    RT_sem.aSubs.perEch.with2first);
avgRT_hdl.LineWidth = lWidth;
xticks(Ech_levels);
xlabel('Effort level');
ylabel('mean RT (s) (with first 2)');
legend_size(pSize);

subplot(1,2,2);
hold on;
sdRT_hdl = errorbar(Ech_levels, RT_avg_sd.aSubs.perEch.with2first,...
    RT_sem_sd.aSubs.perEch.with2first);
sdRT_hdl.LineWidth = lWidth;
xticks(Ech_levels);
xlabel('Effort level');
ylabel('sd RT (s) (with first 2)');
legend_size(pSize);

%% average RT and avg SD per effort level only for Nback replies
fig;

subplot(1,2,1);
hold on;
avgRT_hdl = errorbar(Ech_levels, RT_avg.aSubs.perEch.pureNback,...
    RT_sem.aSubs.perEch.pureNback);
avgRT_hdl.LineWidth = lWidth;
xticks(Ech_levels);
xlabel('Effort level');
ylabel('mean RT (s)');
legend_size(pSize);

subplot(1,2,2);
hold on;
sdRT_hdl = errorbar(Ech_levels, RT_avg_sd.aSubs.perEch.pureNback,...
    RT_sem_sd.aSubs.perEch.pureNback);
sdRT_hdl.LineWidth = lWidth;
xticks(Ech_levels);
xlabel('Effort level');
ylabel('sd RT (s)');
legend_size(pSize);

%% average and SD per effort level and answer
colours = {'k','m','b','r'};
fig;

subplot(1,2,1);
hold on;
for iEch = 1:n_E_ch
    avgRT_Nback_Ech_answer_hdl.(['Ech',num2str(iEch)]) = errorbar(1:n_max_resp,...
        RT_avg.aSubs.perEch_per_n_Nback(iEch,1:n_max_resp),...
        RT_sem.aSubs.perEch_per_n_Nback(iEch,1:n_max_resp));
    avgRT_Nback_Ech_answer_hdl.(['Ech',num2str(iEch)]).LineWidth = lWidth;
    avgRT_Nback_Ech_answer_hdl.(['Ech',num2str(iEch)]).Color = colours{iEch};
end
legend([avgRT_Nback_Ech_answer_hdl.Ech1,...
    avgRT_Nback_Ech_answer_hdl.Ech2,...
    avgRT_Nback_Ech_answer_hdl.Ech3,...
    avgRT_Nback_Ech_answer_hdl.Ech4],...
    {'Ech 0','Ech 1','Ech 2','Ech 3'});
legend('boxoff');
xlabel('Response');
ylabel('mean RT (s)');
legend_size(pSize);

subplot(1,2,2);
hold on;
for iEch = 1:n_E_ch
    sdRT_Nback_Ech_answer_hdl.(['Ech',num2str(iEch)]) = errorbar(1:n_max_resp,...
        RT_avg_sd.aSubs.perEch_per_n_Nback(iEch,1:n_max_resp),...
        RT_sem_sd.aSubs.perEch_per_n_Nback(iEch,1:n_max_resp));
    sdRT_Nback_Ech_answer_hdl.(['Ech',num2str(iEch)]).LineWidth = lWidth;
    sdRT_Nback_Ech_answer_hdl.(['Ech',num2str(iEch)]).Color = colours{iEch};
end
legend([sdRT_Nback_Ech_answer_hdl.Ech1,...
    sdRT_Nback_Ech_answer_hdl.Ech2,...
    sdRT_Nback_Ech_answer_hdl.Ech3,...
    sdRT_Nback_Ech_answer_hdl.Ech4],...
    {'Ech 0','Ech 1','Ech 2','Ech 3'});
xlabel('Response');
ylabel('sd RT (s)');
legend_size(pSize);

%% average and SD per answer
fig;

subplot(1,2,1);
hold on;
avgRT_hdl = errorbar(1:n_max_resp,...
    RT_avg.aSubs.per_n_Nback(1:n_max_resp),...
    RT_sem.aSubs.per_n_Nback(1:n_max_resp));
avgRT_hdl.LineWidth = lWidth;
xlabel('Response');
ylabel('mean RT (s)');
legend_size(pSize);

subplot(1,2,2);
hold on;
avgRT_hdl = errorbar(1:n_max_resp,...
    RT_avg_sd.aSubs.per_n_Nback(1:n_max_resp),...
    RT_sem_sd.aSubs.per_n_Nback(1:n_max_resp));
avgRT_hdl.LineWidth = lWidth;
xlabel('Response');
ylabel('sd RT (s)');
legend_size(pSize);
end % function