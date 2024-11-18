% function[] = choice_Em_f_Em_efficacy()
% choice_Em_f_Em_efficacy aims at assessing whether the choices in the mental
% effort task depend on the efficacy in terms of correct answers and
% time spent solving the N-back task in previous trials.
%
% INPUTS
% 
% OUTPUTS
%

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
n_bins = 9;
[choice_f_efficacy_with2first, efficacy_f_efficacy_with2first,...
    choice_f_efficacy_pureNback, efficacy_f_efficacy_pureNback,...
    choice_f_trialN,...
    nCorrect_f_trialN, nErrors_f_trialN,...
    RT_avg_f_trialN_with2first, RT_avg_f_trialN_pureNback,...
    eff_with2first_f_trialN, eff_pureNback_f_trialN,...
    trialN_f_trialN] = deal(NaN(n_bins, NS));
[r_corr.perSub.choice_f_efficacy_with2first,...
    r_corr.perSub.choice_f_efficacy_pureNback,...
    pval.perSub.choice_f_efficacy_with2first,...
    pval.perSub.choice_f_efficacy_pureNback] = deal(NaN(1, NS));

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
    [RT_avg_run_with2first_tmp,...
        RT_avg_run_pureNback_tmp,...
        nCorrect_avg_run_tmp,...
        nErrors_avg_run_tmp,...
        RTtotal_avg_run_with2first_tmp,...
        RTtotal_avg_run_pureNback_tmp,...
        efficacy_avg_run_with2first_tmp,...
        efficacy_avg_run_pureNback_tmp,...
        efficacy_avg_run_with2first_tmp,...
        efficacy_avg_run_pureNback_tmp,...
        choice_hE_run_tmp,...
        trialN_run_tmp] = deal(NaN(nTotalEmTrials,1));

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
            % choices
            default_LR = mentalE_perf.choiceOptions.default_LR;
            choice_LR = mentalE_perf.choice;
            choice_LR(choice_LR == 2) = 1;
            choice_LR(choice_LR == -2) = -1;
            choice_hE = NaN(1,nTrialsPerRun);
            choice_hE(choice_LR == default_LR) = 0;
            choice_hE(choice_LR == -default_LR) = 1;

            for iTrial = 1:nTrialsPerRun
                jTrial = iTrial + nTrialsPerRun*(jRun - 1);
                rt_tmp = mentalE_perf.perfSummary{1,iTrial}.rt;
                % average RT including 2 first answers (which are just
                % for initialization)
                RT_avg_run_with2first_tmp(jTrial) = mean(rt_tmp, 2, 'omitnan');
                % average RT for N-back answers only
                RT_avg_run_pureNback_tmp(jTrial) = mean(rt_tmp(3:end), 2, 'omitnan');
                % choice for current trial
                choice_hE_run_tmp(jTrial) = choice_hE(iTrial);
                % trial number
                trialN_run_tmp(jTrial) = iTrial;
                % number of errors
                nErrors_avg_run_tmp(jTrial) =  mentalE_perf.perfSummary{1,iTrial}.n_errorsMade;
                % number of corrects
                nCorrect_avg_run_tmp(jTrial) =  mentalE_perf.perfSummary{1,iTrial}.n_correctAnswersProvided;
                % extract total RT
                RTtotal_avg_run_with2first_tmp(jTrial) = mentalE_perf.perfSummary{1,iTrial}.totalTime_success;
                onset_names_tmp = fieldnames(mentalE_perf.perfSummary{1,iTrial}.onsets);
                if ~isempty(onset_names_tmp) && sum(strcmp(onset_names_tmp,'nb_3')) == 1
                    trialStart_tmp = mentalE_perf.perfSummary{1,iTrial}.onsets.nb_3;
                    trialEnd_tmp = mentalE_perf.perfSummary{1,iTrial}.onsets.(['nb_',num2str(length(onset_names_tmp))]) +...
                        rt_tmp(end);
                    RTtotal_avg_run_pureNback_tmp(jTrial) = trialEnd_tmp - trialStart_tmp;
                end
                % efficacy
                efficacy_avg_run_with2first_tmp(jTrial) = (nCorrect_avg_run_tmp(jTrial) - nErrors_avg_run_tmp(jTrial))./RTtotal_avg_run_with2first_tmp(jTrial);
                efficacy_avg_run_pureNback_tmp(jTrial) = (nCorrect_avg_run_tmp(jTrial) - nErrors_avg_run_tmp(jTrial))./RTtotal_avg_run_pureNback_tmp(jTrial);
            end % trial loop
        end % run to include?
    end % run loop
    
    %% perform correlation test
    okTrials_2first_tmp = ~isnan(choice_hE_run_tmp.*efficacy_avg_run_with2first_tmp);
    okTrials_Nback_tmp = ~isnan(choice_hE_run_tmp.*efficacy_avg_run_pureNback_tmp);
    [r_corr.perSub.choice_f_efficacy_with2first(iS),...
        pval.perSub.choice_f_efficacy_with2first(iS)] = corr(choice_hE_run_tmp(okTrials_2first_tmp), efficacy_avg_run_with2first_tmp(okTrials_2first_tmp));
    [r_corr.perSub.choice_f_efficacy_pureNback(iS),...
        pval.perSub.choice_f_efficacy_pureNback(iS)] = corr(choice_hE_run_tmp(okTrials_Nback_tmp), efficacy_avg_run_pureNback_tmp(okTrials_Nback_tmp));
    
    %% average data per subject
    [choice_f_efficacy_with2first(:,iS), efficacy_f_efficacy_with2first(:,iS)] = do_bin2(choice_hE_run_tmp, efficacy_avg_run_with2first_tmp, n_bins,0);
    [choice_f_efficacy_pureNback(:,iS), efficacy_f_efficacy_pureNback(:,iS)] = do_bin2(choice_hE_run_tmp, efficacy_avg_run_pureNback_tmp, n_bins,0);
    [choice_f_trialN(:,iS), trialN_f_trialN(:,iS)] = do_bin2(choice_hE_run_tmp, trialN_run_tmp, n_bins,0);
    [nCorrect_f_trialN(:,iS)] = do_bin2(nCorrect_avg_run_tmp, trialN_run_tmp, n_bins,0);
    [nErrors_f_trialN(:,iS)] = do_bin2(nErrors_avg_run_tmp, trialN_run_tmp, n_bins,0);
    [RT_avg_f_trialN_with2first(:,iS)] = do_bin2(RT_avg_run_with2first_tmp, trialN_run_tmp, n_bins,0);
    [RT_avg_f_trialN_pureNback(:,iS)] = do_bin2(RT_avg_run_pureNback_tmp, trialN_run_tmp, n_bins,0);
    [efficacy_with2first_f_trialN(:,iS)] = do_bin2(efficacy_avg_run_with2first_tmp, trialN_run_tmp, n_bins,0);
    [efficacy_pureNback_f_trialN(:,iS)] = do_bin2(efficacy_avg_run_pureNback_tmp, trialN_run_tmp, n_bins,0);
end % subject loop

%% average across subjects
% correlation coefficients
[r_corr.aSubs.choice_f_efficacy_with2first.mean,...
    r_corr.aSubs.choice_f_efficacy_with2first.sem] = mean_sem_sd(r_corr.perSub.choice_f_efficacy_with2first,2);
[r_corr.aSubs.choice_f_efficacy_pureNback.mean,...
    r_corr.aSubs.choice_f_efficacy_pureNback.sem] = mean_sem_sd(r_corr.perSub.choice_f_efficacy_pureNback,2);
% test if correlation coefficients are significantly different from zero
[~,pval.aSubs.choice_f_efficacy_with2first] = ttest(r_corr.perSub.choice_f_efficacy_with2first);
[~,pval.aSubs.choice_f_efficacy_pureNback] = ttest(r_corr.perSub.choice_f_efficacy_pureNback);
% f efficacy
[choice_f_efficacy_with2first_avg,...
    choice_f_efficacy_with2first_sem] = mean_sem_sd(choice_f_efficacy_with2first,2);
[efficacy_f_efficacy_with2first_avg,...
    efficacy_f_efficacy_with2first_sem] = mean_sem_sd(efficacy_f_efficacy_with2first,2);
[choice_f_efficacy_pureNback_avg,...
    choice_f_efficacy_pureNback_sem] = mean_sem_sd(choice_f_efficacy_pureNback,2);
[efficacy_f_efficacy_pureNback_avg,...
    efficacy_f_efficacy_pureNback_sem] = mean_sem_sd(efficacy_f_efficacy_pureNback,2);
% f trial number
[choice_f_trialN_avg,...
    choice_f_trialN_sem] = mean_sem_sd(choice_f_trialN,2);
[trialN_f_trialN_avg,...
    trialN_f_trialN_sem] = mean_sem_sd(trialN_f_trialN,2);
[nCorrect_f_trialN_avg,...
    nCorrect_f_trialN_sem] = mean_sem_sd(nCorrect_f_trialN,2);
[nErrors_f_trialN_avg,...
    nErrors_f_trialN_sem] = mean_sem_sd(nErrors_f_trialN,2);
[RT_avg_f_trialN_with2first_avg,...
    RT_avg_f_trialN_with2first_sem] = mean_sem_sd(RT_avg_f_trialN_with2first,2);
[RT_avg_f_trialN_pureNback_avg,...
    RT_avg_f_trialN_pureNback_sem] = mean_sem_sd(RT_avg_f_trialN_pureNback,2);
[efficacy_with2first_f_trialN_avg,...
    efficacy_with2first_f_trialN_sem] = mean_sem_sd(efficacy_with2first_f_trialN,2);
[efficacy_pureNback_f_trialN_avg,...
    efficacy_pureNback_f_trialN_sem] = mean_sem_sd(efficacy_pureNback_f_trialN,2);

%% figures
lWidth = 3;
pSize = 40;

%% choice = f(efficacy) with 2 first included
fig;
hold on;
hdl = errorbar(efficacy_f_efficacy_with2first_avg,...
    choice_f_efficacy_with2first_avg,...
    choice_f_efficacy_with2first_sem);
hdl.LineWidth = lWidth;
xlabel('Efficacy (with 2 first included)');
ylabel('Choice hE (%)');
legend_size(pSize);
place_r_and_pval(r_corr.aSubs.choice_f_efficacy_with2first.mean,...
    pval.aSubs.choice_f_efficacy_with2first);

%% choice = f(efficacy) only with relevant times
fig;
hold on;
hdl = errorbar(efficacy_f_efficacy_pureNback_avg,...
    choice_f_efficacy_pureNback_avg,...
    choice_f_efficacy_pureNback_sem);
hdl.LineWidth = lWidth;
xlabel('Efficacy');
ylabel('Choice hE (%)');
legend_size(pSize);
place_r_and_pval(r_corr.aSubs.choice_f_efficacy_pureNback.mean,...
    pval.aSubs.choice_f_efficacy_pureNback);

%% performance = f(trial number) with 2 first included
fig;
hold on;
hdl = errorbar(trialN_f_trialN_avg,...
    efficacy_with2first_f_trialN_avg,...
    efficacy_with2first_f_trialN_sem);
hdl.LineWidth = lWidth;
xticks(trialN_f_trialN_avg);
xLabelNames = cell(1,n_bins);
for iTick = 1:n_bins
    xLabelNames{iTick} = num2str(round(trialN_f_trialN_avg(iTick),0));
end
xticklabels(xLabelNames);
xlabel('Trial number');
ylabel('Efficacy (with 2 first included)');
legend_size(pSize);

%% performance = f(trial number) only with relevant times
fig;
hold on;
hdl = errorbar(trialN_f_trialN_avg,...
    efficacy_pureNback_f_trialN_avg,...
    efficacy_pureNback_f_trialN_sem);
hdl.LineWidth = lWidth;
xticks(trialN_f_trialN_avg);
xLabelNames = cell(1,n_bins);
for iTick = 1:n_bins
    xLabelNames{iTick} = num2str(round(trialN_f_trialN_avg(iTick),0));
end
xticklabels(xLabelNames);
xlabel('Trial number');
ylabel('Efficacy');
legend_size(pSize);

% end % function