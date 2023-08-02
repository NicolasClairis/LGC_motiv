function[b, pval] = Em_efficiency_f_time(study_nm, condition, subject_id)
% [b, pval] = Em_efficiency_f_time(study_nm, condition, subject_id)
% Em_efficiency_f_time will look whether the efficiency in doing the mental
% effort task, defined as (nb correct answer/total time) changed over time.
%
% INPUTS
% study_nm: study name
%
% condition: condition indicating which subjects and runs to include
%
% subject_id: list of subjects
%
% OUTPUTS
% b: structure with betas
%
% pval: structure with p.values

%% subject selection
if ~exist('study_nm','var')
    study_nm = [];
end
if ~exist('condition','var')
    condition = [];
end
if ~exist('subject_id','var')
    subject_id = [];
end
[study_nm, condition, subject_id, NS] = subject_selection(study_nm, condition, subject_id);

%% working directory
computerRoot = LGCM_root_paths;
study_path = [computerRoot, study_nm, filesep];
%% loop through subjects
nTrialsPerRun = 54;
nRunsPerTask = 2;
nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
nBins = 6;
[efficiency.allTrials,...
    trialN.allTrials,...
    run1.allTrials,...
    run2.allTrials,...
    efficiency_fit.allTrials] = deal(NaN(nTrialsPerTask, NS));
[b.efficiency_f_time] = deal(NaN(1,NS));
[efficiency.trial_bins,...
    trialN.trial_bins,...
    efficiency_fit.trial_bins] = deal(NaN(nBins, NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [study_path, 'CID',sub_nm, filesep, 'behavior', filesep];
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    for iR = 1:n_runs
        kRun = runs.runsToKeep(iR);
        run_nm = num2str(kRun);
        task_id = runs.tasks{iR};
        if strcmp(task_id,'Em')
            switch kRun
                case {1,2}
                    jRun = 1;
                case {3,4}
                    jRun = 2;
            end
            run_trials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun - 1);
            [~, ~, ~, ~,...
                ~,...
                ~,...
                ~,...
                efficacy_pureNback_tmp,...
                ~,...
                efficacy_bis_pureNback_tmp] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
            
%             efficiency.allTrials(run_trials_idx, iS) = efficacy_pureNback_tmp.allTrials;
            efficiency.allTrials(run_trials_idx, iS) = efficacy_bis_pureNback_tmp.allTrials;
            trialN.allTrials(run_trials_idx, iS) = 1:nTrialsPerRun;
            switch kRun
                case 1
                    run1.allTrials(run_trials_idx, iS) = 1;
                    run2.allTrials(run_trials_idx, iS) = 0;
                case 2
                    run1.allTrials(run_trials_idx, iS) = 0;
                    run2.allTrials(run_trials_idx, iS) = 1;
            end
        end % run loop
        
        %% perform GLMs
        % general x regressor, including one constant for each run
        x_regs = [run1.allTrials(:, iS),...
            run2.allTrials(:, iS),...
            trialN.allTrials(:, iS)];
        
        % confidence rating = f(trial) with sigmoid (binary output)
        b_efficiency_tmp = glmfit(x_regs, efficiency.allTrials(:,iS),...
            'normal','Constant','off');
        b.efficiency_f_time(iS) = b_efficiency_tmp(3);
        efficiency_fit.allTrials(:,iS) = glmval(b_efficiency_tmp, x_regs,...
            'identity','constant','off');
        
        %% extract bins
        % actual data
        [efficiency.trial_bins(:,iS),...
            trialN.trial_bins(:,iS)] = do_bin2(efficiency.allTrials(:, iS), trialN.allTrials(:, iS), nBins, 0);
        % fitted data
        [efficiency_fit.trial_bins(:,iS)] = do_bin2(efficiency_fit.allTrials(:, iS), trialN.allTrials(:, iS), nBins, 0);
    end % task filter
end % subject loop

%% test at the group-level
[~,pval.efficiency_f_time] = ttest(b.efficiency_f_time);

%% mean and SEM across subjects
[m_efficiency.trial_bins,...
    sem_efficiency.trial_bins] = mean_sem_sd(efficiency.trial_bins, 2);
[m_trialN.trial_bins,...
    sem_trialN.trial_bins] = mean_sem_sd(trialN.trial_bins, 2);

[m_efficiency_fit.trial_bins,...
    sem_efficiency_fit.trial_bins] = mean_sem_sd(efficiency_fit.trial_bins, 2);
[b.mean.efficiency_f_time,...
    b.sem.efficiency_f_time] = mean_sem_sd(b.efficiency_f_time, 2);

%% figure
[pSize, lWidth, col, mSize] = general_fig_prm;

% efficiency = f(time)
fig;
er_hdl = errorbar(m_trialN.trial_bins,...
    m_efficiency.trial_bins,...
    sem_efficiency.trial_bins);
er_hdl.LineStyle = 'none';
er_hdl.Color = col.black;
er_hdl.Marker = 'o';
er_hdl.MarkerSize = mSize;
er_hdl.LineWidth = lWidth;
% add fit
fit_hdl = plot(m_trialN.trial_bins,...
    m_efficiency_fit.trial_bins);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = col.grey;
fit_hdl.LineStyle = '--';
xlabel('Trial');
ylabel('Efficiency');
legend_size(pSize);


end % function