function[b, pval] = conf_f_time(study_nm, condition, subject_id)
% [b, pval] = conf_f_time(study_nm, condition, subject_id)
% conf_f_time will look, for each task (physical/mental) whether
% subjective confidence estimated by the model and by actual ratings
% changes over time.
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

%% bayesian model to use
mdl_nm = 'mdl_3';

%% working directory
computerRoot = LGCM_root_paths;
study_path = [computerRoot, study_nm, filesep];
switch computerRoot
    case ['E:',filesep]
        gitFolder = fullfile('C:','Users','clairis','Desktop');
    case [filesep,filesep,fullfile('svfas5',...
            'sandi-lab','human_data_private','raw_data_subject'),filesep]
        gitFolder = fullfile('C:','Users','Nicolas Clairis','Documents');
end
resultsFolder = [gitFolder, filesep,...
    fullfile('GitHub','LGC_motiv','LGC_Motiv_results',...
    study_nm,'bayesian_modeling'),filesep];
%% loop through subjects
nTrialsPerRun = 54;
nRunsPerTask = 2;
nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
tasks = {'Ep','Em'};
nTasks = length(tasks);
nBins = 6;
for iT = 1:nTasks
    task_nm = tasks{iT};
    [confRtg.(task_nm).allTrials,...
        confInferred.(task_nm).allTrials,...
        trialN.(task_nm).allTrials,...
        run1.(task_nm).allTrials,...
        run2.(task_nm).allTrials,...
        confRtg_fit.(task_nm).allTrials,...
        confInferred_fit.(task_nm).allTrials] = deal(NaN(nTrialsPerTask, NS));
    [b.confRtg_f_time.(task_nm),...
        b.confInferred_f_time.(task_nm)] = deal(NaN(1,NS));
    [confRtg.(task_nm).trial_bins,...
        confInferred.(task_nm).trial_bins,...
        trialN.(task_nm).trial_bins,...
        confRtg_fit.(task_nm).trial_bins,...
        confInferred_fit.(task_nm).trial_bins] = deal(NaN(nBins, NS));
end
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [study_path, 'CID',sub_nm, filesep, 'behavior', filesep];
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    for iR = 1:n_runs
        kRun = runs.runsToKeep(iR);
        run_nm = num2str(kRun);
        task_id = runs.tasks{iR};
        switch task_id
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end
        switch kRun
            case {1,2}
                jRun = 1;
            case {3,4}
                jRun = 2;
        end
        run_trials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun - 1);
        confRtg_tmp = extract_confidence_rating(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [~, ~, confidence_highE_tmp, ~] = extract_bayesian_mdl(resultsFolder, subBehaviorFolder,...
            sub_nm, run_nm, task_fullName, mdl_nm);
        
        confRtg.(task_id).allTrials(run_trials_idx, iS) = confRtg_tmp;
        confInferred.(task_id).allTrials(run_trials_idx, iS) = confidence_highE_tmp;
        trialN.(task_id).allTrials(run_trials_idx, iS) = 1:nTrialsPerRun;
        switch kRun
            case 1
                run1.(task_id).allTrials(run_trials_idx, iS) = 1;
                run2.(task_id).allTrials(run_trials_idx, iS) = 0;
            case 2
                run1.(task_id).allTrials(run_trials_idx, iS) = 0;
                run2.(task_id).allTrials(run_trials_idx, iS) = 1;
        end
    end % run loop
    
    for iT = 1:nTasks
        task_nm = tasks{iT};
        
        %% perform GLMs
        % general x regressor, including one constant for each run
        x_regs = [run1.(task_nm).allTrials(:, iS),...
            run2.(task_nm).allTrials(:, iS),...
            trialN.(task_nm).allTrials(:, iS)];
        
        % confidence rating = f(trial) with sigmoid (binary output)
        b_confRtg_tmp = glmfit(x_regs, confRtg.(task_nm).allTrials(:,iS),...
            'binomial','link','logit','Constant','off');
        b.confRtg_f_time.(task_nm)(iS) = b_confRtg_tmp(3);
        confRtg_fit.(task_nm).allTrials(:,iS) = glmval(b_confRtg_tmp, x_regs,...
            'logit','constant','off');
        
        % confidence inferred = f(trial) (linear)
        b_confInferred_tmp = glmfit(x_regs, confInferred.(task_nm).allTrials(:,iS),...
            'normal','Constant','off');
        b.confInferred_f_time.(task_nm)(iS) = b_confInferred_tmp(3);
        confInferred_fit.(task_nm).allTrials(:,iS) = glmval(b_confInferred_tmp, x_regs,...
            'identity','Constant','off');
        
        %% extract bins
        % actual data
        [confRtg.(task_nm).trial_bins(:,iS),...
            trialN.(task_nm).trial_bins(:,iS)] = do_bin2(confRtg.(task_nm).allTrials(:, iS), trialN.(task_nm).allTrials(:, iS), nBins, 0);
        [confInferred.(task_nm).trial_bins(:,iS)] = do_bin2(confInferred.(task_nm).allTrials(:, iS), trialN.(task_nm).allTrials(:, iS), nBins, 0);
        % fitted data
        [confRtg_fit.(task_nm).trial_bins(:,iS)] = do_bin2(confRtg_fit.(task_nm).allTrials(:, iS), trialN.(task_nm).allTrials(:, iS), nBins, 0);
        [confInferred_fit.(task_nm).trial_bins(:,iS)] = do_bin2(confInferred_fit.(task_nm).allTrials(:, iS), trialN.(task_nm).allTrials(:, iS), nBins, 0);
    end % task loop
end % subject loop

for iT = 1:nTasks
    task_nm = tasks{iT};
    %% test at the group-level
    [~,pval.confRtg_f_time.(task_nm)] = ttest(b.confRtg_f_time.(task_nm));
    [~,pval.confInferred_f_time.(task_nm)] = ttest(b.confInferred_f_time.(task_nm));
    
    %% mean and SEM across subjects
    [m_confRtg.(task_nm).trial_bins,...
        sem_confRtg.(task_nm).trial_bins] = mean_sem_sd(confRtg.(task_nm).trial_bins, 2);
    [m_confInferred.(task_nm).trial_bins,...
        sem_confInferred.(task_nm).trial_bins] = mean_sem_sd(confInferred.(task_nm).trial_bins, 2);
    [m_trialN.(task_nm).trial_bins,...
        sem_trialN.(task_nm).trial_bins] = mean_sem_sd(trialN.(task_nm).trial_bins, 2);
    
    [m_confRtg_fit.(task_nm).trial_bins,...
        sem_confRtg_fit.(task_nm).trial_bins] = mean_sem_sd(confRtg_fit.(task_nm).trial_bins, 2);
    [m_confInferred_fit.(task_nm).trial_bins,...
        sem_confInferred_fit.(task_nm).trial_bins] = mean_sem_sd(confInferred_fit.(task_nm).trial_bins, 2);
    [b.confRtg_f_time.mean.(task_nm),...
        b.confRtg_f_time.sem.(task_nm)] = mean_sem_sd(b.confRtg_f_time.(task_nm), 2);
    [b.confInferred_f_time.mean.(task_nm),...
        b.confInferred_f_time.sem.(task_nm)] = mean_sem_sd(b.confInferred_f_time.(task_nm), 2);
end % task loop

%% figure
[pSize, lWidth, col, mSize] = general_fig_prm;
confRtgFig = fig;
confInfFig = fig;

for iT = 1:nTasks
    task_nm = tasks{iT};
    
    % conf rating = f(time)
    figure(confRtgFig);
    subplot(1,nTasks, iT); hold on;
    er_hdl = errorbar(m_trialN.(task_nm).trial_bins,...
        m_confRtg.(task_nm).trial_bins,...
        sem_confRtg.(task_nm).trial_bins);
    er_hdl.LineStyle = 'none';
    er_hdl.Color = col.black;
    er_hdl.Marker = 'o';
    er_hdl.MarkerSize = mSize;
    er_hdl.LineWidth = lWidth;
    % add fit
    fit_hdl = plot(m_trialN.(task_nm).trial_bins,...
        m_confRtg_fit.(task_nm).trial_bins);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = col.grey;
    fit_hdl.LineStyle = '--';
    xlabel('Trial');
    ylabel('Confidence rating');
    legend_size(pSize);
    
    % conf inferred = f(time)
    figure(confInfFig);
    subplot(1,nTasks, iT); hold on;
    er_hdl = errorbar(m_trialN.(task_nm).trial_bins,...
        m_confInferred.(task_nm).trial_bins,...
        sem_confInferred.(task_nm).trial_bins);
    er_hdl.LineStyle = 'none';
    er_hdl.Color = 'k';
    er_hdl.Marker = 'o';
    er_hdl.MarkerSize = mSize;
    er_hdl.LineWidth = lWidth;
    % add fit
    fit_hdl = plot(m_trialN.(task_nm).trial_bins,...
        m_confInferred_fit.(task_nm).trial_bins);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = col.grey;
    fit_hdl.LineStyle = '--';
    xlabel('Trial');
    ylabel('Confidence');
    legend_size(pSize);
end % task loop


end % function