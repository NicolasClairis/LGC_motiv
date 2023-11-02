% script to test correlation between saliency (|dSV+bias|), subjective
% value (dSV+bias), inferred confidence ([p(chosen)-0.5]²) and rated
% confidence.

%% subject identification
[study_nm, condition, gender, subject_id, NS] = sub_id;

%% working dir
computerRoot = LGCM_root_paths;
studyPath = [computerRoot, study_nm, filesep];
gitResultsFolder = [fullfile('C:','Users','clairis','Desktop','GitHub',...
    'LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];

%% model to use
bayesianMdl_nm = 'mdl_3';

%% initialize vars of interest
nTrialsPerRun = 54;
nRuns = 4;
nTrials = nTrialsPerRun.*nRuns;
[conf_rtg, conf_inferred,...
    saliency, SV,...
    run_cstt.run1, run_cstt.run2,...
    run_cstt.run3, run_cstt.run4] = deal(NaN(nTrials, NS));
[beta_SV_f_saliency, beta_SV_f_conf_rtg,...
    beta_SV_f_conf_inferred, beta_saliency_f_conf_rtg,...
    beta_saliency_f_conf_inferred, beta_conf_rtg_f_conf_inferred] = deal(NaN(nRuns+1,NS));
[r_correl.SV_f_saliency, r_correl.SV_f_conf_rtg,...
    r_correl.SV_f_conf_inferred, r_correl.saliency_f_conf_rtg,...
    r_correl.saliency_f_conf_inferred, r_correl.conf_rtg_f_conf_inferred] = deal(NaN(nRuns+1,NS));
nBins = 6;
[SV_f_saliency, saliency_f_saliency, SV_f_saliency_fit,...
    SV_f_conf_rtg, conf_rtg_f_conf_rtg, SV_f_conf_rtg_fit,...
    SV_f_conf_inferred, conf_inferred_f_conf_inferred, SV_f_conf_inferred_fit,...
    saliency_f_conf_rtg, saliency_f_conf_rtg_fit,...
    saliency_f_conf_inferred, saliency_f_conf_inferred_fit,...
    conf_rtg_f_conf_inferred, conf_rtg_f_conf_inferred_fit] = deal(NaN(nBins,NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subPath = [studyPath, 'CID',sub_nm, filesep, 'behavior', filesep];
    
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    for iR = 1:n_runs
        jR = runs.runsToKeep(iR);
        run_nm = num2str(jR);
        switch runs.tasks{iR}
            case 'Ep'
                task_fullName = 'physical';
                task_result_field = 'physicalPerf';
            case 'Em'
                task_fullName = 'mental';
                task_result_field = 'mentalE_perf';
        end
        
        %% extract relevant data
        % load rated confidence
        [conf_rtg_tmp] = extract_confidence_rating(subPath, sub_nm, run_nm, task_fullName);
        % load saliency, SV and inferred confidence
        [~, ~,...
            conf_inferred_tmp, ~,...
            NV_highE_tmp, NVch_min_NVunch_tmp] = extract_bayesian_mdl(gitResultsFolder, subPath,...
            sub_nm, run_nm, task_fullName, bayesianMdl_nm);
        % index of the trials to consider
        trials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jR - 1);
        conf_rtg(trials_idx,iS) = conf_rtg_tmp;
        conf_inferred(trials_idx,iS) = conf_inferred_tmp;
        saliency(trials_idx,iS) = abs(NV_highE_tmp);
        SV(trials_idx,iS) = NVch_min_NVunch_tmp;
        switch jR
            case 1
                run_cstt.run1(trials_idx,iS) = 1;
                run_cstt.run2(trials_idx,iS) = 0;
                run_cstt.run3(trials_idx,iS) = 0;
                run_cstt.run1(trials_idx,iS) = 0;
            case 2
                run_cstt.run1(trials_idx,iS) = 0;
                run_cstt.run2(trials_idx,iS) = 1;
                run_cstt.run3(trials_idx,iS) = 0;
                run_cstt.run4(trials_idx,iS) = 0;
            case 3
                run_cstt.run1(trials_idx,iS) = 0;
                run_cstt.run2(trials_idx,iS) = 0;
                run_cstt.run3(trials_idx,iS) = 1;
                run_cstt.run4(trials_idx,iS) = 0;
            case 4
                run_cstt.run1(trials_idx,iS) = 0;
                run_cstt.run2(trials_idx,iS) = 0;
                run_cstt.run3(trials_idx,iS) = 0;
                run_cstt.run4(trials_idx,iS) = 1;
        end
    end % run loop
    
    %% perform GLM
    run_cstt_tmp.run1 = run_cstt.run1(:,iS);
    run_cstt_tmp.run2 = run_cstt.run2(:,iS);
    run_cstt_tmp.run3 = run_cstt.run3(:,iS);
    run_cstt_tmp.run4 = run_cstt.run4(:,iS);
    [beta_SV_f_saliency(:,iS), SV_f_saliency_fit_tmp] = glmfit_across_multiple_runs(run_cstt_tmp, saliency(:,iS), SV(:,iS), nRuns);
    [beta_SV_f_conf_rtg(:,iS), SV_f_conf_rtg_fit_tmp] = glmfit_across_multiple_runs(run_cstt_tmp, conf_rtg(:,iS), SV(:,iS), nRuns);
    [beta_SV_f_conf_inferred(:,iS), SV_f_conf_inferred_fit_tmp] = glmfit_across_multiple_runs(run_cstt_tmp, conf_inferred(:,iS), SV(:,iS), nRuns);
    [beta_saliency_f_conf_rtg(:,iS), saliency_f_conf_rtg_fit_tmp] = glmfit_across_multiple_runs(run_cstt_tmp, conf_rtg(:,iS), saliency(:,iS), nRuns);
    [beta_saliency_f_conf_inferred(:,iS), saliency_f_conf_inferred_fit_tmp] = glmfit_across_multiple_runs(run_cstt_tmp, conf_inferred(:,iS), saliency(:,iS), nRuns);
    [beta_conf_rtg_f_conf_inferred(:,iS), conf_rtg_f_conf_inferred_fit_tmp] = glmfit_across_multiple_runs(run_cstt_tmp, conf_inferred(:,iS), conf_rtg(:,iS), nRuns);
    % same with zscored data to obtain R
    [r_correl.SV_f_saliency(:,iS)] = glmfit_across_multiple_runs(run_cstt_tmp, nanzscore(saliency(:,iS)), nanzscore(SV(:,iS)), nRuns);
    [r_correl.SV_f_conf_rtg(:,iS)] = glmfit_across_multiple_runs(run_cstt_tmp, nanzscore(conf_rtg(:,iS)), nanzscore(SV(:,iS)), nRuns);
    [r_correl.SV_f_conf_inferred(:,iS)] = glmfit_across_multiple_runs(run_cstt_tmp, nanzscore(conf_inferred(:,iS)), nanzscore(SV(:,iS)), nRuns);
    [r_correl.saliency_f_conf_rtg(:,iS),] = glmfit_across_multiple_runs(run_cstt_tmp, nanzscore(conf_rtg(:,iS)), nanzscore(saliency(:,iS)), nRuns);
    [r_correl.saliency_f_conf_inferred(:,iS)] = glmfit_across_multiple_runs(run_cstt_tmp, nanzscore(conf_inferred(:,iS)), nanzscore(saliency(:,iS)), nRuns);
    [r_correl.conf_rtg_f_conf_inferred(:,iS)] = glmfit_across_multiple_runs(run_cstt_tmp, nanzscore(conf_inferred(:,iS)), nanzscore(conf_rtg(:,iS)), nRuns);
    
    %% perform split with bins
    [SV_f_saliency(:,iS), saliency_f_saliency(:,iS)] = do_bin2(SV(:,iS), saliency(:,iS), nBins, 0);
    [SV_f_saliency_fit(:,iS)] = do_bin2(SV_f_saliency_fit_tmp, saliency(:,iS), nBins, 0);
    
    [SV_f_conf_rtg(:,iS), conf_rtg_f_conf_rtg(:,iS)] = do_bin2(SV(:,iS), conf_rtg(:,iS), nBins, 0);
    [SV_f_conf_rtg_fit(:,iS)] = do_bin2(SV_f_conf_rtg_fit_tmp, conf_rtg(:,iS), nBins, 0);
    
    [SV_f_conf_inferred(:,iS), conf_inferred_f_conf_inferred(:,iS)] = do_bin2(SV(:,iS), conf_inferred(:,iS), nBins, 0);
    [SV_f_conf_inferred_fit(:,iS)] = do_bin2(SV_f_conf_inferred_fit_tmp, conf_inferred(:,iS), nBins, 0);
    
    saliency_f_conf_rtg(:,iS) = do_bin2(saliency(:,iS), conf_rtg(:,iS), nBins, 0);
    saliency_f_conf_rtg_fit(:,iS) = do_bin2(saliency_f_conf_rtg_fit_tmp, conf_rtg(:,iS), nBins, 0);
    
    saliency_f_conf_inferred(:,iS) = do_bin2(saliency(:,iS), conf_inferred(:,iS), nBins, 0);
    saliency_f_conf_inferred_fit(:,iS) = do_bin2(saliency_f_conf_inferred_fit_tmp, conf_inferred(:,iS), nBins, 0);
    
    conf_rtg_f_conf_inferred(:,iS) = do_bin2(conf_rtg(:,iS), conf_inferred(:,iS), nBins, 0);
    conf_rtg_f_conf_inferred_fit(:,iS) = do_bin2(conf_rtg_f_conf_inferred_fit_tmp, conf_inferred(:,iS), nBins, 0);
end % subject loop

%% average across subjects
% SV = f(saliency)
[m_SV_f_saliency, sem_SV_f_saliency] = mean_sem_sd(SV_f_saliency,2);
[m_saliency_f_saliency, sem_saliency_f_saliency] = mean_sem_sd(saliency_f_saliency,2);
[m_SV_f_saliency_fit, sem_SV_f_saliency_fit] = mean_sem_sd(SV_f_saliency_fit,2);
[~,pval.SV_f_saliency] = ttest(beta_SV_f_saliency(5,:));
r_correl.mean.SV_f_saliency = mean(r_correl.SV_f_saliency(5,:),2,'omitnan');

% SV = f(confidence rating)
[m_SV_f_conf_rtg, sem_SV_f_conf_rtg] = mean_sem_sd(SV_f_conf_rtg,2);
[m_conf_rtg_f_conf_rtg, sem_conf_rtg_f_conf_rtg] = mean_sem_sd(conf_rtg_f_conf_rtg,2);
[m_SV_f_conf_rtg_fit, sem_SV_f_conf_rtg_fit] = mean_sem_sd(SV_f_conf_rtg_fit,2);
[~,pval.SV_f_conf_rtg] = ttest(beta_SV_f_conf_rtg(5,:));
r_correl.mean.SV_f_conf_rtg = mean(r_correl.SV_f_conf_rtg(5,:),2,'omitnan');

% SV = f(confidence inferred)
[m_SV_f_conf_inferred, sem_SV_f_conf_inferred] = mean_sem_sd(SV_f_conf_inferred,2);
[m_conf_inferred_f_conf_inferred, sem_conf_inferred_f_conf_inferred] = mean_sem_sd(conf_inferred_f_conf_inferred,2);
[m_SV_f_conf_inferred_fit, sem_SV_f_conf_inferred_fit] = mean_sem_sd(SV_f_conf_inferred_fit,2);
[~,pval.SV_f_conf_inferred] = ttest(beta_SV_f_conf_inferred(5,:));
r_correl.mean.SV_f_conf_inferred = mean(r_correl.SV_f_conf_inferred(5,:),2,'omitnan');

% saliency = f(confidence rating)
[m_saliency_f_conf_rtg, sem_saliency_f_conf_rtg] = mean_sem_sd(saliency_f_conf_rtg,2);
[m_saliency_f_conf_rtg_fit, sem_saliency_f_conf_rtg_fit] = mean_sem_sd(saliency_f_conf_rtg_fit,2);
[~,pval.saliency_f_conf_rtg] = ttest(beta_saliency_f_conf_rtg(5,:));
r_correl.mean.saliency_f_conf_rtg = mean(r_correl.saliency_f_conf_rtg(5,:),2,'omitnan');

% saliency = f(confidence inferred)
[m_saliency_f_conf_inferred, sem_saliency_f_conf_inferred] = mean_sem_sd(saliency_f_conf_inferred,2);
[m_saliency_f_conf_inferred_fit, sem_saliency_f_conf_inferred_fit] = mean_sem_sd(saliency_f_conf_inferred_fit,2);
[~,pval.saliency_f_conf_inferred] = ttest(beta_saliency_f_conf_inferred(5,:));
r_correl.mean.saliency_f_conf_inferred = mean(r_correl.saliency_f_conf_inferred(5,:),2,'omitnan');

% confidence rating = f(confidence inferred)
[m_conf_rtg_f_conf_inferred, sem_conf_rtg_f_conf_inferred] = mean_sem_sd(conf_rtg_f_conf_inferred,2);
[m_conf_rtg_f_conf_inferred_fit, sem_conf_rtg_f_conf_inferred_fit] = mean_sem_sd(conf_rtg_f_conf_inferred_fit,2);
[~,pval.conf_rtg_f_conf_inferred] = ttest(beta_conf_rtg_f_conf_inferred(5,:));
r_correl.mean.conf_rtg_f_conf_inferred = mean(r_correl.conf_rtg_f_conf_inferred(5,:),2,'omitnan');

%% figures with global results
[pSize, lWidth, col, mSize] = general_fig_prm;

%% main figures
fig;

% SV = f(saliency)
subplot(1,3,1); hold on;
SV_f_saliency_hdl = errorbar(m_saliency_f_saliency,...
    m_SV_f_saliency, sem_SV_f_saliency);
errorbar_hdl_upgrade(SV_f_saliency_hdl);
SV_f_saliency_fit_hdl = plot(m_saliency_f_saliency, m_SV_f_saliency_fit);
fit_hdl_upgrade(SV_f_saliency_fit_hdl);
xlabel('Saliency');
ylabel('SV chosen - unchosen');
legend_size(pSize);

% SV = f(confidence inferred)
subplot(1,3,2); hold on;
SV_f_conf_inferred_hdl = errorbar(m_conf_inferred_f_conf_inferred,...
    m_SV_f_conf_inferred, sem_SV_f_conf_inferred);
errorbar_hdl_upgrade(SV_f_conf_inferred_hdl);
SV_f_inferred_fit_hdl = plot(m_conf_inferred_f_conf_inferred, m_SV_f_conf_inferred_fit);
fit_hdl_upgrade(SV_f_inferred_fit_hdl);
xlabel('Confidence');
ylabel('SV chosen - unchosen');
legend_size(pSize);

% saliency = f(confidence inferred)
subplot(1,3,3); hold on;
saliency_f_conf_hdl = errorbar(m_conf_inferred_f_conf_inferred,...
    m_saliency_f_conf_inferred, sem_saliency_f_conf_inferred);
errorbar_hdl_upgrade(saliency_f_conf_hdl);
saliency_f_conf_fit_hdl = plot(m_conf_inferred_f_conf_inferred, m_saliency_f_conf_inferred_fit);
fit_hdl_upgrade(saliency_f_conf_fit_hdl);
xlabel('Confidence');
ylabel('Saliency');
legend_size(pSize);

%% confidence rating figures
fig;

% saliency = f(confidence rating)
subplot(1,3,1); hold on;
saliency_f_conf_rtg_hdl = errorbar(m_conf_rtg_f_conf_rtg,...
    m_saliency_f_conf_rtg, sem_saliency_f_conf_rtg);
errorbar_hdl_upgrade(saliency_f_conf_rtg_hdl);
saliency_f_conf_rtg_fit_hdl = plot(m_conf_rtg_f_conf_rtg, m_saliency_f_conf_rtg_fit);
fit_hdl_upgrade(saliency_f_conf_rtg_fit_hdl);
xlabel('Confidence rating');
ylabel('Saliency');
legend_size(pSize);

% SV = f(confidence rated)
subplot(1,3,2); hold on;
SV_f_conf_rtg_hdl = errorbar(m_conf_rtg_f_conf_rtg,...
    m_SV_f_conf_rtg, sem_SV_f_conf_rtg);
errorbar_hdl_upgrade(SV_f_conf_rtg_hdl);
SV_f_conf_rtg_fit_hdl = plot(m_conf_rtg_f_conf_rtg, m_SV_f_conf_rtg_fit);
fit_hdl_upgrade(SV_f_conf_rtg_fit_hdl);
xlabel('Confidence rating');
ylabel('SV chosen - unchosen');
legend_size(pSize);

% confidence rating = f(confidence inferred)
subplot(1,3,3); hold on;
conf_rtg_f_conf_hdl = errorbar(m_conf_inferred_f_conf_inferred,...
    m_conf_rtg_f_conf_inferred, sem_conf_rtg_f_conf_inferred);
errorbar_hdl_upgrade(conf_rtg_f_conf_hdl);
conf_rtg_f_conf_fit_hdl = plot(m_conf_inferred_f_conf_inferred, m_conf_rtg_f_conf_inferred_fit);
fit_hdl_upgrade(conf_rtg_f_conf_fit_hdl);
xlabel('Confidence rating');
ylabel('Confidence');
legend_size(pSize);