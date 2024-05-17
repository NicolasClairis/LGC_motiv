function[r_corr, pval, R2] = correl_Ech_SVch_DT(z_perRun)
% [r_corr, pval, R2] = correl_Ech_SVch_DT()
% correl_Ech_SVch_DT will test the correlation between the effort chosen (Ech),
% subjective value of the chosen option (SVch) and deliberation times (DT)
% across subjects in the fMRI GLM to see how correlated these regressors are.
%
% INPUTS
% z_perRun: zscore variables within each run?
%
% OUTPUTS
% r_corr: structure with the correlation coefficients
%
% pval: structure with p.value for each correlation
%
% R2: structure with R2 coefficients indicating proportion of variance
% explained

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;

%% model selection
[~, mdl_n_nm] = which_bayesian_mdl_n;
bayesianMdl_nm = ['mdl_',mdl_n_nm];

%% working directories
% computerRoot = LGCM_root_paths();
computerRoot = ['E:',filesep];
% scripts_folder = fullfile(computer_root,'GitHub','LGC_motiv','LGC_Motiv_analysis','fMRI');
% addpath(scripts_folder);
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computerRoot,'fMRI_pilots');
    case 'study1'
        root = fullfile(computerRoot,'study1');
    case 'study2'
        root = fullfile(computerRoot,'study2');
end
gitResultsFolder = [fullfile('C:','Users','clairis','Desktop',...
            'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];

%% initialize variables of interest
nTrialsPerRun = 54;
nRunsPerTask = 2;
nTotalRuns = 4;
nTotalTrials = nTrialsPerRun.*nTotalRuns;
nTrialsPerTask = nTrialsPerRun.*nRunsPerTask;
[Ech.all, SVch.all, DT.all] = deal(NaN(nTotalTrials, NS));
[Ech.Ep, SVch.Ep, DT.Ep,...
    Ech.Em, SVch.Em, DT.Em] = deal(NaN(nTrialsPerTask, NS));
[r_corr.perSub.all.Ech_vs_SVch,...
    r_corr.perSub.all.Ech_vs_DT,...
    r_corr.perSub.all.SVch_vs_DT,...
    r_corr.perSub.Ep.Ech_vs_SVch,...
    r_corr.perSub.Ep.Ech_vs_DT,...
    r_corr.perSub.Ep.SVch_vs_DT,...
    r_corr.perSub.Em.Ech_vs_SVch,...
    r_corr.perSub.Em.Ech_vs_DT,...
    r_corr.perSub.Em.SVch_vs_DT,...
    pval.perSub.all.Ech_vs_SVch,...
    pval.perSub.all.Ech_vs_DT,...
    pval.perSub.all.SVch_vs_DT,...
    pval.perSub.Ep.Ech_vs_SVch,...
    pval.perSub.Ep.Ech_vs_DT,...
    pval.perSub.Ep.SVch_vs_DT,...
    pval.perSub.Em.Ech_vs_SVch,...
    pval.perSub.Em.Ech_vs_DT,...
    pval.perSub.Em.SVch_vs_DT] = deal(NaN(1, NS));

%% main parameter
if ~exist('z_perRun','var') || isempty(z_perRun)
    z_perRun = 1;
    warning(['z_perRun is set to ',num2str(z_perRun),' by default as not entered in inputs.']);
end

%% extract the data
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    % define working folders
    subj_folder             = [root, filesep, 'CID',sub_nm];
    subj_behavior_folder    = [subj_folder, filesep, 'behavior' filesep];
    
    %% define number of runs
    [~, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    % loop through runs and tasks
    for iRun = 1:n_runs
        % fix run index if some runs were removed
        [~, jRun, run_task_nm] = First_level_subRunFilter(study_nm, sub_nm, [], iRun, condition);
        run_nm = num2str(jRun);
        switch run_task_nm
            case 'Ep'
                run_task_fullName = 'physical';
                task_behavioral_id = 'physicalPerf';
            case 'Em'
                run_task_fullName = 'mental';
                task_behavioral_id = 'mentalE_perf';
        end
        
        % load variables
        currRunBehaviorFileName = ls([subj_behavior_folder,'*_session',run_nm,'_*_task.mat']);
        behavioralDataStruct = load([subj_behavior_folder, currRunBehaviorFileName]);
        % load E chosen
        Ech_tmp = behavioralDataStruct.(task_behavioral_id).E_chosen; % 0 for low E chosen and 1/2/3 for high effort chosen
        % load SV chosen
        [~,~,~,~,...
            ~, SVch_tmp] = extract_bayesian_mdl(gitResultsFolder, subj_behavior_folder,...
            sub_nm, run_nm, run_task_fullName, bayesianMdl_nm);
        % load DT
        T0 = behavioralDataStruct.onsets.T0;
        choiceOnsets            = behavioralDataStruct.(task_behavioral_id).onsets.choice - T0;
        dispChoiceOptionOnsets  = behavioralDataStruct.(task_behavioral_id).onsets.dispChoiceOptions - T0;
        DT_tmp = choiceOnsets - dispChoiceOptionOnsets;
        
        %% remove missed trials
        choiceMissedTrials = isnan(choiceOnsets);
        Ech_tmp(choiceMissedTrials) = NaN;
        SVch_tmp(choiceMissedTrials) = NaN;
        DT_tmp(choiceMissedTrials) = NaN;
        
        %% zscore per run
        if z_perRun == 1
            Ech_tmp = nanzscore(Ech_tmp);
            SVch_tmp = nanzscore(SVch_tmp);
            DT_tmp = nanzscore(DT_tmp);
        end %  zscore
        
        %% upload data in main structures
        % main variable
        general_trial_idx = (1:nTrialsPerRun) + nTrialsPerRun.*(jRun - 1);
        Ech.all(general_trial_idx, iS) = Ech_tmp;
        SVch.all(general_trial_idx, iS) = SVch_tmp;
        DT.all(general_trial_idx, iS) = DT_tmp;
        % upload task-specific variable
        switch jRun
            case {1,2}
                kRun = 1;
            case {3,4}
                kRun = 2;
        end
        task_trial_idx = (1:nTrialsPerRun) + nTrialsPerRun.*(kRun - 1);
        Ech.(run_task_nm)(task_trial_idx, iS) = Ech_tmp;
        SVch.(run_task_nm)(task_trial_idx, iS) = SVch_tmp;
        DT.(run_task_nm)(task_trial_idx, iS) = DT_tmp;
    end % loop over runs
    
    %% perform correlations
    % extract r
    % all trials
    ok_allTrials = ~isnan(Ech.all(:,iS).*SVch.all(:,iS));
    [r_corr.perSub.all.Ech_vs_SVch(iS), pval.perSub.all.Ech_vs_SVch(iS)] = corr(Ech.all(ok_allTrials,iS), SVch.all(ok_allTrials,iS));
    [r_corr.perSub.all.Ech_vs_DT(iS), pval.perSub.all.Ech_vs_DT(iS)] = corr(Ech.all(ok_allTrials,iS), DT.all(ok_allTrials,iS));
    [r_corr.perSub.all.SVch_vs_DT(iS), pval.perSub.all.SVch_vs_DT(iS)] = corr(SVch.all(ok_allTrials,iS), DT.all(ok_allTrials,iS));
    % Ep
    ok_EpTrials = ~isnan(Ech.Ep(:,iS).*SVch.Ep(:,iS));
    [r_corr.perSub.Ep.Ech_vs_SVch(iS), pval.perSub.Ep.Ech_vs_SVch(iS)] = corr(Ech.Ep(ok_EpTrials,iS), SVch.Ep(ok_EpTrials,iS));
    [r_corr.perSub.Ep.Ech_vs_DT(iS), pval.perSub.Ep.Ech_vs_DT(iS)] = corr(Ech.Ep(ok_EpTrials,iS), DT.Ep(ok_EpTrials,iS));
    [r_corr.perSub.Ep.SVch_vs_DT(iS), pval.perSub.Ep.SVch_vs_DT(iS)] = corr(SVch.Ep(ok_EpTrials,iS), DT.Ep(ok_EpTrials,iS));
    % Em
    ok_EmTrials = ~isnan(Ech.Em(:,iS).*SVch.Em(:,iS));
    [r_corr.perSub.Em.Ech_vs_SVch(iS), pval.perSub.Em.Ech_vs_SVch(iS)] = corr(Ech.Em(ok_EmTrials,iS), SVch.Em(ok_EmTrials,iS));
    [r_corr.perSub.Em.Ech_vs_DT(iS), pval.perSub.Em.Ech_vs_DT(iS)] = corr(Ech.Em(ok_EmTrials,iS), DT.Em(ok_EmTrials,iS));
    [r_corr.perSub.Em.SVch_vs_DT(iS), pval.perSub.Em.SVch_vs_DT(iS)] = corr(SVch.Em(ok_EmTrials,iS), DT.Em(ok_EmTrials,iS));
    
    % extract R2
    % all trials
    LM_tmp = fitlm(Ech.all(ok_allTrials,iS), SVch.all(ok_allTrials,iS));
    R2.perSub.all.Ech_vs_SVch(iS) = LM_tmp.Rsquared.Adjusted;
    LM_tmp = fitlm(Ech.all(ok_allTrials,iS), DT.all(ok_allTrials,iS));
    R2.perSub.all.Ech_vs_DT(iS) = LM_tmp.Rsquared.Adjusted;
    LM_tmp = fitlm(SVch.all(ok_allTrials,iS), DT.all(ok_allTrials,iS));
    R2.perSub.all.SVch_vs_DT(iS) = LM_tmp.Rsquared.Adjusted;
    % Ep
    LM_tmp = fitlm(Ech.Ep(ok_EpTrials,iS), SVch.Ep(ok_EpTrials,iS));
    R2.perSub.Ep.Ech_vs_SVch(iS) = LM_tmp.Rsquared.Adjusted;
    LM_tmp = fitlm(Ech.Ep(ok_EpTrials,iS), DT.Ep(ok_EpTrials,iS));
    R2.perSub.Ep.Ech_vs_DT(iS) = LM_tmp.Rsquared.Adjusted;
    LM_tmp = fitlm(SVch.Ep(ok_EpTrials,iS), DT.Ep(ok_EpTrials,iS));
    R2.perSub.Ep.SVch_vs_DT(iS) = LM_tmp.Rsquared.Adjusted;
    % Em
    LM_tmp = fitlm(Ech.Em(ok_EmTrials,iS), SVch.Em(ok_EmTrials,iS));
    R2.perSub.Em.Ech_vs_SVch(iS) = LM_tmp.Rsquared.Adjusted;
    LM_tmp = fitlm(Ech.Em(ok_EmTrials,iS), DT.Em(ok_EmTrials,iS));
    R2.perSub.Em.Ech_vs_DT(iS) = LM_tmp.Rsquared.Adjusted;
    LM_tmp = fitlm(SVch.Em(ok_EmTrials,iS), DT.Em(ok_EmTrials,iS));
    R2.perSub.Em.SVch_vs_DT(iS) = LM_tmp.Rsquared.Adjusted;
end % subject loop

%% average across subjects
% across tasks
[mEch.all, semEch.all] = mean_sem_sd(Ech.all,2); % note: average done on all trials, would be more elegant to bin depending on the Ech/SVch/DT levels
[mSVch.all, semSVch.all] = mean_sem_sd(SVch.all,2);
[mDT.all, semDT.all] = mean_sem_sd(DT.all,2);
r_corr.aSubs.all.Ech_vs_SVch = mean(r_corr.perSub.all.Ech_vs_SVch, 2);
r_corr.aSubs.all.Ech_vs_DT = mean(r_corr.perSub.all.Ech_vs_DT, 2);
r_corr.aSubs.all.SVch_vs_DT = mean(r_corr.perSub.all.SVch_vs_DT, 2);
R2.aSubs.all.Ech_vs_SVch = mean(R2.perSub.all.Ech_vs_SVch, 2);
R2.aSubs.all.Ech_vs_DT = mean(R2.perSub.all.Ech_vs_DT, 2);
R2.aSubs.all.SVch_vs_DT = mean(R2.perSub.all.SVch_vs_DT, 2);
% Ep
[mEch.Ep, semEch.Ep] = mean_sem_sd(Ech.Ep,2);
[mSVch.Ep, semSVch.Ep] = mean_sem_sd(SVch.Ep,2);
[mDT.Ep, semDT.Ep] = mean_sem_sd(DT.Ep,2);
r_corr.aSubs.Ep.Ech_vs_SVch = mean(r_corr.perSub.Ep.Ech_vs_SVch, 2);
r_corr.aSubs.Ep.Ech_vs_DT = mean(r_corr.perSub.Ep.Ech_vs_DT, 2);
r_corr.aSubs.Ep.SVch_vs_DT = mean(r_corr.perSub.Ep.SVch_vs_DT, 2);
R2.aSubs.Ep.Ech_vs_SVch = mean(R2.perSub.Ep.Ech_vs_SVch, 2);
R2.aSubs.Ep.Ech_vs_DT = mean(R2.perSub.Ep.Ech_vs_DT, 2);
R2.aSubs.Ep.SVch_vs_DT = mean(R2.perSub.Ep.SVch_vs_DT, 2);
% Em
[mEch.Em, semEch.Em] = mean_sem_sd(Ech.Em,2);
[mSVch.Em, semSVch.Em] = mean_sem_sd(SVch.Em,2);
[mDT.Em, semDT.Em] = mean_sem_sd(DT.Em,2);
r_corr.aSubs.Em.Ech_vs_SVch = mean(r_corr.perSub.Em.Ech_vs_SVch, 2);
r_corr.aSubs.Em.Ech_vs_DT = mean(r_corr.perSub.Em.Ech_vs_DT, 2);
r_corr.aSubs.Em.SVch_vs_DT = mean(r_corr.perSub.Em.SVch_vs_DT, 2);
R2.aSubs.Em.Ech_vs_SVch = mean(R2.perSub.Em.Ech_vs_SVch, 2);
R2.aSubs.Em.Ech_vs_DT = mean(R2.perSub.Em.Ech_vs_DT, 2);
R2.aSubs.Em.SVch_vs_DT = mean(R2.perSub.Em.SVch_vs_DT, 2);

% test if correlation coefficients are significantly different from zero
% across tasks
[~,pval.aSubs.all.Ech_vs_SVch] = ttest(r_corr.perSub.all.Ech_vs_SVch);
[~,pval.aSubs.all.Ech_vs_DT] = ttest(r_corr.perSub.all.Ech_vs_DT);
[~,pval.aSubs.all.SVch_vs_DT] = ttest(r_corr.perSub.all.SVch_vs_DT);
% Ep
[~,pval.aSubs.Ep.Ech_vs_SVch] = ttest(r_corr.perSub.Ep.Ech_vs_SVch);
[~,pval.aSubs.Ep.Ech_vs_DT] = ttest(r_corr.perSub.Ep.Ech_vs_DT);
[~,pval.aSubs.Ep.SVch_vs_DT] = ttest(r_corr.perSub.Ep.SVch_vs_DT);
% Em
[~,pval.aSubs.Em.Ech_vs_SVch] = ttest(r_corr.perSub.Em.Ech_vs_SVch);
[~,pval.aSubs.Em.Ech_vs_DT] = ttest(r_corr.perSub.Em.Ech_vs_DT);
[~,pval.aSubs.Em.SVch_vs_DT] = ttest(r_corr.perSub.Em.SVch_vs_DT);

% show data
% across all trials
fig;

subplot(1,3,1);
scat_hdl = scatter(mEch.all, mSVch.all);
scat_hdl_upgrade(scat_hdl);
xlabel('Ech');
ylabel('SVch');
place_r_and_pval(r_corr.aSubs.all.Ech_vs_SVch,...
    pval.aSubs.all.Ech_vs_SVch);

subplot(1,3,2);
scat_hdl = scatter(mEch.all, mDT.all);
scat_hdl_upgrade(scat_hdl);
xlabel('Ech');
ylabel('DT');
place_r_and_pval(r_corr.aSubs.all.Ech_vs_DT,...
    pval.aSubs.all.Ech_vs_DT);

subplot(1,3,3);
scat_hdl = scatter(mSVch.all, mDT.all);
scat_hdl_upgrade(scat_hdl);
xlabel('SVch');
ylabel('DT');
place_r_and_pval(r_corr.aSubs.all.SVch_vs_DT,...
    pval.aSubs.all.SVch_vs_DT);

% Ep
fig;

subplot(1,3,1);
scat_hdl = scatter(mEch.Ep, mSVch.Ep);
scat_hdl_upgrade(scat_hdl);
xlabel('Ech');
ylabel('SVch');
place_r_and_pval(r_corr.aSubs.Ep.Ech_vs_SVch,...
    pval.aSubs.Ep.Ech_vs_SVch);

subplot(1,3,2);
scat_hdl = scatter(mEch.Ep, mDT.Ep);
scat_hdl_upgrade(scat_hdl);
xlabel('Ech');
ylabel('DT');
place_r_and_pval(r_corr.aSubs.Ep.Ech_vs_DT,...
    pval.aSubs.Ep.Ech_vs_DT);

subplot(1,3,3);
scat_hdl = scatter(mSVch.Ep, mDT.Ep);
scat_hdl_upgrade(scat_hdl);
xlabel('SVch');
ylabel('DT');
place_r_and_pval(r_corr.aSubs.Ep.SVch_vs_DT,...
    pval.aSubs.Ep.SVch_vs_DT);

% Em
fig;

subplot(1,3,1);
scat_hdl = scatter(mEch.Em, mSVch.Em);
scat_hdl_upgrade(scat_hdl);
xlabel('Ech');
ylabel('SVch');
place_r_and_pval(r_corr.aSubs.Em.Ech_vs_SVch,...
    pval.aSubs.Em.Ech_vs_SVch);

subplot(1,3,2);
scat_hdl = scatter(mEch.Em, mDT.Em);
scat_hdl_upgrade(scat_hdl);
xlabel('Ech');
ylabel('DT');
place_r_and_pval(r_corr.aSubs.Em.Ech_vs_DT,...
    pval.aSubs.Em.Ech_vs_DT);

subplot(1,3,3);
scat_hdl = scatter(mSVch.Em, mDT.Em);
scat_hdl_upgrade(scat_hdl);
xlabel('SVch');
ylabel('DT');
place_r_and_pval(r_corr.aSubs.Em.SVch_vs_DT,...
    pval.aSubs.Em.SVch_vs_DT);
end % function