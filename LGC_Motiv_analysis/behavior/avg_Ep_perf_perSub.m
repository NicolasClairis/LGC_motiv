function [latency, AUC, forcePeak, AUC_overshoot, AUC_N, forcePeak_N, AUC_overshoot_N] = avg_Ep_perf_perSub(study_nm, sub_nm, condition)
% [latency, AUC, forcePeak, AUC_overshoot, AUC_N, forcePeak_N, AUC_overshoot_N] = avg_Ep_perf_perSub(study_nm, sub_nm, condition)
% avg_Ep_perf_perSub will extract average physical performance
% characteristics for the subject defined in input. Each variable can then
% be dissected according to effort level chosen, reward level chosen, run,
% etc.
%
% INPUTS
% study_nm: study name
%
% sub_nm: subject name
%
% condition: will define which runs should be considered
%
% OUTPUTS
% latency: latency to start squeezing
%
% AUC: area under the curve for force performed
%
% forcePeak: peak force (in % of Fmax)
%
% AUC_overshoot: area under the curve overshoot above the threshold of 
% force to produce to do the trial
%
% AUC_N: area under the curve (based on force in Newtons)
%
% forcePeak_N: peak force in Newtons
%
% AUC_overshoot_N: area under the curve overshoot above the threshold of 
% force to produce to do the trial (based on force in Newtons)

%% working directories
% computerRoot = LGCM_root_paths;
computerRoot = ['E:',filesep];
dataRoot = [computerRoot, filesep, study_nm, filesep];

%% extract effort produced
Ech_levels = 0:3;
n_hE_levels = 3;
nTrialsPerRun = 54;
n_Ep_runs = 2;
n_total_Ep_trials = nTrialsPerRun*n_Ep_runs;
[latency_perTrial,...
    AUC_perTrial, forcePeak_perTrial, AUC_overshoot_perTrial,...
    AUC_N_perTrial, forcePeak_N_perTrial, AUC_overshoot_N_perTrial,...
    RP_perTrial, monetary_inc_perTrial,...
    Ech_perTrial, hE_level_perTrial, choice_highE_perTrial] = deal(NaN(1,n_total_Ep_trials));

task_fullName = 'physical';
subBehaviorFolder = [dataRoot, filesep, 'CID',sub_nm, filesep, 'behavior',filesep];
runs = runs_definition(study_nm, sub_nm, condition);
for iRun = 1:runs.nb_runs.Ep
    jRun = runs.Ep.runsToKeep(iRun);
    run_nm = num2str(jRun);
    switch jRun
        case {1,2}
            kRun = 1;
        case {3,4}
            kRun = 2;
    end
    run_nm_bis = ['run',num2str(kRun)];
    run_trials_idx = (1:nTrialsPerRun) + nTrialsPerRun.*(kRun - 1);
    
    % extract force
    [latency_tmp,...
        AUC_tmp, forcePeak_tmp, AUC_overshoot_tmp,...
        AUC_N_tmp, forcePeak_N_tmp, AUC_overshoot_N_tmp] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
    % extract info for all trials across all runs
    latency_perTrial(run_trials_idx) = latency_tmp.allTrials;
    AUC_perTrial(run_trials_idx) = AUC_tmp.allTrials;
    forcePeak_perTrial(run_trials_idx) = forcePeak_tmp.allTrials;
    AUC_overshoot_perTrial(run_trials_idx) = AUC_overshoot_tmp.allTrials;
    AUC_N_perTrial(run_trials_idx) = AUC_N_tmp.allTrials;
    forcePeak_N_perTrial(run_trials_idx) = forcePeak_N_tmp.allTrials;
    AUC_overshoot_N_perTrial(run_trials_idx) = AUC_overshoot_N_tmp.allTrials;
    % extract other variables of interest for spliting
    [Ech_perTrial(run_trials_idx)] = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
    [~, monetary_inc_perTrial(run_trials_idx)] = extract_money_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
    [choice_highE_perTrial(run_trials_idx)] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
    [hE_level_perTrial(run_trials_idx)] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
    [RP_perTrial(run_trials_idx)] = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
    
    % extract trials separately for each run
    latency.(run_nm_bis).allTrials = latency_tmp.allTrials;
    AUC.(run_nm_bis).allTrials = AUC_tmp.allTrials;
    forcePeak.(run_nm_bis).allTrials = forcePeak_tmp.allTrials;
    AUC_overshoot.(run_nm_bis).allTrials = AUC_overshoot_tmp.allTrials;
    AUC_N.(run_nm_bis).allTrials = AUC_N_tmp.allTrials;
    forcePeak_N.(run_nm_bis).allTrials = forcePeak_N_tmp.allTrials;
    AUC_overshoot_N.(run_nm_bis).allTrials = AUC_overshoot_N_tmp.allTrials;
    % extract trials split
    latency.(run_nm_bis).per_Ech = latency_tmp.per_Ech;
    AUC.(run_nm_bis).per_Ech = AUC_tmp.per_Ech;
    forcePeak.(run_nm_bis).per_Ech = forcePeak_tmp.per_Ech;
    AUC_overshoot.(run_nm_bis).per_Ech = AUC_overshoot_tmp.per_Ech;
    AUC_N.(run_nm_bis).per_Ech = AUC_N_tmp.per_Ech;
    forcePeak_N.(run_nm_bis).per_Ech = forcePeak_N_tmp.per_Ech;
    AUC_overshoot_N.(run_nm_bis).per_Ech = AUC_overshoot_N_tmp.per_Ech;
    % extract trials split per Effort level*choice
    latency.(run_nm_bis).per_hE.choice_lowE = latency_tmp.per_hE.choice_lowE;
    AUC.(run_nm_bis).per_hE.choice_lowE = AUC_tmp.per_hE.choice_lowE;
    forcePeak.(run_nm_bis).per_hE.choice_lowE = forcePeak_tmp.per_hE.choice_lowE;
    AUC_overshoot.(run_nm_bis).per_hE.choice_lowE = AUC_overshoot_tmp.per_hE.choice_lowE;
    AUC_N.(run_nm_bis).per_hE.choice_lowE = AUC_N_tmp.per_hE.choice_lowE;
    forcePeak_N.(run_nm_bis).per_hE.choice_lowE = forcePeak_N_tmp.per_hE.choice_lowE;
    AUC_overshoot_N.(run_nm_bis).per_hE.choice_lowE = AUC_overshoot_N_tmp.per_hE.choice_lowE;
    latency.(run_nm_bis).per_hE.choice_highE = latency_tmp.per_hE.choice_highE;
    AUC.(run_nm_bis).per_hE.choice_highE = AUC_tmp.per_hE.choice_highE;
    forcePeak.(run_nm_bis).per_hE.choice_highE = forcePeak_tmp.per_hE.choice_highE;
    AUC_overshoot.(run_nm_bis).per_hE.choice_highE = AUC_overshoot_tmp.per_hE.choice_highE;
    AUC_N.(run_nm_bis).per_hE.choice_highE = AUC_N_tmp.per_hE.choice_highE;
    forcePeak_N.(run_nm_bis).per_hE.choice_highE = forcePeak_N_tmp.per_hE.choice_highE;
    AUC_overshoot_N.(run_nm_bis).per_hE.choice_highE = AUC_overshoot_N_tmp.per_hE.choice_highE;
    
end % run loop

%% average data per variable of interest
latency.allRuns.allTrials           = mean(latency_perTrial,2,'omitnan');
AUC.allRuns.allTrials               = mean(AUC_perTrial,2,'omitnan');
forcePeak.allRuns.allTrials         = mean(forcePeak_perTrial,2,'omitnan');
AUC_overshoot.allRuns.allTrials     = mean(AUC_overshoot_perTrial,2,'omitnan');
AUC_N.allRuns.allTrials             = mean(AUC_N_perTrial,2,'omitnan');
forcePeak_N.allRuns.allTrials       = mean(forcePeak_N_perTrial,2,'omitnan');
AUC_overshoot_N.allRuns.allTrials   = mean(AUC_overshoot_N_perTrial,2,'omitnan');

%% split per effort chosen
for iEch = Ech_levels
    Ech_idx = Ech_perTrial == iEch;
    Ech_nm = ['Ech',num2str(iEch)];
    
    latency.allRuns.(Ech_nm) = mean(latency_perTrial(Ech_idx),2,'omitnan');
    AUC.allRuns.(Ech_nm) = mean(AUC_perTrial(Ech_idx),2,'omitnan');
    forcePeak.allRuns.(Ech_nm) = mean(forcePeak_perTrial(Ech_idx),2,'omitnan');
    AUC_overshoot.allRuns.(Ech_nm) = mean(AUC_overshoot_perTrial(Ech_idx),2,'omitnan');
    AUC_N.allRuns.(Ech_nm) = mean(AUC_N_perTrial(Ech_idx),2,'omitnan');
    forcePeak_N.allRuns.(Ech_nm) = mean(forcePeak_N_perTrial(Ech_idx),2,'omitnan');
    AUC_overshoot_N.allRuns.(Ech_nm) = mean(AUC_overshoot_N_perTrial(Ech_idx),2,'omitnan');
end % effort level chosen loop

%% split per level of money chosen (pooling R and P)
inc_ch_levels = 0:3;
for iInc = inc_ch_levels
    inc_idx = monetary_inc_perTrial == iInc;
    inc_nm = ['inc',num2str(iInc)];
    
    latency.allRuns.(inc_nm) = mean(latency_perTrial(inc_idx),2,'omitnan');
    AUC.allRuns.(inc_nm) = mean(AUC_perTrial(inc_idx),2,'omitnan');
    forcePeak.allRuns.(inc_nm) = mean(forcePeak_perTrial(inc_idx),2,'omitnan');
    AUC_overshoot.allRuns.(inc_nm) = mean(AUC_overshoot_perTrial(inc_idx),2,'omitnan');
    AUC_N.allRuns.(inc_nm) = mean(AUC_N_perTrial(inc_idx),2,'omitnan');
    forcePeak_N.allRuns.(inc_nm) = mean(forcePeak_N_perTrial(inc_idx),2,'omitnan');
    AUC_overshoot_N.allRuns.(inc_nm) = mean(AUC_overshoot_N_perTrial(inc_idx),2,'omitnan');
end

%% split per reward/punishment
for iRP = 1:2
    switch iRP
        case 1
            RP_nm = 'R';
            RP_idx = RP_perTrial == 1;
        case 2
            RP_nm = 'P';
            RP_idx = RP_perTrial == -1;
    end
    latency.allRuns.(RP_nm) = mean(latency_perTrial(RP_idx),2,'omitnan');
    AUC.allRuns.(RP_nm) = mean(AUC_perTrial(RP_idx),2,'omitnan');
    forcePeak.allRuns.(RP_nm) = mean(forcePeak_perTrial(RP_idx),2,'omitnan');
    AUC_overshoot.allRuns.(RP_nm) = mean(AUC_overshoot_perTrial(RP_idx),2,'omitnan');
    AUC_N.allRuns.(RP_nm) = mean(AUC_N_perTrial(RP_idx),2,'omitnan');
    forcePeak_N.allRuns.(RP_nm) = mean(forcePeak_N_perTrial(RP_idx),2,'omitnan');
    AUC_overshoot_N.allRuns.(RP_nm) = mean(AUC_overshoot_N_perTrial(RP_idx),2,'omitnan');
end % R/P trials

%% split per effort level proposed * choice made
choice_hE_idx = choice_highE_perTrial == 1;
choice_lE_idx = choice_highE_perTrial == 0;
[latency.allRuns.per_hE.choice_lowE,...
    AUC.allRuns.per_hE.choice_lowE,...
    forcePeak.allRuns.per_hE.choice_lowE,...
    AUC_overshoot.allRuns.per_hE.choice_lowE,...
    AUC_N.allRuns.per_hE.choice_lowE,...
    forcePeak_N.allRuns.per_hE.choice_lowE,...
    AUC_overshoot_N.allRuns.per_hE.choice_lowE,...
    latency.allRuns.per_hE.choice_highE,...
    AUC.allRuns.per_hE.choice_highE,...
    forcePeak.allRuns.per_hE.choice_highE,...
    AUC_overshoot.allRuns.per_hE.choice_highE,...
    AUC_N.allRuns.per_hE.choice_highE,...
    forcePeak_N.allRuns.per_hE.choice_highE,...
    AUC_overshoot_N.allRuns.per_hE.choice_highE] = deal(NaN(1,n_hE_levels));
for iE = 1:n_hE_levels
    hE_idx = hE_level_perTrial == iE;
    hEch_idx = (hE_idx.*choice_hE_idx) == 1;
    lEch_idx = (hE_idx.*choice_lE_idx) == 1;
    % extract data
    % low effort chosen
    latency.allRuns.per_hE.choice_lowE(iE) = mean(latency_perTrial(lEch_idx),2,'omitnan');
    AUC.allRuns.per_hE.choice_lowE(iE) = mean(AUC_perTrial(lEch_idx),2,'omitnan');
    forcePeak.allRuns.per_hE.choice_lowE(iE) = mean(forcePeak_perTrial(lEch_idx),2,'omitnan');
    AUC_overshoot.allRuns.per_hE.choice_lowE(iE) = mean(AUC_overshoot_perTrial(lEch_idx),2,'omitnan');
    AUC_N.allRuns.per_hE.choice_lowE(iE) = mean(AUC_N_perTrial(lEch_idx),2,'omitnan');
    forcePeak_N.allRuns.per_hE.choice_lowE(iE) = mean(forcePeak_N_perTrial(lEch_idx),2,'omitnan');
    AUC_overshoot_N.allRuns.per_hE.choice_lowE(iE) = mean(AUC_overshoot_N_perTrial(lEch_idx),2,'omitnan');
    % high effort chosen
    latency.allRuns.per_hE.choice_highE(iE) = mean(latency_perTrial(hEch_idx),2,'omitnan');
    AUC.allRuns.per_hE.choice_highE(iE) = mean(AUC_perTrial(hEch_idx),2,'omitnan');
    forcePeak.allRuns.per_hE.choice_highE(iE) = mean(forcePeak_perTrial(hEch_idx),2,'omitnan');
    AUC_overshoot.allRuns.per_hE.choice_highE(iE) = mean(AUC_overshoot_perTrial(hEch_idx),2,'omitnan');
    AUC_N.allRuns.per_hE.choice_highE(iE) = mean(AUC_N_perTrial(hEch_idx),2,'omitnan');
    forcePeak_N.allRuns.per_hE.choice_highE(iE) = mean(forcePeak_N_perTrial(hEch_idx),2,'omitnan');
    AUC_overshoot_N.allRuns.per_hE.choice_highE(iE) = mean(AUC_overshoot_N_perTrial(hEch_idx),2,'omitnan');
end % effort level

%% extract intercept and slope for each variable
% latency high effort choice
b_latency_highE = glmfit(hE_level_perTrial(choice_hE_idx),...
    latency_perTrial(choice_hE_idx), 'normal');
latency.intercept.choice_highE = b_latency_highE(1);
latency.slope.choice_highE = b_latency_highE(2);
% latency low effort choice
b_latency_lowE = glmfit(hE_level_perTrial(choice_lE_idx),...
    latency_perTrial(choice_lE_idx), 'normal');
latency.intercept.choice_lowE = b_latency_lowE(1);
latency.slope.choice_lowE = b_latency_lowE(2);

% AUC high effort choice
b_AUC_highE = glmfit(hE_level_perTrial(choice_hE_idx),...
    AUC_perTrial(choice_hE_idx), 'normal');
AUC.intercept.choice_highE = b_AUC_highE(1);
AUC.slope.choice_highE = b_AUC_highE(2);
% AUC low effort choice
b_AUC_lowE = glmfit(hE_level_perTrial(choice_lE_idx),...
    AUC_perTrial(choice_lE_idx), 'normal');
AUC.intercept.choice_lowE = b_AUC_lowE(1);
AUC.slope.choice_lowE = b_AUC_lowE(2);

% forcePeak high effort choice
b_forcePeak_highE = glmfit(hE_level_perTrial(choice_hE_idx),...
    forcePeak_perTrial(choice_hE_idx), 'normal');
forcePeak.intercept.choice_highE = b_forcePeak_highE(1);
forcePeak.slope.choice_highE = b_forcePeak_highE(2);
% forcePeak low effort choice
b_forcePeak_lowE = glmfit(hE_level_perTrial(choice_lE_idx),...
    forcePeak_perTrial(choice_lE_idx), 'normal');
forcePeak.intercept.choice_lowE = b_forcePeak_lowE(1);
forcePeak.slope.choice_lowE = b_forcePeak_lowE(2);

% AUC_overshoot high effort choice
b_AUC_overshoot_highE = glmfit(hE_level_perTrial(choice_hE_idx),...
    AUC_overshoot_perTrial(choice_hE_idx), 'normal');
AUC_overshoot.intercept.choice_highE = b_AUC_overshoot_highE(1);
AUC_overshoot.slope.choice_highE = b_AUC_overshoot_highE(2);
% AUC_overshoot low effort choice
b_AUC_overshoot_lowE = glmfit(hE_level_perTrial(choice_lE_idx),...
    AUC_overshoot_perTrial(choice_lE_idx), 'normal');
AUC_overshoot.intercept.choice_lowE = b_AUC_overshoot_lowE(1);
AUC_overshoot.slope.choice_lowE = b_AUC_overshoot_lowE(2);

% AUC_N high effort choice
b_AUC_N_highE = glmfit(hE_level_perTrial(choice_hE_idx),...
    AUC_N_perTrial(choice_hE_idx), 'normal');
AUC_N.intercept.choice_highE = b_AUC_N_highE(1);
AUC_N.slope.choice_highE = b_AUC_N_highE(2);
% AUC_N low effort choice
b_AUC_N_lowE = glmfit(hE_level_perTrial(choice_lE_idx),...
    AUC_N_perTrial(choice_lE_idx), 'normal');
AUC_N.intercept.choice_lowE = b_AUC_N_lowE(1);
AUC_N.slope.choice_lowE = b_AUC_N_lowE(2);

% forcePeak_N high effort choice
b_forcePeak_N_highE = glmfit(hE_level_perTrial(choice_hE_idx),...
    forcePeak_N_perTrial(choice_hE_idx), 'normal');
forcePeak_N.intercept.choice_highE = b_forcePeak_N_highE(1);
forcePeak_N.slope.choice_highE = b_forcePeak_N_highE(2);
% forcePeak_N low effort choice
b_forcePeak_N_lowE = glmfit(hE_level_perTrial(choice_lE_idx),...
    forcePeak_N_perTrial(choice_lE_idx), 'normal');
forcePeak_N.intercept.choice_lowE = b_forcePeak_N_lowE(1);
forcePeak_N.slope.choice_lowE = b_forcePeak_N_lowE(2);

% AUC_overshoot_N high effort choice
b_AUC_overshoot_N_highE = glmfit(hE_level_perTrial(choice_hE_idx),...
    AUC_overshoot_N_perTrial(choice_hE_idx), 'normal');
AUC_overshoot_N.intercept.choice_highE = b_AUC_overshoot_N_highE(1);
AUC_overshoot_N.slope.choice_highE = b_AUC_overshoot_N_highE(2);
% AUC_overshoot_N low effort choice
b_AUC_overshoot_N_lowE = glmfit(hE_level_perTrial(choice_lE_idx),...
    AUC_overshoot_N_perTrial(choice_lE_idx), 'normal');
AUC_overshoot_N.intercept.choice_lowE = b_AUC_overshoot_N_lowE(1);
AUC_overshoot_N.slope.choice_lowE = b_AUC_overshoot_N_lowE(2);

end % function