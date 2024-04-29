%% extract the level of activity of a given ROI in function of each effort level chosen
% and split subjects with a median split based on the metabolites. the same
% script will also look at the choice depending on dmPFC activity and
% effort level proposed for the high effort option
%
% Same as ROI_f_Rch_Pch__metabolites_mSplit_per_task.m but spliting data
% between mental and physical effort task

%% general parameters
% display figures?
dispFig = true;

%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directory
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% subject definition
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load ROI
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);
% define which ROI, and which time period is of interest to you
% define ROI
ROI_names = fieldnames(ROI_trial_b_trial);
ROI_subList = ROI_trial_b_trial.subject_id;
ROI_names(strcmp(ROI_names,'subject_id')) = [];
if length(ROI_names) > 1
    error(['There should be only 1 ROI selected, not ',num2str(length(ROI_names))])
else
    fMRI_ROI_nm = ROI_names;
end
ROI_short_nm = inputdlg('ROI short name?');
ROI_short_nm = ROI_short_nm{1};
% define task
task_names = {'Ep','Em'};
nTasks = length(task_names);
% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(fMRI_ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};

%% load behavior 
nRunsPerTask = 2;
nTrialsPerRun = 54;
nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
% main variables of interest
% reward
Rchosen_possible = 0:3;
n_Rch = length(Rchosen_possible);
hR_levels = 1:3;
n_hR_levels = length(hR_levels);
% punishment
Pchosen_possible = 0:3;
n_Pch = length(Pchosen_possible);
hP_levels = 1:3;
n_hP_levels = length(hP_levels);
% incentives
Ichosen_possible = 0:3;
n_Ich = length(Ichosen_possible);
hI_levels = 1:3;
n_hI_levels = length(hI_levels);

% data for all trials
[fMRI_ROI.Ep.allTrials, choice_hE.Ep.allTrials,...
    RP.Ep.allTrials,...
    hR_level.Ep.allTrials, Rchosen.Ep.allTrials,...
    hP_level.Ep.allTrials, Pchosen.Ep.allTrials,...
    hI_level.Ep.allTrials, Ichosen.Ep.allTrials,...
    fMRI_ROI.Em.allTrials, choice_hE.Em.allTrials,...
    RP.Em.allTrials,...
    hR_level.Em.allTrials, Rchosen.Em.allTrials,...
    hP_level.Em.allTrials, Pchosen.Em.allTrials,...
    hI_level.Em.allTrials, Ichosen.Em.allTrials] = deal(NaN(nTrialsPerTask, NS));
% data split by effort chosen
[fMRI_ROI.Ep.Rch, fMRI_ROI.Em.Rch] = deal(NaN(n_Rch,NS));
[fMRI_ROI.Ep.Pch, fMRI_ROI.Em.Pch] = deal(NaN(n_Pch,NS));
[fMRI_ROI.Ep.Ich, fMRI_ROI.Em.Ich] = deal(NaN(n_Ich,NS));
% data split by level of incentive proposed for the high effort option
[fMRI_ROI.Ep.hR_level, choice_hE.Ep.hR_level,...
    fMRI_ROI.Em.hR_level, choice_hE.Em.hR_level] = deal(NaN(n_hR_levels, NS));
[fMRI_ROI.Ep.hP_level, choice_hE.Ep.hP_level,...
    fMRI_ROI.Em.hP_level, choice_hE.Em.hP_level] = deal(NaN(n_hP_levels, NS));
[fMRI_ROI.Ep.hI_level, choice_hE.Ep.hI_level,...
    fMRI_ROI.Em.hI_level, choice_hE.Em.hI_level] = deal(NaN(n_hI_levels, NS));
[fMRI_ROI.Ep.choice_low.hR_level,...
    fMRI_ROI.Ep.choice_high.hR_level,...
    fMRI_ROI.Ep.choice_low.hP_level,...
    fMRI_ROI.Ep.choice_high.hP_level,...
    fMRI_ROI.Ep.choice_low.hI_level,...
    fMRI_ROI.Ep.choice_high.hI_level,...
    fMRI_ROI.Em.choice_low.hR_level,...
    fMRI_ROI.Em.choice_high.hR_level,...
    fMRI_ROI.Em.choice_low.hP_level,...
    fMRI_ROI.Em.choice_high.hP_level,...
    fMRI_ROI.Em.choice_low.hI_level,...
    fMRI_ROI.Em.choice_high.hI_level] = deal(NaN(n_hR_levels, NS));

% slopes
[b_fMRI_f_Rch.Ep.allSubs,...
    b_fMRI_f_R.Ep.lE_chosen.allSubs,...
    b_fMRI_f_R.Ep.hE_chosen.allSubs,...
    b_fMRI_f_Rch.Em.allSubs,...
    b_fMRI_f_R.Em.lE_chosen.Em.allSubs,...
    b_fMRI_f_R.Em.hE_chosen.Em.allSubs,...
    b_fMRI_f_Pch.Ep.allSubs,...
    b_fMRI_f_P.Ep.lE_chosen.allSubs,...
    b_fMRI_f_P.Ep.hE_chosen.allSubs,...
    b_fMRI_f_Pch.Em.allSubs,...
    b_fMRI_f_P.Em.lE_chosen.Em.allSubs,...
    b_fMRI_f_P.Em.hE_chosen.Em.allSubs,...
    b_fMRI_f_Ich.Ep.allSubs,...
    b_fMRI_f_I.Ep.lE_chosen.allSubs,...
    b_fMRI_f_I.Ep.hE_chosen.allSubs,...
    b_fMRI_f_Ich.Em.allSubs,...
    b_fMRI_f_I.Em.lE_chosen.Em.allSubs,...
    b_fMRI_f_I.Em.hE_chosen.Em.allSubs] = deal(NaN(2,NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];

    % extract runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    jRun = 0;
    for iRun = 1:length(okRuns)
        kRun = okRuns(iRun);
        run_nm = num2str(kRun);
        jRun = jRun + 1;
        task_nm_tmp = taskNames{jRun};
        switch kRun
            case {1,2}
                runTrials_idx = 1:nTrialsPerRun;
            case {3,4}
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun;
        end
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        % define which task session it is
        switch kRun
            case {1,2}
                run_nm_bis = ['run',num2str(1)];
            case {3,4}
                run_nm_bis = ['run',num2str(2)];
        end
        
        %% extract relevant variables
        % R/P trials
        RP_var_tmp = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        R_trials = RP_var_tmp == 1;
        P_trials = RP_var_tmp == -1;
        RP.(task_nm_tmp).allTrials(runTrials_idx,iS) = RP_var_tmp;
        % choice
        choice_hE.(task_nm_tmp).allTrials(runTrials_idx,iS) = extract_choice_hE(subBehaviorFolder,sub_nm,run_nm,task_fullName);
        % extract R and P level for high effort option and for choice
        hI_tmp = extract_hI_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        hI_level.(task_nm_tmp).allTrials(runTrials_idx,iS) = hI_tmp;
        hR_level.(task_nm_tmp).allTrials(runTrials_idx(R_trials),iS) = hI_tmp(R_trials);
        hP_level.(task_nm_tmp).allTrials(runTrials_idx(P_trials),iS) = hI_tmp(P_trials);
        Ich_tmp = extract_I_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        I_chosen.(task_nm_tmp).allTrials(runTrials_idx,iS) = Ich_tmp;
        Rchosen.(task_nm_tmp).allTrials(runTrials_idx(R_trials),iS) = Ich_tmp(R_trials);
        Pchosen.(task_nm_tmp).allTrials(runTrials_idx(P_trials),iS) = Ich_tmp(P_trials);
        
        %% extract fMRI ROI mediator
        fMRI_ROI.(task_nm_tmp).allTrials(runTrials_idx, iS) = ROI_trial_b_trial.(fMRI_ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
    end % run loop
    
    %% split the data
    for iT = 1:nTasks
        task_nm = task_names{iT};
        %% extract data per incentive level and effort chosen
        for iInc_level = 1:n_hI_levels
            jI_idx = hI_level.(task_nm).allTrials(:, iS) == iInc_level;
            % fMRI ROI depending on effort level
            fMRI_ROI.(task_nm).hI_level(iInc_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jI_idx, iS),1,'omitnan');
            choice_hE.(task_nm).hI_level(iInc_level,iS) = mean(choice_hE.(task_nm).allTrials(jI_idx, iS),1,'omitnan');

            % fMRI ROI depending on effort level AND choice (like in Kurniawan
            % et al, 2021)
            jI_highEchosen_idx = (hI_level.(task_nm).allTrials(:, iS) == iInc_level).*(choice_hE.(task_nm).allTrials(:, iS) == 1) == 1;
            jI_lowEchosen_idx = (hI_level.(task_nm).allTrials(:, iS) == iInc_level).*(choice_hE.(task_nm).allTrials(:, iS) == 0) == 1;
            fMRI_ROI.(task_nm).choice_high.hI_level(iInc_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jI_highEchosen_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm).choice_low.hI_level(iInc_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jI_lowEchosen_idx, iS),1,'omitnan');
        end % incentive level
        
        for iR_level = 1:n_hR_levels
            jR_idx = hR_level.(task_nm).allTrials(:, iS) == iR_level;
            % fMRI ROI depending on effort level
            fMRI_ROI.(task_nm).hR_level(iR_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jR_idx, iS),1,'omitnan');
            choice_hE.(task_nm).hR_level(iR_level,iS) = mean(choice_hE.(task_nm).allTrials(jR_idx, iS),1,'omitnan');

            % fMRI ROI depending on effort level AND choice (like in Kurniawan
            % et al, 2021)
            jR_highEchosen_idx = (hR_level.(task_nm).allTrials(:, iS) == iR_level).*(choice_hE.(task_nm).allTrials(:, iS) == 1) == 1;
            jR_lowEchosen_idx = (hR_level.(task_nm).allTrials(:, iS) == iR_level).*(choice_hE.(task_nm).allTrials(:, iS) == 0) == 1;
            fMRI_ROI.(task_nm).choice_high.hR_level(iR_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jR_highEchosen_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm).choice_low.hR_level(iR_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jR_lowEchosen_idx, iS),1,'omitnan');
        end % incentive level
        
        for iP_level = 1:n_hP_levels
            jP_idx = hP_level.(task_nm).allTrials(:, iS) == iP_level;
            % fMRI ROI depending on effort level
            fMRI_ROI.(task_nm).hP_level(iP_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jP_idx, iS),1,'omitnan');
            choice_hE.(task_nm).hP_level(iP_level,iS) = mean(choice_hE.(task_nm).allTrials(jP_idx, iS),1,'omitnan');

            % fMRI ROI depending on effort level AND choice (like in Kurniawan
            % et al, 2021)
            jP_highEchosen_idx = (hP_level.(task_nm).allTrials(:, iS) == iP_level).*(choice_hE.(task_nm).allTrials(:, iS) == 1) == 1;
            jP_lowEchosen_idx = (hP_level.(task_nm).allTrials(:, iS) == iP_level).*(choice_hE.(task_nm).allTrials(:, iS) == 0) == 1;
            fMRI_ROI.(task_nm).choice_high.hP_level(iP_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jP_highEchosen_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm).choice_low.hP_level(iP_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jP_lowEchosen_idx, iS),1,'omitnan');
        end % incentive level
        %% incentive chosen
        for iI_ch = 1:n_Ich
            jI_ch_idx = Ichosen.(task_nm).allTrials(:, iS) == (iI_ch - 1);
            fMRI_ROI.(task_nm).Ich(iI_ch,iS) = mean(fMRI_ROI.(task_nm).allTrials(jI_ch_idx, iS),1,'omitnan');
        end % effort chosen
        % reward chosen
        for iR_ch = 1:n_Rch
            jR_ch_idx = Rchosen.(task_nm).allTrials(:, iS) == (iR_ch - 1);
            fMRI_ROI.(task_nm).Ich(iR_ch,iS) = mean(fMRI_ROI.(task_nm).allTrials(jR_ch_idx, iS),1,'omitnan');
        end % effort chosen
        % punishment chosen
        for iP_ch = 1:n_Pch
            jP_ch_idx = Pchosen.(task_nm).allTrials(:, iS) == (iP_ch - 1);
            fMRI_ROI.(task_nm).Pch(iP_ch,iS) = mean(fMRI_ROI.(task_nm).allTrials(jP_ch_idx, iS),1,'omitnan');
        end % effort chosen

        %% test whether slopes are different (similar to what the mediation tests)
        % incentive
        b_fMRI_f_Ich.(task_nm).allTrials.allSubs(:,iS) = glmfit(Ichosen_possible, fMRI_ROI.(task_nm).Ich(:,iS), 'normal');
        b_fMRI_f_I.(task_nm).lE_chosen.allSubs(:,iS) = glmfit(hI_levels, fMRI_ROI.(task_nm).choice_low.hI_level(:,iS), 'normal');
        b_fMRI_f_I.(task_nm).hE_chosen.allSubs(:,iS) = glmfit(hI_levels, fMRI_ROI.(task_nm).choice_high.hI_level(:,iS), 'normal');
        % Reward
        b_fMRI_f_Rch.(task_nm).allTrials.allSubs(:,iS) = glmfit(Rchosen_possible, fMRI_ROI.(task_nm).Rch(:,iS), 'normal');
        b_fMRI_f_R.(task_nm).lE_chosen.allSubs(:,iS) = glmfit(hR_levels, fMRI_ROI.(task_nm).choice_low.hR_level(:,iS), 'normal');
        b_fMRI_f_R.(task_nm).hE_chosen.allSubs(:,iS) = glmfit(hR_levels, fMRI_ROI.(task_nm).choice_high.hR_level(:,iS), 'normal');
        % Punishment
        b_fMRI_f_Pch.(task_nm).allTrials.allSubs(:,iS) = glmfit(Pchosen_possible, fMRI_ROI.(task_nm).Pch(:,iS), 'normal');
        b_fMRI_f_P.(task_nm).lE_chosen.allSubs(:,iS) = glmfit(hP_levels, fMRI_ROI.(task_nm).choice_low.hP_level(:,iS), 'normal');
        b_fMRI_f_P.(task_nm).hE_chosen.allSubs(:,iS) = glmfit(hP_levels, fMRI_ROI.(task_nm).choice_high.hP_level(:,iS), 'normal');
    end % task loop
end % subject loop

%% median split based on metabolites
% extract of subjects based on the median split
[low_met_subs, high_met_subs, metabolite_nm, MRS_ROI_nm] = medSplit_metabolites(study_nm, subject_id);

%% perform extraction for each task independently
for iTask = 1:nTasks
    task_nm = task_names{iTask};
    %% average data independent of metabolites
    [fMRI_ROI.(task_nm).allSubs.choice_low.hI_level.mean,...
        fMRI_ROI.(task_nm).allSubs.choice_low.hI_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hI_level,2);
    [fMRI_ROI.(task_nm).allSubs.choice_high.hI_level.mean,...
        fMRI_ROI.(task_nm).allSubs.choice_high.hI_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hI_level,2);
    [fMRI_ROI.(task_nm).allSubs.choice_low.hR_level.mean,...
        fMRI_ROI.(task_nm).allSubs.choice_low.hR_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hR_level,2);
    [fMRI_ROI.(task_nm).allSubs.choice_high.hR_level.mean,...
        fMRI_ROI.(task_nm).allSubs.choice_high.hR_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hR_level,2);
    [fMRI_ROI.(task_nm).allSubs.choice_low.hP_level.mean,...
        fMRI_ROI.(task_nm).allSubs.choice_low.hP_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hP_level,2);
    [fMRI_ROI.(task_nm).allSubs.choice_high.hP_level.mean,...
        fMRI_ROI.(task_nm).allSubs.choice_high.hP_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hP_level,2);
    
    % perform the median split
    % median split on choices and fMRI for high effort level
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).hI_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).hI_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).hI_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).hI_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).hI_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).hI_level(:,high_met_subs),2);
    [choice_hE.(task_nm).(['l',metabolite_nm]).hI_level.mean,...
        choice_hE.(task_nm).(['l',metabolite_nm]).hI_level.sem] = mean_sem_sd(choice_hE.(task_nm).hI_level(:,low_met_subs),2);
    [choice_hE.(task_nm).(['h',metabolite_nm]).hI_level.mean,...
        choice_hE.(task_nm).(['h',metabolite_nm]).hI_level.sem] = mean_sem_sd(choice_hE.(task_nm).hI_level(:,high_met_subs),2);
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).hR_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).hR_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).hR_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).hR_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).hR_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).hR_level(:,high_met_subs),2);
    [choice_hE.(task_nm).(['l',metabolite_nm]).hR_level.mean,...
        choice_hE.(task_nm).(['l',metabolite_nm]).hR_level.sem] = mean_sem_sd(choice_hE.(task_nm).hR_level(:,low_met_subs),2);
    [choice_hE.(task_nm).(['h',metabolite_nm]).hR_level.mean,...
        choice_hE.(task_nm).(['h',metabolite_nm]).hR_level.sem] = mean_sem_sd(choice_hE.(task_nm).hR_level(:,high_met_subs),2);
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).hP_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).hP_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).hP_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).hP_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).hP_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).hP_level(:,high_met_subs),2);
    [choice_hE.(task_nm).(['l',metabolite_nm]).hP_level.mean,...
        choice_hE.(task_nm).(['l',metabolite_nm]).hP_level.sem] = mean_sem_sd(choice_hE.(task_nm).hP_level(:,low_met_subs),2);
    [choice_hE.(task_nm).(['h',metabolite_nm]).hP_level.mean,...
        choice_hE.(task_nm).(['h',metabolite_nm]).hP_level.sem] = mean_sem_sd(choice_hE.(task_nm).hP_level(:,high_met_subs),2);
    
    % median split on choices and fMRI for effort chosen
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).Ich.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).Ich.sem] = mean_sem_sd(fMRI_ROI.(task_nm).Ich(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).Ich.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).Ich.sem] = mean_sem_sd(fMRI_ROI.(task_nm).Ich(:,high_met_subs),2);
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).Rch.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).Rch.sem] = mean_sem_sd(fMRI_ROI.(task_nm).Rch(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).Rch.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).Rch.sem] = mean_sem_sd(fMRI_ROI.(task_nm).Rch(:,high_met_subs),2);
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).Pch.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).Pch.sem] = mean_sem_sd(fMRI_ROI.(task_nm).Pch(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).Pch.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).Pch.sem] = mean_sem_sd(fMRI_ROI.(task_nm).Pch(:,high_met_subs),2);
    % median split on fMRI for effort level and choice made
    % incentive
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hI_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hI_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hI_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hI_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hI_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hI_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hI_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hI_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hI_level(:,high_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hI_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hI_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hI_level(:,high_met_subs),2);
    % reward
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hR_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hR_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hR_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hR_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hR_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hR_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hR_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hR_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hR_level(:,high_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hR_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hR_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hR_level(:,high_met_subs),2);
    % punishment
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hP_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hP_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hP_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hP_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hP_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hP_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hP_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hP_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hP_level(:,high_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hP_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hP_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hP_level(:,high_met_subs),2);

    % betas
    % global betas with effort chosen
    b_fMRI_f_Ich.(task_nm).(['low_',metabolite_nm]) = b_fMRI_f_Ich.(task_nm).allSubs(:,low_met_subs);
    b_fMRI_f_Ich.(task_nm).(['high_',metabolite_nm]) = b_fMRI_f_Ich.(task_nm).allSubs(:,high_met_subs);
    b_fMRI_f_Rch.(task_nm).(['low_',metabolite_nm]) = b_fMRI_f_Rch.(task_nm).allSubs(:,low_met_subs);
    b_fMRI_f_Rch.(task_nm).(['high_',metabolite_nm]) = b_fMRI_f_Rch.(task_nm).allSubs(:,high_met_subs);
    b_fMRI_f_Pch.(task_nm).(['low_',metabolite_nm]) = b_fMRI_f_Pch.(task_nm).allSubs(:,low_met_subs);
    b_fMRI_f_Pch.(task_nm).(['high_',metabolite_nm]) = b_fMRI_f_Pch.(task_nm).allSubs(:,high_met_subs);
    
    % betas split by effort level and choice
    % incentive
    b_fMRI_f_I.(task_nm).lE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_I.(task_nm).lE_chosen.allSubs(:,low_met_subs);
    b_fMRI_f_I.(task_nm).lE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_I.(task_nm).lE_chosen.allSubs(:,high_met_subs);
    b_fMRI_f_I.(task_nm).hE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_I.(task_nm).hE_chosen.allSubs(:,low_met_subs);
    b_fMRI_f_I.(task_nm).hE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_I.(task_nm).hE_chosen.allSubs(:,high_met_subs);
    % reward
    b_fMRI_f_R.(task_nm).lE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_R.(task_nm).lE_chosen.allSubs(:,low_met_subs);
    b_fMRI_f_R.(task_nm).lE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_R.(task_nm).lE_chosen.allSubs(:,high_met_subs);
    b_fMRI_f_R.(task_nm).hE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_R.(task_nm).hE_chosen.allSubs(:,low_met_subs);
    b_fMRI_f_R.(task_nm).hE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_R.(task_nm).hE_chosen.allSubs(:,high_met_subs);
    % punishment
    b_fMRI_f_P.(task_nm).lE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_P.(task_nm).lE_chosen.allSubs(:,low_met_subs);
    b_fMRI_f_P.(task_nm).lE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_P.(task_nm).lE_chosen.allSubs(:,high_met_subs);
    b_fMRI_f_P.(task_nm).hE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_P.(task_nm).hE_chosen.allSubs(:,low_met_subs);
    b_fMRI_f_P.(task_nm).hE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_P.(task_nm).hE_chosen.allSubs(:,high_met_subs);
    
    %% ttest
    % incentive
    [~,pval.(task_nm).choice_hE_vs_lE.allSubs] = ttest(b_fMRI_f_I.(task_nm).lE_chosen.allSubs(2,:),...
        b_fMRI_f_I.(task_nm).hE_chosen.allSubs(2,:));
    [~,pval.(task_nm).choice_hE_vs_choice_hE_slope.(['low_vs_high_',metabolite_nm])] = ttest2(b_fMRI_f_I.(task_nm).hE_chosen.(['low_',metabolite_nm])(2,:),...
        b_fMRI_f_I.(task_nm).hE_chosen.(['high_',metabolite_nm])(2,:));
    % reward
    [~,pval.(task_nm).choice_hE_vs_lE.allSubs] = ttest(b_fMRI_f_R.(task_nm).lE_chosen.allSubs(2,:),...
        b_fMRI_f_R.(task_nm).hE_chosen.allSubs(2,:));
    [~,pval.(task_nm).choice_hE_vs_choice_hE_slope.(['low_vs_high_',metabolite_nm])] = ttest2(b_fMRI_f_R.(task_nm).hE_chosen.(['low_',metabolite_nm])(2,:),...
        b_fMRI_f_R.(task_nm).hE_chosen.(['high_',metabolite_nm])(2,:));
    % punishment
    [~,pval.(task_nm).choice_hE_vs_lE.allSubs] = ttest(b_fMRI_f_P.(task_nm).lE_chosen.allSubs(2,:),...
        b_fMRI_f_P.(task_nm).hE_chosen.allSubs(2,:));
    [~,pval.(task_nm).choice_hE_vs_choice_hE_slope.(['low_vs_high_',metabolite_nm])] = ttest2(b_fMRI_f_P.(task_nm).hE_chosen.(['low_',metabolite_nm])(2,:),...
        b_fMRI_f_P.(task_nm).hE_chosen.(['high_',metabolite_nm])(2,:));
    
    %% figure
    if dispFig == true
        % general parameters
        lWidth = 3;
        black = [0 0 0];
        purple = [148 0 211]./255;
        orange = [241 163 64]./255;
        blue = [51 153 255]./255;
        pSize = 50;

        %% fMRI = f(E chosen) and metabolites
        fig;
        % incentive
        subplot(1,3,1);
        hold on;
        fMRI_Ich_low = errorbar(Rchosen_possible-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).Rch.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).Rch.sem);
        fMRI_Ich_high = errorbar(Rchosen_possible+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).Rch.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).Rch.sem);
        % curve parameters
        fMRI_Ich_low.LineStyle = '--';
        fMRI_Ich_low.LineWidth = lWidth;
        fMRI_Ich_low.MarkerEdgeColor = blue;
        fMRI_Ich_high.LineStyle = '-';
        fMRI_Ich_high.LineWidth = lWidth;
        fMRI_Ich_high.MarkerEdgeColor = purple;
        legend([fMRI_Ich_high, fMRI_Ich_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(Rchosen_possible);
        xlabel('Incentive chosen');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);
        
        % reward
        subplot(1,3,2);
        hold on;
        fMRI_Rch_low = errorbar(Rchosen_possible-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).Rch.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).Rch.sem);
        fMRI_Rch_high = errorbar(Rchosen_possible+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).Rch.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).Rch.sem);
        % curve parameters
        fMRI_Rch_low.LineStyle = '--';
        fMRI_Rch_low.LineWidth = lWidth;
        fMRI_Rch_low.MarkerEdgeColor = blue;
        fMRI_Rch_high.LineStyle = '-';
        fMRI_Rch_high.LineWidth = lWidth;
        fMRI_Rch_high.MarkerEdgeColor = purple;
        legend([fMRI_Rch_high, fMRI_Rch_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(Rchosen_possible);
        xlabel('R chosen');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);
        
        % punishment
        subplot(1,3,3);
        hold on;
        fMRI_Pch_low = errorbar(Rchosen_possible-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).Rch.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).Rch.sem);
        fMRI_Pch_high = errorbar(Rchosen_possible+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).Rch.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).Rch.sem);
        % curve parameters
        fMRI_Pch_low.LineStyle = '--';
        fMRI_Pch_low.LineWidth = lWidth;
        fMRI_Pch_low.MarkerEdgeColor = blue;
        fMRI_Pch_high.LineStyle = '-';
        fMRI_Pch_high.LineWidth = lWidth;
        fMRI_Pch_high.MarkerEdgeColor = purple;
        legend([fMRI_Pch_high, fMRI_Pch_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(Pchosen_possible);
        xlabel('P chosen');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);
        
        %% fMRI and choice = f(high E proposed) and metabolites
        fig;
        % fMRI
        subplot(3,2,1);
        hold on;
        fMRI_hI_low = errorbar(hI_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).hI_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).hI_level.sem);
        fMRI_hI_high = errorbar(hI_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).hI_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).hI_level.sem);
        % curve parameters
        fMRI_hI_low.LineStyle = '--';
        fMRI_hI_low.LineWidth = lWidth;
        fMRI_hI_low.MarkerEdgeColor = blue;
        fMRI_hI_high.LineStyle = '-';
        fMRI_hI_high.LineWidth = lWidth;
        fMRI_hI_high.MarkerEdgeColor = purple;
        legend([fMRI_hI_high, fMRI_hI_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xlabel('high Incentive level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);

        % choice
        subplot(3,2,2);
        hold on;
        choice_hE_low = errorbar(hI_levels-0.01,...
            choice_hE.(task_nm).(['l',metabolite_nm]).hI_level.mean,...
            choice_hE.(task_nm).(['l',metabolite_nm]).hI_level.sem);
        choice_hE_high = errorbar(hI_levels+0.01,...
            choice_hE.(task_nm).(['h',metabolite_nm]).hI_level.mean,...
            choice_hE.(task_nm).(['h',metabolite_nm]).hI_level.sem);
        % curve parameters
        choice_hE_low.LineStyle = '--';
        choice_hE_low.LineWidth = lWidth;
        choice_hE_low.MarkerEdgeColor = blue;
        choice_hE_high.LineStyle = '-';
        choice_hE_high.LineWidth = lWidth;
        choice_hE_high.MarkerEdgeColor = purple;
        legend([choice_hE_high, choice_hE_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xlabel('high Incentive level');
        ylabel(['Choice = high effort ',task_nm,' (%)']);
        legend_size(pSize);
        
        
        % fMRI
        subplot(3,2,3);
        hold on;
        fMRI_hI_low = errorbar(hI_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).hR_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).hR_level.sem);
        fMRI_hI_high = errorbar(hI_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).hR_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).hR_level.sem);
        % curve parameters
        fMRI_hI_low.LineStyle = '--';
        fMRI_hI_low.LineWidth = lWidth;
        fMRI_hI_low.MarkerEdgeColor = blue;
        fMRI_hI_high.LineStyle = '-';
        fMRI_hI_high.LineWidth = lWidth;
        fMRI_hI_high.MarkerEdgeColor = purple;
        legend([fMRI_hI_high, fMRI_hI_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xlabel('high Reward level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);

        % choice
        subplot(3,2,4);
        hold on;
        choice_hE_low = errorbar(hI_levels-0.01,...
            choice_hE.(task_nm).(['l',metabolite_nm]).hR_level.mean,...
            choice_hE.(task_nm).(['l',metabolite_nm]).hR_level.sem);
        choice_hE_high = errorbar(hI_levels+0.01,...
            choice_hE.(task_nm).(['h',metabolite_nm]).hR_level.mean,...
            choice_hE.(task_nm).(['h',metabolite_nm]).hR_level.sem);
        % curve parameters
        choice_hE_low.LineStyle = '--';
        choice_hE_low.LineWidth = lWidth;
        choice_hE_low.MarkerEdgeColor = blue;
        choice_hE_high.LineStyle = '-';
        choice_hE_high.LineWidth = lWidth;
        choice_hE_high.MarkerEdgeColor = purple;
        legend([choice_hE_high, choice_hE_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xlabel('high Reward level');
        ylabel(['Choice = high effort ',task_nm,' (%)']);
        legend_size(pSize);
        
        
        % fMRI
        subplot(3,2,5);
        hold on;
        fMRI_hI_low = errorbar(hR_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).hP_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).hP_level.sem);
        fMRI_hI_high = errorbar(hR_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).hP_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).hP_level.sem);
        % curve parameters
        fMRI_hI_low.LineStyle = '--';
        fMRI_hI_low.LineWidth = lWidth;
        fMRI_hI_low.MarkerEdgeColor = blue;
        fMRI_hI_high.LineStyle = '-';
        fMRI_hI_high.LineWidth = lWidth;
        fMRI_hI_high.MarkerEdgeColor = purple;
        legend([fMRI_hI_high, fMRI_hI_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xlabel('high Reward level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);

        % choice
        subplot(3,2,6);
        hold on;
        choice_hE_low = errorbar(hP_levels-0.01,...
            choice_hE.(task_nm).(['l',metabolite_nm]).hP_level.mean,...
            choice_hE.(task_nm).(['l',metabolite_nm]).hP_level.sem);
        choice_hE_high = errorbar(hP_levels+0.01,...
            choice_hE.(task_nm).(['h',metabolite_nm]).hP_level.mean,...
            choice_hE.(task_nm).(['h',metabolite_nm]).hP_level.sem);
        % curve parameters
        choice_hE_low.LineStyle = '--';
        choice_hE_low.LineWidth = lWidth;
        choice_hE_low.MarkerEdgeColor = blue;
        choice_hE_high.LineStyle = '-';
        choice_hE_high.LineWidth = lWidth;
        choice_hE_high.MarkerEdgeColor = purple;
        legend([choice_hE_high, choice_hE_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xlabel('high Reward level');
        ylabel(['Choice = high effort ',task_nm,' (%)']);
        legend_size(pSize);

        %% fMRI = f(E level) according to choice made independent of metabolites
        fig;
        
        % incentive
        subplot(1,3,1);
        hold on;
        fMRI_lowEch = errorbar(hI_levels-0.01,...
            fMRI_ROI.(task_nm).allSubs.choice_low.hI_level.mean,...
            fMRI_ROI.(task_nm).allSubs.choice_low.hI_level.sem);
        fMRI_highEch = errorbar(hI_levels+0.01,...
            fMRI_ROI.(task_nm).allSubs.choice_high.hI_level.mean,...
            fMRI_ROI.(task_nm).allSubs.choice_high.hI_level.sem);
        % curve parameters
        fMRI_lowEch.LineStyle = '--';
        fMRI_lowEch.LineWidth = lWidth;
        fMRI_lowEch.Color = blue;
        fMRI_highEch.LineStyle = '-';
        fMRI_highEch.LineWidth = lWidth;
        fMRI_highEch.Color = purple;
        legend([fMRI_highEch, fMRI_lowEch],...
            {'high E chosen','low E chosen'});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(hI_levels);
        xlabel('high Incentive level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);
        
        % reward
        subplot(1,3,2);
        hold on;
        fMRI_lowEch = errorbar(hR_levels-0.01,...
            fMRI_ROI.(task_nm).allSubs.choice_low.hR_level.mean,...
            fMRI_ROI.(task_nm).allSubs.choice_low.hR_level.sem);
        fMRI_highEch = errorbar(hI_levels+0.01,...
            fMRI_ROI.(task_nm).allSubs.choice_high.hR_level.mean,...
            fMRI_ROI.(task_nm).allSubs.choice_high.hR_level.sem);
        % curve parameters
        fMRI_lowEch.LineStyle = '--';
        fMRI_lowEch.LineWidth = lWidth;
        fMRI_lowEch.Color = blue;
        fMRI_highEch.LineStyle = '-';
        fMRI_highEch.LineWidth = lWidth;
        fMRI_highEch.Color = purple;
        legend([fMRI_highEch, fMRI_lowEch],...
            {'high E chosen','low E chosen'});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(hR_levels);
        xlabel('high Reward level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);
        
        % punishment
        subplot(1,3,3);
        hold on;
        fMRI_lowEch = errorbar(hP_levels-0.01,...
            fMRI_ROI.(task_nm).allSubs.choice_low.hP_level.mean,...
            fMRI_ROI.(task_nm).allSubs.choice_low.hP_level.sem);
        fMRI_highEch = errorbar(hP_levels+0.01,...
            fMRI_ROI.(task_nm).allSubs.choice_high.hP_level.mean,...
            fMRI_ROI.(task_nm).allSubs.choice_high.hP_level.sem);
        % curve parameters
        fMRI_lowEch.LineStyle = '--';
        fMRI_lowEch.LineWidth = lWidth;
        fMRI_lowEch.Color = blue;
        fMRI_highEch.LineStyle = '-';
        fMRI_highEch.LineWidth = lWidth;
        fMRI_highEch.Color = purple;
        legend([fMRI_highEch, fMRI_lowEch],...
            {'high E chosen','low E chosen'});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(hP_levels);
        xlabel('high Punishment level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);

        %% fMRI = f(E level) according to choice made
        fig;
        
        % incentive
        subplot(1,3,1);
        hold on;
        fMRI_lowEch_lowMet = errorbar(hI_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hI_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hI_level.sem);
        fMRI_lowEch_highMet = errorbar(hI_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hI_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hI_level.sem);
        fMRI_highEch_lowMet = errorbar(hI_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hI_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hI_level.sem);
        fMRI_highEch_highMet = errorbar(hI_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hI_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hI_level.sem);
        % curve parameters
        fMRI_lowEch_lowMet.LineStyle = '--';
        fMRI_lowEch_lowMet.LineWidth = lWidth;
        fMRI_lowEch_lowMet.Color = blue;
        fMRI_lowEch_highMet.LineStyle = '-';
        fMRI_lowEch_highMet.LineWidth = lWidth;
        fMRI_lowEch_highMet.Color = blue;
        fMRI_highEch_lowMet.LineStyle = '--';
        fMRI_highEch_lowMet.LineWidth = lWidth;
        fMRI_highEch_lowMet.Color = purple;
        fMRI_highEch_highMet.LineStyle = '-';
        fMRI_highEch_highMet.LineWidth = lWidth;
        fMRI_highEch_highMet.Color = purple;
        legend([fMRI_highEch_highMet, fMRI_highEch_lowMet,...
            fMRI_lowEch_highMet, fMRI_lowEch_lowMet],...
            {['high ',metabolite_nm, '- high E chosen'],['low ', metabolite_nm, '- high E chosen'],...
            ['high ',metabolite_nm, '- low E chosen'],['low ', metabolite_nm, '- low E chosen']});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(hI_levels);
        xlabel('high Incentive level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);
        
        % reward
        subplot(1,3,2);
        hold on;
        fMRI_lowEch_lowMet = errorbar(hR_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hR_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hR_level.sem);
        fMRI_lowEch_highMet = errorbar(hR_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hR_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hR_level.sem);
        fMRI_highEch_lowMet = errorbar(hR_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hR_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hR_level.sem);
        fMRI_highEch_highMet = errorbar(hR_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hR_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hR_level.sem);
        % curve parameters
        fMRI_lowEch_lowMet.LineStyle = '--';
        fMRI_lowEch_lowMet.LineWidth = lWidth;
        fMRI_lowEch_lowMet.Color = blue;
        fMRI_lowEch_highMet.LineStyle = '-';
        fMRI_lowEch_highMet.LineWidth = lWidth;
        fMRI_lowEch_highMet.Color = blue;
        fMRI_highEch_lowMet.LineStyle = '--';
        fMRI_highEch_lowMet.LineWidth = lWidth;
        fMRI_highEch_lowMet.Color = purple;
        fMRI_highEch_highMet.LineStyle = '-';
        fMRI_highEch_highMet.LineWidth = lWidth;
        fMRI_highEch_highMet.Color = purple;
        legend([fMRI_highEch_highMet, fMRI_highEch_lowMet,...
            fMRI_lowEch_highMet, fMRI_lowEch_lowMet],...
            {['high ',metabolite_nm, '- high E chosen'],['low ', metabolite_nm, '- high E chosen'],...
            ['high ',metabolite_nm, '- low E chosen'],['low ', metabolite_nm, '- low E chosen']});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(hR_levels);
        xlabel('high Reward level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);
        
        % incentive
        subplot(1,3,1);
        hold on;
        fMRI_lowEch_lowMet = errorbar(hP_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hI_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hI_level.sem);
        fMRI_lowEch_highMet = errorbar(hP_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hI_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hI_level.sem);
        fMRI_highEch_lowMet = errorbar(hP_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hI_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hI_level.sem);
        fMRI_highEch_highMet = errorbar(hP_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hI_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hI_level.sem);
        % curve parameters
        fMRI_lowEch_lowMet.LineStyle = '--';
        fMRI_lowEch_lowMet.LineWidth = lWidth;
        fMRI_lowEch_lowMet.Color = blue;
        fMRI_lowEch_highMet.LineStyle = '-';
        fMRI_lowEch_highMet.LineWidth = lWidth;
        fMRI_lowEch_highMet.Color = blue;
        fMRI_highEch_lowMet.LineStyle = '--';
        fMRI_highEch_lowMet.LineWidth = lWidth;
        fMRI_highEch_lowMet.Color = purple;
        fMRI_highEch_highMet.LineStyle = '-';
        fMRI_highEch_highMet.LineWidth = lWidth;
        fMRI_highEch_highMet.Color = purple;
        legend([fMRI_highEch_highMet, fMRI_highEch_lowMet,...
            fMRI_lowEch_highMet, fMRI_lowEch_lowMet],...
            {['high ',metabolite_nm, '- high E chosen'],['low ', metabolite_nm, '- high E chosen'],...
            ['high ',metabolite_nm, '- low E chosen'],['low ', metabolite_nm, '- low E chosen']});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(hP_levels);
        xlabel('high Punishment level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);
    end % figure display
end % task loop