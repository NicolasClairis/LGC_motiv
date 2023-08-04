%% extract the level of activity of a given ROI (extracted trial by trial 
% with extract_ROI_betas_onsets_only_bis.m) in function of the varying
% reward, punishment and the level of effort chosen

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% working directory
computerRoot = LGCM_root_paths;
studyFolder = [computerRoot, study_nm, filesep];

%% load ROI
[ROI_trial_b_trial, ROI_subList,...
    ROI_nm, ROI_short_nm,...
    ~, timePeriod_nm] = extract_ROI_betas_onsets_only_bis(computerRoot,...
    study_nm, subject_id, condition);

%% main parameters of interest
nRunsPerTask = 2;
nTrialsPerRun = 54;
nTrialsPerTask = nTrialsPerRun*nRunsPerTask;

hR_levels = 1:3;
n_hR_levels = length(hR_levels);
hP_levels = 1:3;
n_hP_levels = length(hP_levels);
hE_levels = 1:3;
n_hE_levels = length(hE_levels);

tasks = {'Ep','Em'};
nTasks = length(tasks);

for iT = 1:nTasks
    task_nm = tasks{iT};
    % data for all trials
    [fMRI_ROI.(task_nm).allTrials, choice_hE.(task_nm).allTrials,...
        RP.(task_nm).allTrials,...
        hR_level.(task_nm).allTrials, Rchosen.(task_nm).allTrials,...
        hP_level.(task_nm).allTrials, Pchosen.(task_nm).allTrials,...
        hE_level.(task_nm).allTrials, Echosen.(task_nm).allTrials] = deal(NaN(nTrialsPerTask, NS));
    % data split by level of incentive proposed for high effort option
    [fMRI_ROI.(task_nm).hR_level, choice_hE.(task_nm).hR_level] = deal(NaN(n_hR_levels, NS));
    [fMRI_ROI.(task_nm).hP_level, choice_hE.(task_nm).hP_level] = deal(NaN(n_hP_levels, NS));
    
    % data split by level of effort proposed for high effort option
    [fMRI_ROI.(task_nm).hE_level, choice_hE.(task_nm).hE_level] = deal(NaN(n_hE_levels, NS));
    
    % data split by incentive*choice
    [fMRI_ROI.(task_nm).lEch.hR_level, fMRI_ROI.(task_nm).hEch.hR_level] = deal(NaN(n_hR_levels, NS));
    [fMRI_ROI.(task_nm).lEch.hP_level, fMRI_ROI.(task_nm).hEch.hP_level] = deal(NaN(n_hP_levels, NS));
    % data split by effort*choice
    [fMRI_ROI.(task_nm).lEch.hE_level, fMRI_ROI.(task_nm).hEch.hE_level] = deal(NaN(n_hP_levels, NS));
end % task loop

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyFolder, 'CID',sub_nm, filesep, 'behavior',filesep];

    % extract runs
    [runsStruct, n_runs] = runs_definition(study_nm, sub_nm, 'behavior');
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    for iRun = 1:n_runs
        kRun = okRuns(iRun);
        run_nm = num2str(kRun);
        task_nm_tmp = taskNames{kRun};
        switch kRun
            case {1,2}
                jRun = 1;
            case {3,4}
                jRun = 2;
        end
        run_nm_bis = ['run',num2str(jRun)];
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun - 1);
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        
        %% extract relevant variables
        % R/P trials
        RP_var_tmp = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        R_trials = RP_var_tmp == 1;
        P_trials = RP_var_tmp == -1;
        RP.(task_nm_tmp).allTrials(runTrials_idx,iS) = RP_var_tmp;
        % choice
        choice_hE.(task_nm_tmp).allTrials(runTrials_idx,iS) = extract_choice_hE(subBehaviorFolder,sub_nm,run_nm,task_fullName);
        % extract R, P and E level for high effort option and for choice
        hR_level.(task_nm_tmp).allTrials(runTrials_idx,iS) = extract_hR_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        hP_level.(task_nm_tmp).allTrials(runTrials_idx,iS) = extract_hR_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        hE_level.(task_nm_tmp).allTrials(runTrials_idx,iS) = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        Rchosen.(task_nm_tmp).allTrials(runTrials_idx,iS) = extract_R_level_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        Pchosen.(task_nm_tmp).allTrials(runTrials_idx,iS) = extract_P_level_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        Echosen.(task_nm_tmp).allTrials(runTrials_idx,iS) = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        
        %% extract fMRI ROI mediator
        fMRI_ROI.(task_nm_tmp).allTrials(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
    end % run loop
    
    %% now average the data according to R/P/E and choice
    for iT = 1:nTasks
        task_nm = tasks{iT};
        
        %% extract data per high effort option features
        % extract data per reward level
        for iR = 1:n_hR_levels
            R_idx = hR_level.(task_nm).allTrials(:,iS) == iR;
            fMRI_ROI.(task_nm).hR_level(iR,iS) = mean(fMRI_ROI.(task_nm_tmp).allTrials(R_idx, iS),1,'omitnan');
            choice_hE.(task_nm).hR_level(iR,iS) = mean(choice_hE.(task_nm_tmp).allTrials(R_idx, iS),1,'omitnan');
        end % reward

        % extract data per punishment level
        for iP = 1:n_hP_levels
            P_idx = hP_level.(task_nm).allTrials(:,iS) == iP;
            fMRI_ROI.(task_nm).hP_level(iP,iS) = mean(fMRI_ROI.(task_nm_tmp).allTrials(P_idx, iS),1,'omitnan');
            choice_hE.(task_nm).hP_level(iP,iS) = mean(choice_hE.(task_nm_tmp).allTrials(P_idx, iS),1,'omitnan');
        end % punishment

        % extract data per effort level
        for iE = 1:n_hE_levels
            E_idx = hE_level.(task_nm).allTrials(:,iS) == iE;
            fMRI_ROI.(task_nm).hE_level(iE,iS) = mean(fMRI_ROI.(task_nm_tmp).allTrials(E_idx, iS),1,'omitnan');
            choice_hE.(task_nm).hE_level(iE,iS) = mean(choice_hE.(task_nm_tmp).allTrials(E_idx, iS),1,'omitnan');
        end % effort

        %% extract data depending on choice (high/low) and high effort option features
        choice_highE_idx_tmp    = choice_hE.(task_nm_tmp).allTrials(:, iS) == 1;

        % reward
        for iR = 1:n_hR_levels
            R_hE_idx = (hR_level.(task_nm).allTrials(:,iS) == iR).*(choice_highE_idx_tmp == 1) == 1;
            R_lE_idx = (hR_level.(task_nm).allTrials(:,iS) == iR).*(choice_highE_idx_tmp == 0) == 1;
            fMRI_ROI.(task_nm).hEch.hR_level(iR,iS) = mean(fMRI_ROI.(task_nm_tmp).allTrials(R_hE_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm).lEch.hR_level(iR,iS) = mean(fMRI_ROI.(task_nm_tmp).allTrials(R_lE_idx, iS),1,'omitnan');
        end

        % punishment
        for iP = 1:n_hP_levels
            P_hE_idx = (hP_level.(task_nm).allTrials(:,iS) == iP).*(choice_highE_idx_tmp == 1) == 1;
            P_lE_idx = (hP_level.(task_nm).allTrials(:,iS) == iP).*(choice_highE_idx_tmp == 0) == 1;
            fMRI_ROI.(task_nm).hEch.hP_level(iP,iS) = mean(fMRI_ROI.(task_nm_tmp).allTrials(P_hE_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm).lEch.hP_level(iP,iS) = mean(fMRI_ROI.(task_nm_tmp).allTrials(P_lE_idx, iS),1,'omitnan');
        end

        % effort
        for iE = 1:n_hE_levels
            E_hE_idx = (hP_level.(task_nm).allTrials(:,iS) == iE).*(choice_highE_idx_tmp == 1) == 1;
            E_lE_idx = (hP_level.(task_nm).allTrials(:,iS) == iE).*(choice_highE_idx_tmp == 0) == 1;
            fMRI_ROI.(task_nm).hEch.hE_level(iE,iS) = mean(fMRI_ROI.(task_nm_tmp).allTrials(E_hE_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm).lEch.hE_level(iE,iS) = mean(fMRI_ROI.(task_nm_tmp).allTrials(E_lE_idx, iS),1,'omitnan');
        end
    end % task loop
end % subject loop

%% figure parameter
[pSize, lWidth, col, mSize] = general_fig_prm;
figHighE = fig;
figChoice_x_highE = fig;

%% loop through tasks for average + figure
for iT = 1:nTasks
    task_nm = tasks{iT};

    %% average across subjects
    % data split by level of incentive proposed for high effort option
    [mean_fMRI_ROI.(task_nm).hR_level,...
        sem_fMRI_ROI.(task_nm).hR_level] = mean_sem_sd(fMRI_ROI.(task_nm).hR_level, 2);
    [mean_choice_hE.(task_nm).hR_level,...
        sem_choice_hE.(task_nm).hR_level] = mean_sem_sd(choice_hE.(task_nm).hR_level,2);
    [mean_fMRI_ROI.(task_nm).hP_level,...
        sem_fMRI_ROI.(task_nm).hP_level] = mean_sem_sd(fMRI_ROI.(task_nm).hP_level,2);
    [mean_choice_hE.(task_nm).hP_level,...
        sem_choice_hE.(task_nm).hP_level] = mean_sem_sd(choice_hE.(task_nm).hP_level,2);
    
    % data split by level of effort proposed for high effort option
    [mean_fMRI_ROI.(task_nm).hE_level,...
        sem_fMRI_ROI.(task_nm).hE_level] = mean_sem_sd(fMRI_ROI.(task_nm).hE_level,2);
    [mean_choice_hE.(task_nm).hE_level,...
        sem_choice_hE.(task_nm).hE_level] = mean_sem_sd(choice_hE.(task_nm).hE_level,2);
    
    % data split by incentive*choice
    [mean_fMRI_ROI.(task_nm).lEch.hR_level,...
        sem_fMRI_ROI.(task_nm).lEch.hR_level] = mean_sem_sd(fMRI_ROI.(task_nm).lEch.hR_level,2);
    [mean_fMRI_ROI.(task_nm).hEch.hR_level,...
        sem_fMRI_ROI.(task_nm).hEch.hR_level] = mean_sem_sd(fMRI_ROI.(task_nm).hEch.hR_level,2);
    [mean_fMRI_ROI.(task_nm).lEch.hP_level,...
        sem_fMRI_ROI.(task_nm).lEch.hP_level] = mean_sem_sd(fMRI_ROI.(task_nm).lEch.hP_level,2);
    [mean_fMRI_ROI.(task_nm).hEch.hP_level,...
        sem_fMRI_ROI.(task_nm).hEch.hP_level] = mean_sem_sd(fMRI_ROI.(task_nm).hEch.hP_level,2);

    % data split by effort*choice
    [mean_fMRI_ROI.(task_nm).lEch.hE_level,...
        sem_fMRI_ROI.(task_nm).lEch.hE_level] = mean_sem_sd(fMRI_ROI.(task_nm).lEch.hE_level,2);
    [mean_fMRI_ROI.(task_nm).hEch.hE_level,...
        sem_fMRI_ROI.(task_nm).hEch.hE_level] = mean_sem_sd(fMRI_ROI.(task_nm).hEch.hE_level,2);

    %% figures

    %% figure focused on varying option
    figure(figHighE);

    % reward
    subplot(nTasks,3,1 + 3*(iT-1)); hold on;
    R_er_hdl = errorbar(1:n_hR_levels,...
        mean_fMRI_ROI.(task_nm).hR_level,...
        sem_fMRI_ROI.(task_nm).hR_level);
    errorbar_hdl_upgrade(R_er_hdl);
    xlabel('Reward');
    ylabel(ROI_short_nm);
    legend_size(pSize);

    % punishment
    subplot(nTasks,3,2 + 3*(iT-1)); hold on;
    P_er_hdl = errorbar(1:n_hP_levels,...
        mean_fMRI_ROI.(task_nm).hP_level,...
        sem_fMRI_ROI.(task_nm).hP_level);
    errorbar_hdl_upgrade(P_er_hdl);
    xlabel('Punishment');
    ylabel(ROI_short_nm);
    legend_size(pSize);

    % effort
    subplot(nTasks,3,3 + 3*(iT-1)); hold on;
    E_er_hdl = errorbar(1:n_hE_levels,...
        mean_fMRI_ROI.(task_nm).hE_level,...
        sem_fMRI_ROI.(task_nm).hE_level);
    errorbar_hdl_upgrade(E_er_hdl);
    xlabel('Effort');
    ylabel(ROI_short_nm);
    legend_size(pSize);


    %% figure for choice*high effort option
    figure(figChoice_x_highE);

    % reward
    subplot(nTasks,3,1 + 3*(iT-1)); hold on;
    R_hE_er_hdl = errorbar(1:n_hR_levels,...
        mean_fMRI_ROI.(task_nm).hEch.hR_level,...
        sem_fMRI_ROI.(task_nm).hEch.hR_level);
    R_lE_er_hdl = errorbar(1:n_hR_levels,...
        mean_fMRI_ROI.(task_nm).lEch.hR_level,...
        sem_fMRI_ROI.(task_nm).lEch.hR_level);
    errorbar_hdl_upgrade(R_hE_er_hdl, col.orange);
    errorbar_hdl_upgrade(R_lE_er_hdl, col.blue);
    xlabel('Reward');
    ylabel(ROI_short_nm);
    legend_size(pSize);

    % punishment
    subplot(nTasks,3,2 + 3*(iT-1)); hold on;
    P_hE_er_hdl = errorbar(1:n_hP_levels,...
        mean_fMRI_ROI.(task_nm).hEch.hP_level,...
        sem_fMRI_ROI.(task_nm).hEch.hP_level);
    P_lE_er_hdl = errorbar(1:n_hP_levels,...
        mean_fMRI_ROI.(task_nm).lEch.hP_level,...
        sem_fMRI_ROI.(task_nm).lEch.hP_level);
    errorbar_hdl_upgrade(P_hE_er_hdl, col.orange);
    errorbar_hdl_upgrade(P_lE_er_hdl, col.blue);
    xlabel('Punishment');
    ylabel(ROI_short_nm);
    legend_size(pSize);

    % effort
    subplot(nTasks,3,3 + 3*(iT-1)); hold on;
    E_hE_er_hdl = errorbar(1:n_hE_levels,...
        mean_fMRI_ROI.(task_nm).hEch.hE_level,...
        sem_fMRI_ROI.(task_nm).hEch.hE_level);
    E_lE_er_hdl = errorbar(1:n_hE_levels,...
        mean_fMRI_ROI.(task_nm).lEch.hE_level,...
        sem_fMRI_ROI.(task_nm).lEch.hE_level);
    errorbar_hdl_upgrade(E_hE_er_hdl, col.orange);
    errorbar_hdl_upgrade(E_lE_er_hdl, col.blue);
    xlabel('Effort');
    ylabel(ROI_short_nm);
    legend_size(pSize);
end % task loop