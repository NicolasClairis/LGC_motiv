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
    task_to_look, timePeriod_nm] = extract_ROI_betas_onsets_only_bis(computerRoot,...
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
    [fMRI_ROI.(task_nm).hE_level, choice_hE.(task_nm).hE_level] = deal(NaN(n_hI_levels, NS));
    
    % data split by incentive*choice
    [fMRI_ROI.(task_nm).lEch.hR_level, fMRI_ROI.(task_nm).hEch.hR_level] = deal(NaN(n_hR_levels, NS));
    [fMRI_ROI.(task_nm).lEch.hP_level, fMRI_ROI.(task_nm).hEch.hP_level] = deal(NaN(n_hP_levels, NS));
    % data split by effort*choice
    [fMRI_ROI.(task_nm).lEch.hE_level, fMRI_ROI.(task_nm).hEch.hE_level] = deal(NaN(n_hP_levels, NS));
end % task loop

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];

    % extract runs
    [runsStruct, n_runs] = runs_definition(study_nm, sub_nm, 'behavior');
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    for iRun = 1:n_runs
        kRun = okRuns(iRun);
        run_nm = num2str(kRun);
        task_nm_tmp = taskNames{jRun};
        switch kRun
            case {1,2}
                jRun = 1;
            case {3,4}
                jRun = 2;
        end
        run_nm_bis = num2str(jRun);
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
        fMRI_ROI.(task_nm_tmp).allTrials(runTrials_idx, iS) = ROI_trial_b_trial.(fMRI_ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
    end % run loop
    
    %% now average the data according to R/P/E and choice
    for iT = 1:nTasks
        task_nm = tasks{iT};
        
        % extract data per reward level
        for iR = 1:n_hR_levels
            
        end % reward
    end % task loop
end % subject loop