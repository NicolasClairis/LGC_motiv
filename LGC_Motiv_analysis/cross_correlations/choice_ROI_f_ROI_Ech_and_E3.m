% script aiming at looking at choices=f(E) depending on the ROI slope for
% effort chosen or the delta between Effort level 3 chosen and non-chosen
% (ie BOLD threshold for decision).


%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% working directories
computerRoot = LGCM_root_paths;
studyFolder = [computerRoot, filesep, study_nm, filesep];

%% main parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
% which ROI contrast do you want to base the analysis on?
contrasts = {'Ech_slope','E3_delta'};
which_con = listdlg('PromptString','which contrast?',....
    'ListString',contrasts);
con_to_use = contrasts{which_con};

% orthogonalize ROI to RT or not?
if ~exist('RT_orth','var') || isempty(RT_orth)
    RT_orth = 0;
end

% display figure?
figDisp = 1;

% load ROI for the corresponding contrast
switch con_to_use
    case 'Ech_slope'
        [con_vec_all,...
            ~, ~, ~,...
            con_names,...
            ROI_coords, ~] = ROI_extraction_group(study_nm, [],...
            subject_id, condition, 0);
        % extract slope for each task
        Ep_Ech = strcmp(con_names,'Ep REG choice RP E: effort chosen');
        Em_Ech = strcmp(con_names,'Em REG choice RP E: effort chosen');
        if sum(Ep_Ech) ~= 1 || sum(Em_Ech) ~= 1
            error('problem in identifying Ech slope contrast in the selected GLM');
        end
        
        %% extract corresponding values for each subject
        Ech_slope.Ep = con_vec_all(Ep_Ech, :);
        Ech_slope.Em = con_vec_all(Em_Ech, :);
        
        %% perform median split
        % extract median
        med_Ech_slope.Ep = median(Ech_slope.Ep,2,'omitnan');
        med_Ech_slope.Em = median(Ech_slope.Em,2,'omitnan');
        
        % split between low and high
        % physical
        low_con.Ep = Ech_slope.Ep <= med_Ech_slope.Ep;
        high_con.Ep = Ech_slope.Ep > med_Ech_slope.Ep;
        % mental
        low_con.Em = Ech_slope.Em <= med_Ech_slope.Em;
        high_con.Em = Ech_slope.Em > med_Ech_slope.Em;
        
        %% extract corresponding activity trial/trial
        [ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
            study_nm, subject_id, condition);
        find_ROI_name = fieldnames(ROI_trial_b_trial);
        fMRI_ROI_trialPerTrial_idx = ~strcmp(find_ROI_name,'subject_id');
        fMRI_ROI_trialPerTrial_name = find_ROI_name{fMRI_ROI_trialPerTrial_idx};
        timePeriod_nm = 'choice';
    case 'E3_delta'
        [ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
            study_nm, subject_id, condition);
        timePeriod_nm = 'choice';
        
        % extract delta for each task (remove subjects where data is
        % missing)
        
end

nTrialsPerRun = 54;
nRunsPerTask = 2;
nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
[fMRI_allTrials.Ep, choice_hE_allTrials.Ep, E_level.Ep, choice_hE_fit_allTrials.Ep, RT_allTrials.Ep,...
    fMRI_allTrials.Em, choice_hE_allTrials.Em, E_level.Em, choice_hE_fit_allTrials.Em, RT_allTrials.Em] = deal(NaN(nTrialsPerTask, NS));
n_E_levels = 3;
[choice_hE.Ep, choice_hE.Em,...
    fMRI_choice_highE.Ep, fMRI_choice_lowE.Ep,...
    fMRI_choice_highE.Em, fMRI_choice_lowE.Em] = deal(NaN(n_E_levels, NS));

% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
    % extract runs
    [runsStruct] = runs_definition(study_nm, sub_nm, condition);
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    for iRun = 1:length(okRuns)
        kRun = okRuns(iRun);
        task_nm_tmp = taskNames{iRun};
        run_nm = num2str(kRun);
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        % define which task session it is
        switch kRun
            case {1,2}
                taskRun_idx = 1;
            case {3,4}
                taskRun_idx = 2;
        end
        run_nm_bis = ['run',num2str(taskRun_idx)];
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(taskRun_idx-1);
        runTrials_EpEmPool_idx = (1:nTrialsPerRun) + nTrialsPerRun*(kRun-1);
        %% extract choices for the current session
        choice_hE_allTrials.(task_nm_tmp)(runTrials_idx, iS) = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        E_level.(task_nm_tmp)(runTrials_idx, iS) = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        RT_allTrials.(task_nm_tmp)(runTrials_idx, iS) = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        %% extract fMRI ROI
        fMRI_allTrials.(task_nm_tmp)(runTrials_idx, iS) = ROI_trial_b_trial.(fMRI_ROI_trialPerTrial_name).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
    end % run loop
    
    % loop through tasks
    for iT = 1:nTasks
        task_nm_tmp = tasks{iT};
        %% remove RT effect from fMRI ROI
        if RT_orth == 1
            % estimate how ROI is impacted by RT
            beta_tmp.(task_nm_tmp) = glmfit(RT_allTrials.(task_nm_tmp)(:,iS), fMRI_allTrials.(task_nm_tmp)(:,iS),'normal');
            % then orthogonalize ROI activity to RT
            fMRI_allTrials.(task_nm_tmp)(:,iS) = fMRI_allTrials.(task_nm_tmp)(:,iS) - beta_tmp.(task_nm_tmp)(2).*RT_allTrials.(task_nm_tmp)(:,iS);
        end
        
        %% extract the bins split per effort level
        for iE = 1:n_E_levels
            % choice=f(E)
            E_lvl_idx = E_level.(task_nm_tmp)(:, iS) == iE;
            choice_hE.(task_nm_tmp)(iE, iS) = mean(choice_hE_allTrials.(task_nm_tmp)(E_lvl_idx, iS),1,'omitnan');
            % fMRI=f(E*choice high E)
            choice_highE_E_lvl_idx = (E_lvl_idx.*(choice_hE_allTrials.(task_nm_tmp)(:, iS) == 1)) == 1;
            fMRI_choice_highE.(task_nm_tmp)(iE, iS) = mean(fMRI_allTrials.(task_nm_tmp)(choice_highE_E_lvl_idx, iS),1,'omitnan');
            % fMRI=f(E*choice low E)
            choice_lowE_E_lvl_idx = (E_lvl_idx.*(choice_hE_allTrials.(task_nm_tmp)(:, iS) == 0)) == 1;
            fMRI_choice_lowE.(task_nm_tmp)(iE, iS) = mean(fMRI_allTrials.(task_nm_tmp)(choice_lowE_E_lvl_idx, iS),1,'omitnan');
        end
    end % task loop
end % subject loop

%% perform the median split
for iT = 1:nTasks
    task_nm_tmp = tasks{iT};
    % extract choice = f(E) for each group
    [choice_hE.low_con_m_subs.(task_nm_tmp),...
        choice_hE.low_con_sem_subs.(task_nm_tmp)] = mean_sem_sd(choice_hE.(task_nm_tmp)(:,low_con.(task_nm_tmp)),2);
    [choice_hE.high_con_m_subs.(task_nm_tmp),...
        choice_hE.high_con_sem_subs.(task_nm_tmp)] = mean_sem_sd(choice_hE.(task_nm_tmp)(:,high_con.(task_nm_tmp)),2);
    % extract ROI activity = f(E) for each group
    [fMRI_choice_highE.low_con_m_subs.(task_nm_tmp),...
        fMRI_choice_highE.low_con_sem_subs.(task_nm_tmp)] = mean_sem_sd(fMRI_choice_highE.(task_nm_tmp)(:,low_con.(task_nm_tmp)),2);
    [fMRI_choice_highE.high_con_m_subs.(task_nm_tmp),...
        fMRI_choice_highE.high_con_sem_subs.(task_nm_tmp)] = mean_sem_sd(fMRI_choice_highE.(task_nm_tmp)(:,high_con.(task_nm_tmp)),2);
    [fMRI_choice_lowE.low_con_m_subs.(task_nm_tmp),...
        fMRI_choice_lowE.low_con_sem_subs.(task_nm_tmp)] = mean_sem_sd(fMRI_choice_lowE.(task_nm_tmp)(:,low_con.(task_nm_tmp)),2);
    [fMRI_choice_lowE.high_con_m_subs.(task_nm_tmp),...
        fMRI_choice_lowE.high_con_sem_subs.(task_nm_tmp)] = mean_sem_sd(fMRI_choice_lowE.(task_nm_tmp)(:,high_con.(task_nm_tmp)),2);
    
    if figDisp == 1
        %% general figure parameters
        pSize = 30;
        lWidth = 3;
        low_col = [241 163 64]./255;
        high_col = [153 142 195]./255;
        con_to_use_bis = strrep(con_to_use,'_',' ');
        
        fig;
        %% choice = f(E)
        subplot(1,2,1);
        hold on;
        low_con_hdl = errorbar(1:n_E_levels,...
            choice_hE.low_con_m_subs.(task_nm_tmp),...
            choice_hE.low_con_sem_subs.(task_nm_tmp));
        low_con_hdl.Color = low_col;
        low_con_hdl.LineStyle = '-';
        low_con_hdl.LineWidth = lWidth;
        high_con_hdl = errorbar(1:n_E_levels,...
            choice_hE.high_con_m_subs.(task_nm_tmp),...
            choice_hE.high_con_sem_subs.(task_nm_tmp));
        high_con_hdl.Color = high_col;
        high_con_hdl.LineWidth = lWidth;
        legend([low_con_hdl, high_con_hdl],...
            ['low ',con_to_use_bis],['high ',con_to_use_bis]);
        legend('boxoff');
        xlabel('high E level');
        ylabel([task_nm_tmp,' choices (%)']);
        legend_size(pSize);
        
        %% fMRI = f(E)
        subplot(1,2,2);
        hold on;
        % low contrast, high effort chosen
        low_con_highEch_hdl = errorbar(1:n_E_levels,...
            fMRI_choice_highE.low_con_m_subs.(task_nm_tmp),...
            fMRI_choice_highE.low_con_sem_subs.(task_nm_tmp));
        low_con_highEch_hdl.Color = low_col;
        low_con_highEch_hdl.LineStyle = '-';
        low_con_highEch_hdl.LineWidth = lWidth;
        % low contrast, low effort chosen
        low_con_lowEch_hdl = errorbar(1:n_E_levels,...
            fMRI_choice_lowE.low_con_m_subs.(task_nm_tmp),...
            fMRI_choice_lowE.low_con_sem_subs.(task_nm_tmp));
        low_con_lowEch_hdl.Color = low_col;
        low_con_lowEch_hdl.LineStyle = '--';
        low_con_lowEch_hdl.LineWidth = lWidth;
        % high contrast, high effort chosen
        high_con_highEch_hdl = errorbar(1:n_E_levels,...
            fMRI_choice_highE.high_con_m_subs.(task_nm_tmp),...
            fMRI_choice_highE.high_con_sem_subs.(task_nm_tmp));
        high_con_highEch_hdl.Color = high_col;
        high_con_highEch_hdl.LineStyle = '-';
        high_con_highEch_hdl.LineWidth = lWidth;
        % high contrast, low effort chosen
        high_con_lowEch_hdl = errorbar(1:n_E_levels,...
            fMRI_choice_lowE.high_con_m_subs.(task_nm_tmp),...
            fMRI_choice_lowE.high_con_sem_subs.(task_nm_tmp));
        high_con_lowEch_hdl.Color = high_col;
        high_con_lowEch_hdl.LineStyle = '--';
        high_con_lowEch_hdl.LineWidth = lWidth;
        legend([low_con_highEch_hdl, low_con_lowEch_hdl,...
            high_con_highEch_hdl, high_con_lowEch_hdl],...
            ['low ',con_to_use_bis,' high E ch'],...
            ['low ',con_to_use_bis,' low E ch'],...
            ['high ',con_to_use_bis,' high E ch'],...
            ['high ',con_to_use_bis,' low E ch']);
        legend('boxoff');
        xlabel('high E level');
        ylabel([task_nm_tmp,' ',...
            ROI_coords.ROI_nm.ROI_1_shortName,' BOLD']);
        legend_size(pSize);
        
    end % figure display
    
end % task loop