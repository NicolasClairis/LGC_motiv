%% check whether the ROI activity predicts the performance, 
% within each effort level

%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% selection of participants
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end

%% load ROI
[ROI_trial_b_trial, ROI_subList,...
    ROI_nm, ROI_short_nm,...
    task_to_look, timePeriod_nm] = extract_ROI_betas_onsets_only_bis(computerRoot,...
    study_nm, subject_id, condition);

%% general parameters
nRuns = 4;
nTrialsPerRun = 54;
nTrials = nTrialsPerRun*nRuns;
task_names = {'physical','mental'};
nTasks = length(task_names);
[ROI_activity.allTrials,...
    perf.latency.allTrials,...
    perf.AUC.allTrials,...
    perf.forcePeak.allTrials,...
    perf.AUC_overshoot.allTrials,...
    perf.n_errors.allTrials,...
    perf.successSpeed.allTrials,...
    choice_hE,...
    hE_level,...
    E_chosen,...
    task_cstt.physical, task_cstt.mental] = deal(NaN(nTrials, NS));
n_hE_levels = 3;
hE_levels = 1:3;
n_Ech = 4;
Ech_levels = 0:3;
for iTask = 1:nTasks
    task_nm = task_names{iTask};
    [ROI_activity.(task_nm).perEch,...
        perf.(task_nm).latency.perEch.(['low_',ROI_short_nm]),...
        perf.(task_nm).AUC.perEch.(['low_',ROI_short_nm]),...
        perf.(task_nm).forcePeak.perEch.(['low_',ROI_short_nm]),...
        perf.(task_nm).AUC_overshoot.perEch.(['low_',ROI_short_nm]),...
        perf.(task_nm).n_errors.perEch.(['low_',ROI_short_nm]),...
        perf.(task_nm).successSpeed.perEch.(['low_',ROI_short_nm]),...
        perf.(task_nm).latency.perEch.(['high_',ROI_short_nm]),...
        perf.(task_nm).AUC.perEch.(['high_',ROI_short_nm]),...
        perf.(task_nm).forcePeak.perEch.(['high_',ROI_short_nm]),...
        perf.(task_nm).AUC_overshoot.perEch.(['high_',ROI_short_nm]),...
        perf.(task_nm).n_errors.perEch.(['high_',ROI_short_nm]),...
        perf.(task_nm).successSpeed.perEch.(['high_',ROI_short_nm])] = deal(NaN(n_Ech,NS));
    [ROI_activity.(task_nm).choice_highE.perHighElevel,...
        perf.(task_nm).latency.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).AUC.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).forcePeak.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).AUC_overshoot.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).n_errors.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).successSpeed.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).latency.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).AUC.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).forcePeak.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).AUC_overshoot.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).n_errors.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).successSpeed.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
        ROI_activity.(task_nm).choice_lowE.perHighElevel,...
        perf.(task_nm).latency.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).AUC.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).forcePeak.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).AUC_overshoot.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).n_errors.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).successSpeed.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
        perf.(task_nm).latency.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).AUC.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).forcePeak.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).AUC_overshoot.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).n_errors.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
        perf.(task_nm).successSpeed.choice_lowE.perHighElevel.(['high_',ROI_short_nm])] = deal(NaN(n_hE_levels,NS));
end % task loop

%% load performance and extract corresponding ROI activity
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
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(kRun-1);
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
        
        %% load choice informations
        choice_hE(runTrials_idx, iS) = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        hE_level(runTrials_idx, iS) = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        E_chosen(runTrials_idx, iS) = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        task_cstt.(task_fullName)(runTrials_idx,iS) = 1;
        %% load the performance data
        switch task_nm_tmp
            case 'Ep'
                [latency, AUC, forcePeak, AUC_overshoot] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
                perf.latency.allTrials(runTrials_idx,iS) = latency.allTrials;
                perf.AUC.allTrials(runTrials_idx,iS) = AUC.allTrials;
                perf.forcePeak.allTrials(runTrials_idx,iS) = forcePeak.allTrials;
                perf.AUC_overshoot.allTrials(runTrials_idx,iS) = AUC_overshoot.allTrials;
                
            case 'Em'
                [successSpeed, n_errors, RT_avg] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
                perf.successSpeed.allTrials(runTrials_idx,iS) = successSpeed.allTrials;
                perf.n_errors.allTrials(runTrials_idx,iS) = n_errors.allTrials;
                perf.RT_avg.allTrials(runTrials_idx,iS) = RT_avg.allTrials;
        end
        
        %% extract fMRI ROI mediator
        if strcmp(task_to_look,'EpEmPool') ||...
                (strcmp(task_to_look, task_nm_tmp))
            ROI_activity.allTrials(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
        end
    end % run loop
    
    %% split data by effort levels and effort chosen and ROI activity
    for iTask = 1:nTasks
        task_nm = task_names{iTask};
        task_idx = task_cstt.(task_nm)(:,iS);
        for iEch = 1:n_Ech
            jEch_idx = (E_chosen(:,iS) == (iEch - 1)).*(task_idx == 1) == 1;
            ROI_activity.(task_nm).perEch(iEch, iS)       = mean(ROI_activity.allTrials(jEch_idx,iS),1,'omitnan');
            median_ROI_Ech = median(ROI_activity.allTrials(jEch_idx, iS),1,'omitnan');
            low_ROI_Ech = ((ROI_activity.allTrials(:,iS) <= median_ROI_Ech).*jEch_idx) == 1;
            high_ROI_Ech = ((ROI_activity.allTrials(:,iS) > median_ROI_Ech).*jEch_idx) == 1;
            switch task_nm
                case 'physical'
                    % low ROI level
                    perf.(task_nm).latency.perEch.(['low_',ROI_short_nm])(iEch, iS)       = mean(perf.latency.allTrials(low_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).AUC.perEch.(['low_',ROI_short_nm])(iEch, iS)           = mean(perf.AUC.allTrials(low_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).forcePeak.perEch.(['low_',ROI_short_nm])(iEch, iS)     = mean(perf.forcePeak.allTrials(low_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).AUC_overshoot.perEch.(['low_',ROI_short_nm])(iEch, iS) = mean(perf.AUC_overshoot.allTrials(low_ROI_Ech,iS),1,'omitnan');
                    % high ROI level
                    perf.(task_nm).latency.perEch.(['high_',ROI_short_nm])(iEch, iS)       = mean(perf.latency.allTrials(high_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).AUC.perEch.(['high_',ROI_short_nm])(iEch, iS)           = mean(perf.AUC.allTrials(high_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).forcePeak.perEch.(['high_',ROI_short_nm])(iEch, iS)     = mean(perf.forcePeak.allTrials(high_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).AUC_overshoot.perEch.(['high_',ROI_short_nm])(iEch, iS) = mean(perf.AUC_overshoot.allTrials(high_ROI_Ech,iS),1,'omitnan');
                case 'mental'
                    % low ROI level
                    perf.(task_nm).n_errors.perEch.(['low_',ROI_short_nm])(iEch, iS)        = mean(perf.n_errors.allTrials(low_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).successSpeed.perEch.(['low_',ROI_short_nm])(iEch, iS)    = mean(perf.successSpeed.allTrials(low_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).RT_avg.perEch.(['low_',ROI_short_nm])(iEch, iS)          = mean(perf.RT_avg.allTrials(low_ROI_Ech,iS),1,'omitnan');
                    % high ROI level
                    perf.(task_nm).n_errors.perEch.(['high_',ROI_short_nm])(iEch, iS)       = mean(perf.n_errors.allTrials(high_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).successSpeed.perEch.(['high_',ROI_short_nm])(iEch, iS)   = mean(perf.successSpeed.allTrials(high_ROI_Ech,iS),1,'omitnan');
                    perf.(task_nm).RT_avg.perEch.(['high_',ROI_short_nm])(iEch, iS)         = mean(perf.RT_avg.allTrials(high_ROI_Ech,iS),1,'omitnan');
            end
        end % effort chosen
        
        for iEff = 1:n_hE_levels
            jEff_lowEchoice_idx = ((hE_level(:,iS) == iEff).*(choice_hE(:,iS) == 0).*(task_idx == 1)) == 1;
            jEff_highEchoice_idx = ((hE_level(:,iS) == iEff).*(choice_hE(:,iS) == 1).*(task_idx == 1)) == 1;
            ROI_activity.(task_nm).choice_lowE.perHighElevel(iEff, iS)       = mean(ROI_activity.allTrials(jEff_lowEchoice_idx,iS),1,'omitnan');
            ROI_activity.(task_nm).choice_highE.perHighElevel(iEff, iS)       = mean(ROI_activity.allTrials(jEff_highEchoice_idx,iS),1,'omitnan');
            
            median_ROI_hE_lowCh = median(ROI_activity.allTrials(jEff_lowEchoice_idx, iS),1,'omitnan');
            low_ROI_hE_lowCh = ((ROI_activity.allTrials(:,iS) <= median_ROI_hE_lowCh).*jEff_lowEchoice_idx) == 1;
            high_ROI_hE_lowCh = ((ROI_activity.allTrials(:,iS) > median_ROI_hE_lowCh).*jEff_lowEchoice_idx) == 1;
            
            median_ROI_hE_highCh = median(ROI_activity.allTrials(jEff_highEchoice_idx, iS),1,'omitnan');
            low_ROI_hE_highCh = ((ROI_activity.allTrials(:,iS) <= median_ROI_hE_highCh).*jEff_highEchoice_idx) == 1;
            high_ROI_hE_highCh = ((ROI_activity.allTrials(:,iS) > median_ROI_hE_highCh).*jEff_highEchoice_idx) == 1;
            switch task_nm
                case 'physical'
                    % low ROI - low effort chosen
                    perf.(task_nm).latency.choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)       = mean(perf.latency.allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).AUC.choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)           = mean(perf.AUC.allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).forcePeak.choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)     = mean(perf.forcePeak.allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).AUC_overshoot.choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS) = mean(perf.AUC_overshoot.allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                    % high ROI - low effort chosen
                    perf.(task_nm).latency.choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)       = mean(perf.latency.allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).AUC.choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)           = mean(perf.AUC.allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).forcePeak.choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)     = mean(perf.forcePeak.allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).AUC_overshoot.choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS) = mean(perf.AUC_overshoot.allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                    % low ROI - high effort chosen
                    perf.(task_nm).latency.choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)      = mean(perf.latency.allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).AUC.choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)          = mean(perf.AUC.allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).forcePeak.choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)    = mean(perf.forcePeak.allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).AUC_overshoot.choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)= mean(perf.AUC_overshoot.allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                    % high ROI - high effort chosen
                    perf.(task_nm).latency.choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)      = mean(perf.latency.allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).AUC.choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)          = mean(perf.AUC.allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).forcePeak.choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)    = mean(perf.forcePeak.allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).AUC_overshoot.choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)= mean(perf.AUC_overshoot.allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                case 'mental'
                    % low ROI - low effort chosen
                    perf.(task_nm).n_errors.choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)         = mean(perf.n_errors.allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).successSpeed.choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)     = mean(perf.successSpeed.allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).RT_avg.choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)           = mean(perf.RT_avg.allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                    % high ROI - low effort chosen
                    perf.(task_nm).n_errors.choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)        = mean(perf.n_errors.allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).successSpeed.choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)    = mean(perf.successSpeed.allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                    perf.(task_nm).RT_avg.choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)          = mean(perf.RT_avg.allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                    % low ROI - high effort chosen
                    perf.(task_nm).n_errors.choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)        = mean(perf.n_errors.allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).successSpeed.choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)    = mean(perf.successSpeed.allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).RT_avg.choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)          = mean(perf.RT_avg.allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                    % high ROI - high effort chosen
                    perf.(task_nm).n_errors.choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)       = mean(perf.n_errors.allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).successSpeed.choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)   = mean(perf.successSpeed.allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                    perf.(task_nm).RT_avg.choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)         = mean(perf.RT_avg.allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
            end
        end % effort level
        
        %% statistical test and fit eventually
    end % task loop
    
end % subject loop

%% average the data across subjects
for iTask = 1:nTasks
    task_nm = task_names{iTask};
    % effort chosen
    [m_ROI_activity.(task_nm).perEch,...
        sem_ROI_activity.(task_nm).perEch] = mean_sem_sd(ROI_activity.(task_nm).perEch, 2);
    [m_ROI_activity.(task_nm).choice_highE.perHighElevel,...
         sem_ROI_activity.(task_nm).choice_highE.perHighElevel] = mean_sem_sd(ROI_activity.(task_nm).choice_highE.perHighElevel, 2);
     [m_ROI_activity.(task_nm).choice_lowE.perHighElevel,...
         sem_ROI_activity.(task_nm).choice_lowE.perHighElevel] = mean_sem_sd(ROI_activity.(task_nm).choice_lowE.perHighElevel, 2);
end % task loop
     
%% loop through levels of activity of the ROI
ROI_level = {'low','high'};
for iROI_lvl = 1:length(ROI_level)
    ROI_info_nm = [ROI_level{iROI_lvl},'_',ROI_short_nm];
    
    % E chosen
    [m_perf.physical.latency.perEch.(ROI_info_nm),...
        sem_perf.physical.latency.perEch.(ROI_info_nm)] = mean_sem_sd(perf.physical.latency.perEch.(ROI_info_nm), 2);
    [m_perf.physical.AUC.perEch.(ROI_info_nm),...
        sem_perf.physical.AUC.perEch.(ROI_info_nm)] = mean_sem_sd(perf.physical.AUC.perEch.(ROI_info_nm), 2);
    [m_perf.physical.forcePeak.perEch.(ROI_info_nm),...
        sem_perf.physical.forcePeak.perEch.(ROI_info_nm)] = mean_sem_sd(perf.physical.forcePeak.perEch.(ROI_info_nm), 2);
    [m_perf.physical.AUC_overshoot.perEch.(ROI_info_nm),...
        sem_perf.physical.AUC_overshoot.perEch.(ROI_info_nm)] = mean_sem_sd(perf.physical.AUC_overshoot.perEch.(ROI_info_nm), 2);
    [m_perf.mental.n_errors.perEch.(ROI_info_nm),...
        sem_perf.mental.n_errors.perEch.(ROI_info_nm)] = mean_sem_sd(perf.mental.n_errors.perEch.(ROI_info_nm), 2);
    [m_perf.mental.successSpeed.perEch.(ROI_info_nm),...
        sem_perf.mental.successSpeed.perEch.(ROI_info_nm)] = mean_sem_sd(perf.mental.successSpeed.perEch.(ROI_info_nm), 2);
    [m_perf.mental.RT_avg.perEch.(ROI_info_nm),...
        sem_perf.mental.RT_avg.perEch.(ROI_info_nm)] = mean_sem_sd(perf.mental.RT_avg.perEch.(ROI_info_nm), 2);
    
    % low effort chosen
    [m_perf.physical.latency.choice_lowE.perHighElevel.(ROI_info_nm),...
        sem_perf.physical.latency.choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.latency.choice_lowE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.physical.AUC.choice_lowE.perHighElevel.(ROI_info_nm),...
        sem_perf.physical.AUC.choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.AUC.choice_lowE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.physical.forcePeak.choice_lowE.perHighElevel.(ROI_info_nm),...
        sem_perf.physical.forcePeak.choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.forcePeak.choice_lowE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.physical.AUC_overshoot.choice_lowE.perHighElevel.(ROI_info_nm),...
        sem_perf.physical.AUC_overshoot.choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.AUC_overshoot.choice_lowE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.mental.n_errors.choice_lowE.perHighElevel.(ROI_info_nm),...
        sem_perf.mental.n_errors.choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.mental.n_errors.choice_lowE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.mental.successSpeed.choice_lowE.perHighElevel.(ROI_info_nm),...
        sem_perf.mental.successSpeed.choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.mental.successSpeed.choice_lowE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.mental.RT_avg.choice_lowE.perHighElevel.(ROI_info_nm),...
        sem_perf.mental.RT_avg.choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.mental.RT_avg.choice_lowE.perHighElevel.(ROI_info_nm), 2);
    % high effort chosen
    [m_perf.physical.latency.choice_highE.perHighElevel.(ROI_info_nm),...
        sem_perf.physical.latency.choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.latency.choice_highE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.physical.AUC.choice_highE.perHighElevel.(ROI_info_nm),...
        sem_perf.physical.AUC.choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.AUC.choice_highE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.physical.forcePeak.choice_highE.perHighElevel.(ROI_info_nm),...
        sem_perf.physical.forcePeak.choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.forcePeak.choice_highE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.physical.AUC_overshoot.choice_highE.perHighElevel.(ROI_info_nm),...
        sem_perf.physical.AUC_overshoot.choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.AUC_overshoot.choice_highE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.mental.n_errors.choice_highE.perHighElevel.(ROI_info_nm),...
        sem_perf.mental.n_errors.choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.mental.n_errors.choice_highE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.mental.successSpeed.choice_highE.perHighElevel.(ROI_info_nm),...
        sem_perf.mental.successSpeed.choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.mental.successSpeed.choice_highE.perHighElevel.(ROI_info_nm), 2);
    [m_perf.mental.RT_avg.choice_highE.perHighElevel.(ROI_info_nm),...
        sem_perf.mental.RT_avg.choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.mental.RT_avg.choice_highE.perHighElevel.(ROI_info_nm), 2);
end % ROI activity

%% display figures
pSize = 30;
lWidth = 1;
blue = [0 143 255]./255;
purple = [128 0 128]./255;
%% E chosen - grip
fig;

% latency = f(ROI)
subplot(2,2,1);
hold on;
Ech_lat_low_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.latency.perEch.(['low_',ROI_short_nm]),...
    sem_perf.physical.latency.perEch.(['low_',ROI_short_nm]));
Ech_lat_high_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.latency.perEch.(['high_',ROI_short_nm]),...
    sem_perf.physical.latency.perEch.(['high_',ROI_short_nm]));
% modify curve properties
Ech_lat_low_ROI_hdl.LineStyle = '--';
Ech_lat_low_ROI_hdl.LineWidth = lWidth;
Ech_lat_low_ROI_hdl.Color = blue;
Ech_lat_high_ROI_hdl.LineStyle = '-';
Ech_lat_high_ROI_hdl.LineWidth = lWidth;
Ech_lat_high_ROI_hdl.Color = purple;
legend([Ech_lat_high_ROI_hdl, Ech_lat_low_ROI_hdl],...
    {['high ',ROI_short_nm],['low ',ROI_short_nm]});
legend('boxoff');
xticks(Ech_levels);
xlabel('E chosen');
ylabel('Latency to squeeze (s)');
legend_size(pSize);

% force peak = f(ROI)
subplot(2,2,2);
hold on;
Ech_Fpeak_low_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.forcePeak.perEch.(['low_',ROI_short_nm]),...
    sem_perf.physical.forcePeak.perEch.(['low_',ROI_short_nm]));
Ech_Fpeak_high_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.forcePeak.perEch.(['high_',ROI_short_nm]),...
    sem_perf.physical.forcePeak.perEch.(['high_',ROI_short_nm]));
% modify curve properties
Ech_Fpeak_low_ROI_hdl.LineStyle = '--';
Ech_Fpeak_low_ROI_hdl.LineWidth = lWidth;
Ech_Fpeak_low_ROI_hdl.Color = blue;
Ech_Fpeak_high_ROI_hdl.LineStyle = '-';
Ech_Fpeak_high_ROI_hdl.LineWidth = lWidth;
Ech_Fpeak_high_ROI_hdl.Color = purple;
legend([Ech_Fpeak_high_ROI_hdl, Ech_Fpeak_low_ROI_hdl],...
    {['high ',ROI_short_nm],['low ',ROI_short_nm]});
legend('boxoff');
xticks(Ech_levels);
xlabel('E chosen');
ylabel('Peak force (%)');
legend_size(pSize);

% AUC force = f(ROI)
subplot(2,2,3);
hold on;
Ech_AUC_low_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.AUC.perEch.(['low_',ROI_short_nm]),...
    sem_perf.physical.AUC.perEch.(['low_',ROI_short_nm]));
Ech_AUC_high_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.AUC.perEch.(['high_',ROI_short_nm]),...
    sem_perf.physical.AUC.perEch.(['high_',ROI_short_nm]));
% modify curve properties
Ech_AUC_low_ROI_hdl.LineStyle = '--';
Ech_AUC_low_ROI_hdl.LineWidth = lWidth;
Ech_AUC_low_ROI_hdl.Color = blue;
Ech_AUC_high_ROI_hdl.LineStyle = '-';
Ech_AUC_high_ROI_hdl.LineWidth = lWidth;
Ech_AUC_high_ROI_hdl.Color = purple;
legend([Ech_AUC_high_ROI_hdl, Ech_AUC_low_ROI_hdl],...
    {['high ',ROI_short_nm],['low ',ROI_short_nm]});
legend('boxoff');
xticks(Ech_levels);
xlabel('E chosen');
ylabel('AUC force');
legend_size(pSize);

% AUC overshoot force = f(ROI)
subplot(2,2,4);
hold on;
Ech_AUC_overshoot_low_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.AUC_overshoot.perEch.(['low_',ROI_short_nm]),...
    sem_perf.physical.AUC_overshoot.perEch.(['low_',ROI_short_nm]));
Ech_AUC_overshoot_high_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.AUC_overshoot.perEch.(['high_',ROI_short_nm]),...
    sem_perf.physical.AUC_overshoot.perEch.(['high_',ROI_short_nm]));
% modify curve properties
Ech_AUC_overshoot_low_ROI_hdl.LineStyle = '--';
Ech_AUC_overshoot_low_ROI_hdl.LineWidth = lWidth;
Ech_AUC_overshoot_low_ROI_hdl.Color = blue;
Ech_AUC_overshoot_high_ROI_hdl.LineStyle = '-';
Ech_AUC_overshoot_high_ROI_hdl.LineWidth = lWidth;
Ech_AUC_overshoot_high_ROI_hdl.Color = purple;
legend([Ech_AUC_overshoot_high_ROI_hdl, Ech_AUC_overshoot_low_ROI_hdl],...
    {['high ',ROI_short_nm],['low ',ROI_short_nm]});
legend('boxoff');
xticks(Ech_levels);
xlabel('E chosen');
ylabel('AUC force overshoot');
legend_size(pSize);

%% E chosen - 2-back (not ready yet)
fig;

% success speed = f(ROI)
subplot(1,3,1);
hold on;
Ech_successSpeed_low_ROI_hdl = errorbar(Ech_levels,...
    m_perf.mental.successSpeed.perEch.(['low_',ROI_short_nm]),...
    sem_perf.mental.successSpeed.perEch.(['low_',ROI_short_nm]));
Ech_successSpeed_high_ROI_hdl = errorbar(Ech_levels,...
    m_perf.mental.successSpeed.perEch.(['high_',ROI_short_nm]),...
    sem_perf.mental.successSpeed.perEch.(['high_',ROI_short_nm]));
% modify curve properties
Ech_successSpeed_low_ROI_hdl.LineStyle = '--';
Ech_successSpeed_low_ROI_hdl.LineWidth = lWidth;
Ech_successSpeed_low_ROI_hdl.Color = blue;
Ech_successSpeed_high_ROI_hdl.LineStyle = '-';
Ech_successSpeed_high_ROI_hdl.LineWidth = lWidth;
Ech_successSpeed_high_ROI_hdl.Color = purple;
legend([Ech_successSpeed_high_ROI_hdl, Ech_successSpeed_low_ROI_hdl],...
    {['high ',ROI_short_nm],['low ',ROI_short_nm]});
legend('boxoff');
xticks(Ech_levels);
xlabel('E chosen');
ylabel('success speed (s)');
legend_size(pSize);

% number of errors made = f(ROI)
subplot(1,3,2);
hold on;
Ech_nErrors_low_ROI_hdl = errorbar(Ech_levels,...
    m_perf.mental.n_errors.perEch.(['low_',ROI_short_nm]),...
    sem_perf.mental.n_errors.perEch.(['low_',ROI_short_nm]));
Ech_nErrors_high_ROI_hdl = errorbar(Ech_levels,...
    m_perf.mental.n_errors.perEch.(['high_',ROI_short_nm]),...
    sem_perf.mental.n_errors.perEch.(['high_',ROI_short_nm]));
% modify curve properties
Ech_nErrors_low_ROI_hdl.LineStyle = '--';
Ech_nErrors_low_ROI_hdl.LineWidth = lWidth;
Ech_nErrors_low_ROI_hdl.Color = blue;
Ech_nErrors_high_ROI_hdl.LineStyle = '-';
Ech_nErrors_high_ROI_hdl.LineWidth = lWidth;
Ech_nErrors_high_ROI_hdl.Color = purple;
legend([Ech_nErrors_high_ROI_hdl, Ech_nErrors_low_ROI_hdl],...
    {['high ',ROI_short_nm],['low ',ROI_short_nm]});
legend('boxoff');
xticks(Ech_levels);
xlabel('E chosen');
ylabel('Number of errors');
legend_size(pSize);

% average(RT) = f(ROI)
subplot(1,3,3);
hold on;
Ech_RT_avg_low_ROI_hdl = errorbar(Ech_levels,...
    m_perf.mental.RT_avg.perEch.(['low_',ROI_short_nm]),...
    sem_perf.mental.RT_avg.perEch.(['low_',ROI_short_nm]));
Ech_RT_avg_high_ROI_hdl = errorbar(Ech_levels,...
    m_perf.mental.RT_avg.perEch.(['high_',ROI_short_nm]),...
    sem_perf.mental.RT_avg.perEch.(['high_',ROI_short_nm]));
% modify curve properties
Ech_RT_avg_low_ROI_hdl.LineStyle = '--';
Ech_RT_avg_low_ROI_hdl.LineWidth = lWidth;
Ech_RT_avg_low_ROI_hdl.Color = blue;
Ech_RT_avg_high_ROI_hdl.LineStyle = '-';
Ech_RT_avg_high_ROI_hdl.LineWidth = lWidth;
Ech_RT_avg_high_ROI_hdl.Color = purple;
legend([Ech_RT_avg_high_ROI_hdl, Ech_RT_avg_low_ROI_hdl],...
    {['high ',ROI_short_nm],['low ',ROI_short_nm]});
legend('boxoff');
xticks(Ech_levels);
xlabel('E chosen');
ylabel('avg(RT) (s)');
legend_size(pSize);

%% effort level and choice - grip
fig;

% latency = f(ROI)
subplot(2,2,1);
hold on;
hE_lowEchoice_lat_lowROI_hdl = errorbar(hE_levels,...
    m_perf.physical.latency.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.physical.latency.choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
hE_lowEchoice_lat_highROI_hdl = errorbar(hE_levels,...
    m_perf.physical.latency.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.physical.latency.choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
hE_highEchoice_lat_lowROI_hdl = errorbar(hE_levels,...
    m_perf.physical.latency.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.physical.latency.choice_highE.perHighElevel.(['low_',ROI_short_nm]));
hE_highEchoice_lat_highROI_hdl = errorbar(hE_levels,...
    m_perf.physical.latency.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.physical.latency.choice_highE.perHighElevel.(['high_',ROI_short_nm]));
% modify curve properties
hE_lowEchoice_lat_lowROI_hdl.LineStyle = '--';
hE_lowEchoice_lat_lowROI_hdl.LineWidth = lWidth;
hE_lowEchoice_lat_lowROI_hdl.Color = blue;
hE_lowEchoice_lat_highROI_hdl.LineStyle = '--';
hE_lowEchoice_lat_highROI_hdl.LineWidth = lWidth;
hE_lowEchoice_lat_highROI_hdl.Color = purple;
hE_highEchoice_lat_lowROI_hdl.LineStyle = '-';
hE_highEchoice_lat_lowROI_hdl.LineWidth = lWidth;
hE_highEchoice_lat_lowROI_hdl.Color = blue;
hE_highEchoice_lat_highROI_hdl.LineStyle = '-';
hE_highEchoice_lat_highROI_hdl.LineWidth = lWidth;
hE_highEchoice_lat_highROI_hdl.Color = purple;
legend([hE_highEchoice_lat_highROI_hdl,...
    hE_highEchoice_lat_lowROI_hdl,...
    hE_lowEchoice_lat_highROI_hdl,...
    hE_lowEchoice_lat_lowROI_hdl],...
    {['hEch - high ',ROI_short_nm],...
    ['hEch - low ',ROI_short_nm],...
    ['lowEch - high ',ROI_short_nm],...
    ['lowEch - low ',ROI_short_nm]});
legend('boxoff');
xticks(hE_levels);
xlabel('high effort level');
ylabel('Latency to squeeze (s)');
legend_size(pSize);

% force peak = f(ROI)
subplot(2,2,2);
hold on;
hE_lowEchoice_Fpeak_lowROI_hdl = errorbar(hE_levels,...
    m_perf.physical.forcePeak.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.physical.forcePeak.choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
hE_lowEchoice_Fpeak_highROI_hdl = errorbar(hE_levels,...
    m_perf.physical.forcePeak.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.physical.forcePeak.choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
hE_highEchoice_Fpeak_lowROI_hdl = errorbar(hE_levels,...
    m_perf.physical.forcePeak.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.physical.forcePeak.choice_highE.perHighElevel.(['low_',ROI_short_nm]));
hE_highEchoice_Fpeak_highROI_hdl = errorbar(hE_levels,...
    m_perf.physical.forcePeak.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.physical.forcePeak.choice_highE.perHighElevel.(['high_',ROI_short_nm]));
% modify curve properties
hE_lowEchoice_Fpeak_lowROI_hdl.LineStyle = '--';
hE_lowEchoice_Fpeak_lowROI_hdl.LineWidth = lWidth;
hE_lowEchoice_Fpeak_lowROI_hdl.Color = blue;
hE_lowEchoice_Fpeak_highROI_hdl.LineStyle = '--';
hE_lowEchoice_Fpeak_highROI_hdl.LineWidth = lWidth;
hE_lowEchoice_Fpeak_highROI_hdl.Color = purple;
hE_highEchoice_Fpeak_lowROI_hdl.LineStyle = '-';
hE_highEchoice_Fpeak_lowROI_hdl.LineWidth = lWidth;
hE_highEchoice_Fpeak_lowROI_hdl.Color = blue;
hE_highEchoice_Fpeak_highROI_hdl.LineStyle = '-';
hE_highEchoice_Fpeak_highROI_hdl.LineWidth = lWidth;
hE_highEchoice_Fpeak_highROI_hdl.Color = purple;
legend([hE_highEchoice_Fpeak_highROI_hdl,...
    hE_highEchoice_Fpeak_lowROI_hdl,...
    hE_lowEchoice_Fpeak_highROI_hdl,...
    hE_lowEchoice_Fpeak_lowROI_hdl],...
    {['hEch - high ',ROI_short_nm],...
    ['hEch - low ',ROI_short_nm],...
    ['lowEch - high ',ROI_short_nm],...
    ['lowEch - low ',ROI_short_nm]});
legend('boxoff');
xticks(hE_levels);
xlabel('high effort level');
ylabel('Peak force (%)');
legend_size(pSize);

% AUC force = f(ROI)
subplot(2,2,3);
hold on;
hE_lowEchoice_AUC_lowROI_hdl = errorbar(hE_levels,...
    m_perf.physical.AUC.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.physical.AUC.choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
hE_lowEchoice_AUC_highROI_hdl = errorbar(hE_levels,...
    m_perf.physical.AUC.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.physical.AUC.choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
hE_highEchoice_AUC_lowROI_hdl = errorbar(hE_levels,...
    m_perf.physical.AUC.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.physical.AUC.choice_highE.perHighElevel.(['low_',ROI_short_nm]));
hE_highEchoice_AUC_highROI_hdl = errorbar(hE_levels,...
    m_perf.physical.AUC.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.physical.AUC.choice_highE.perHighElevel.(['high_',ROI_short_nm]));
% modify curve properties
hE_lowEchoice_AUC_lowROI_hdl.LineStyle = '--';
hE_lowEchoice_AUC_lowROI_hdl.LineWidth = lWidth;
hE_lowEchoice_AUC_lowROI_hdl.Color = blue;
hE_lowEchoice_AUC_highROI_hdl.LineStyle = '--';
hE_lowEchoice_AUC_highROI_hdl.LineWidth = lWidth;
hE_lowEchoice_AUC_highROI_hdl.Color = purple;
hE_highEchoice_AUC_lowROI_hdl.LineStyle = '-';
hE_highEchoice_AUC_lowROI_hdl.LineWidth = lWidth;
hE_highEchoice_AUC_lowROI_hdl.Color = blue;
hE_highEchoice_AUC_highROI_hdl.LineStyle = '-';
hE_highEchoice_AUC_highROI_hdl.LineWidth = lWidth;
hE_highEchoice_AUC_highROI_hdl.Color = purple;
legend([hE_highEchoice_AUC_highROI_hdl,...
    hE_highEchoice_AUC_lowROI_hdl,...
    hE_lowEchoice_AUC_highROI_hdl,...
    hE_lowEchoice_AUC_lowROI_hdl],...
    {['hEch - high ',ROI_short_nm],...
    ['hEch - low ',ROI_short_nm],...
    ['lowEch - high ',ROI_short_nm],...
    ['lowEch - low ',ROI_short_nm]});
legend('boxoff');
xticks(hE_levels);
xlabel('high effort level');
ylabel('AUC force');
legend_size(pSize);

% AUC overshoot force = f(ROI)
subplot(2,2,4);
hold on;
hE_lowEchoice_AUC_overshoot_lowROI_hdl = errorbar(hE_levels,...
    m_perf.physical.AUC_overshoot.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.physical.AUC_overshoot.choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
hE_lowEchoice_AUC_overshoot_highROI_hdl = errorbar(hE_levels,...
    m_perf.physical.AUC_overshoot.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.physical.AUC_overshoot.choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
hE_highEchoice_AUC_overshoot_lowROI_hdl = errorbar(hE_levels,...
    m_perf.physical.AUC_overshoot.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.physical.AUC_overshoot.choice_highE.perHighElevel.(['low_',ROI_short_nm]));
hE_highEchoice_AUC_overshoot_highROI_hdl = errorbar(hE_levels,...
    m_perf.physical.AUC_overshoot.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.physical.AUC_overshoot.choice_highE.perHighElevel.(['high_',ROI_short_nm]));
% modify curve properties
hE_lowEchoice_AUC_overshoot_lowROI_hdl.LineStyle = '--';
hE_lowEchoice_AUC_overshoot_lowROI_hdl.LineWidth = lWidth;
hE_lowEchoice_AUC_overshoot_lowROI_hdl.Color = blue;
hE_lowEchoice_AUC_overshoot_highROI_hdl.LineStyle = '--';
hE_lowEchoice_AUC_overshoot_highROI_hdl.LineWidth = lWidth;
hE_lowEchoice_AUC_overshoot_highROI_hdl.Color = purple;
hE_highEchoice_AUC_overshoot_lowROI_hdl.LineStyle = '-';
hE_highEchoice_AUC_overshoot_lowROI_hdl.LineWidth = lWidth;
hE_highEchoice_AUC_overshoot_lowROI_hdl.Color = blue;
hE_highEchoice_AUC_overshoot_highROI_hdl.LineStyle = '-';
hE_highEchoice_AUC_overshoot_highROI_hdl.LineWidth = lWidth;
hE_highEchoice_AUC_overshoot_highROI_hdl.Color = purple;
legend([hE_highEchoice_AUC_overshoot_highROI_hdl,...
    hE_highEchoice_AUC_overshoot_lowROI_hdl,...
    hE_lowEchoice_AUC_overshoot_highROI_hdl,...
    hE_lowEchoice_AUC_overshoot_lowROI_hdl],...
    {['hEch - high ',ROI_short_nm],...
    ['hEch - low ',ROI_short_nm],...
    ['lowEch - high ',ROI_short_nm],...
    ['lowEch - low ',ROI_short_nm]});
legend('boxoff');
xticks(hE_levels);
xlabel('high effort level');
ylabel('AUC force overshoot');
legend_size(pSize);

%% effort level and choice - 2-back
fig;

% success speed = f(ROI)
subplot(1,3,1);
hold on;
hE_lowEchoice_successSpeed_lowROI_hdl = errorbar(hE_levels,...
    m_perf.mental.successSpeed.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.mental.successSpeed.choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
hE_lowEchoice_successSpeed_highROI_hdl = errorbar(hE_levels,...
    m_perf.mental.successSpeed.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.mental.successSpeed.choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
hE_highEchoice_successSpeed_lowROI_hdl = errorbar(hE_levels,...
    m_perf.mental.successSpeed.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.mental.successSpeed.choice_highE.perHighElevel.(['low_',ROI_short_nm]));
hE_highEchoice_successSpeed_highROI_hdl = errorbar(hE_levels,...
    m_perf.mental.successSpeed.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.mental.successSpeed.choice_highE.perHighElevel.(['high_',ROI_short_nm]));
% modify curve properties
hE_lowEchoice_successSpeed_lowROI_hdl.LineStyle = '--';
hE_lowEchoice_successSpeed_lowROI_hdl.LineWidth = lWidth;
hE_lowEchoice_successSpeed_lowROI_hdl.Color = blue;
hE_lowEchoice_successSpeed_highROI_hdl.LineStyle = '--';
hE_lowEchoice_successSpeed_highROI_hdl.LineWidth = lWidth;
hE_lowEchoice_successSpeed_highROI_hdl.Color = purple;
hE_highEchoice_successSpeed_lowROI_hdl.LineStyle = '-';
hE_highEchoice_successSpeed_lowROI_hdl.LineWidth = lWidth;
hE_highEchoice_successSpeed_lowROI_hdl.Color = blue;
hE_highEchoice_successSpeed_highROI_hdl.LineStyle = '-';
hE_highEchoice_successSpeed_highROI_hdl.LineWidth = lWidth;
hE_highEchoice_successSpeed_highROI_hdl.Color = purple;
legend([hE_highEchoice_successSpeed_highROI_hdl,...
    hE_highEchoice_successSpeed_lowROI_hdl,...
    hE_lowEchoice_successSpeed_highROI_hdl,...
    hE_lowEchoice_successSpeed_lowROI_hdl],...
    {['hEch - high ',ROI_short_nm],...
    ['hEch - low ',ROI_short_nm],...
    ['lowEch - high ',ROI_short_nm],...
    ['lowEch - low ',ROI_short_nm]});
legend('boxoff');
xticks(hE_levels);
xlabel('high effort level');
ylabel('successSpeed to squeeze (s)');
legend_size(pSize);

% number of errors = f(ROI)
subplot(1,3,2);
hold on;
hE_lowEchoice_nErrors_lowROI_hdl = errorbar(hE_levels,...
    m_perf.mental.n_errors.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.mental.n_errors.choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
hE_lowEchoice_nErrors_highROI_hdl = errorbar(hE_levels,...
    m_perf.mental.n_errors.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.mental.n_errors.choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
hE_highEchoice_nErrors_lowROI_hdl = errorbar(hE_levels,...
    m_perf.mental.n_errors.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.mental.n_errors.choice_highE.perHighElevel.(['low_',ROI_short_nm]));
hE_highEchoice_nErrors_highROI_hdl = errorbar(hE_levels,...
    m_perf.mental.n_errors.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.mental.n_errors.choice_highE.perHighElevel.(['high_',ROI_short_nm]));
% modify curve properties
hE_lowEchoice_nErrors_lowROI_hdl.LineStyle = '--';
hE_lowEchoice_nErrors_lowROI_hdl.LineWidth = lWidth;
hE_lowEchoice_nErrors_lowROI_hdl.Color = blue;
hE_lowEchoice_nErrors_highROI_hdl.LineStyle = '--';
hE_lowEchoice_nErrors_highROI_hdl.LineWidth = lWidth;
hE_lowEchoice_nErrors_highROI_hdl.Color = purple;
hE_highEchoice_nErrors_lowROI_hdl.LineStyle = '-';
hE_highEchoice_nErrors_lowROI_hdl.LineWidth = lWidth;
hE_highEchoice_nErrors_lowROI_hdl.Color = blue;
hE_highEchoice_nErrors_highROI_hdl.LineStyle = '-';
hE_highEchoice_nErrors_highROI_hdl.LineWidth = lWidth;
hE_highEchoice_nErrors_highROI_hdl.Color = purple;
legend([hE_highEchoice_nErrors_highROI_hdl,...
    hE_highEchoice_nErrors_lowROI_hdl,...
    hE_lowEchoice_nErrors_highROI_hdl,...
    hE_lowEchoice_nErrors_lowROI_hdl],...
    {['hEch - high ',ROI_short_nm],...
    ['hEch - low ',ROI_short_nm],...
    ['lowEch - high ',ROI_short_nm],...
    ['lowEch - low ',ROI_short_nm]});
legend('boxoff');
xticks(hE_levels);
xlabel('high effort level');
ylabel('Peak force (%)');
legend_size(pSize);

% RT_avg force = f(ROI)
subplot(1,3,3);
hold on;
hE_lowEchoice_RT_avg_lowROI_hdl = errorbar(hE_levels,...
    m_perf.mental.RT_avg.choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.mental.RT_avg.choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
hE_lowEchoice_RT_avg_highROI_hdl = errorbar(hE_levels,...
    m_perf.mental.RT_avg.choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.mental.RT_avg.choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
hE_highEchoice_RT_avg_lowROI_hdl = errorbar(hE_levels,...
    m_perf.mental.RT_avg.choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
    sem_perf.mental.RT_avg.choice_highE.perHighElevel.(['low_',ROI_short_nm]));
hE_highEchoice_RT_avg_highROI_hdl = errorbar(hE_levels,...
    m_perf.mental.RT_avg.choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
    sem_perf.mental.RT_avg.choice_highE.perHighElevel.(['high_',ROI_short_nm]));
% modify curve properties
hE_lowEchoice_RT_avg_lowROI_hdl.LineStyle = '--';
hE_lowEchoice_RT_avg_lowROI_hdl.LineWidth = lWidth;
hE_lowEchoice_RT_avg_lowROI_hdl.Color = blue;
hE_lowEchoice_RT_avg_highROI_hdl.LineStyle = '--';
hE_lowEchoice_RT_avg_highROI_hdl.LineWidth = lWidth;
hE_lowEchoice_RT_avg_highROI_hdl.Color = purple;
hE_highEchoice_RT_avg_lowROI_hdl.LineStyle = '-';
hE_highEchoice_RT_avg_lowROI_hdl.LineWidth = lWidth;
hE_highEchoice_RT_avg_lowROI_hdl.Color = blue;
hE_highEchoice_RT_avg_highROI_hdl.LineStyle = '-';
hE_highEchoice_RT_avg_highROI_hdl.LineWidth = lWidth;
hE_highEchoice_RT_avg_highROI_hdl.Color = purple;
legend([hE_highEchoice_RT_avg_highROI_hdl,...
    hE_highEchoice_RT_avg_lowROI_hdl,...
    hE_lowEchoice_RT_avg_highROI_hdl,...
    hE_lowEchoice_RT_avg_lowROI_hdl],...
    {['hEch - high ',ROI_short_nm],...
    ['hEch - low ',ROI_short_nm],...
    ['lowEch - high ',ROI_short_nm],...
    ['lowEch - low ',ROI_short_nm]});
legend('boxoff');
xticks(hE_levels);
xlabel('high effort level');
ylabel('RT_avg force');
legend_size(pSize);