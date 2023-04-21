% perf_and_efficacy_f_RT_choice will test whether there is a direct
% correlation between reaction times during choices and subsequent
% performance (0-100%) and efficacy (AUC/latency/peakForce for grip and
% (nCorrect-nErrors)/time for N-back), split by effort chosen (and money
% chosen)


%% subject identification
[study_nm,~,~,subject_id,NS] = sub_id;
%% working directories
computerRoot = LGCM_root_paths;

%% main parameters
task_names = {'Ep','Em'};
nTasks = length(task_names);
nTrialsPerRun = 54;
nRunsPerTask = 2;
nTotalTrialsPerTask = nTrialsPerRun.*nRunsPerTask;
Ech_levels = 0:3;
absMoney_levels = 0:3;
nBins = 6;
[perf.Ep.all, AUC.Ep.all, peakF.Ep.all,...
    latency.Ep.all, AUC_overshoot.Ep.all,...
    RT.Ep.all,...
    Ech.Ep.all, RP.Ep.all, money_level_ch.Ep.all,...
    perf.Em.all, nCorrect.Em.all, nErrors.Em.all,...
    totalTime.Em.all,...
    efficacy_with2first.Em.all,...
    efficacy_pureNback.Em.all,...
    RT.Em.all,...
    Ech.Em.all, RP.Em.all, money_level_ch.Em.all] = deal(NaN(nTotalTrialsPerTask, NS));
[perf.Ep.RTbins, AUC.Ep.RTbins, peakF.Ep.RTbins,...
    latency.Ep.RTbins, AUC_overshoot.Ep.RTbins,...
    RT.Ep.RTbins,...
    perf.Em.RTbins, nCorrect.Em.RTbins, nErrors.Em.RTbins, totalTime.Em.RTbins,...
    efficacy_with2first.Em.RTbins,...
    efficacy_pureNback.Em.RTbins,...
    RT.Em.RTbins] = deal(NaN(nBins, NS));
for iEch = Ech_levels
    Ech_nm_init = ['Ech',num2str(iEch)];
    [perf.Ep.(Ech_nm_init).RTbins, AUC.Ep.(Ech_nm_init).RTbins,...
        peakF.Ep.(Ech_nm_init).RTbins,...
        latency.Ep.(Ech_nm_init).RTbins,...
        AUC_overshoot.Ep.(Ech_nm_init).RTbins,...
        RT.Ep.(Ech_nm_init).RTbins,...
        perf.Em.(Ech_nm_init).RTbins,...
        nCorrect.Em.(Ech_nm_init).RTbins,...
        nErrors.Em.(Ech_nm_init).RTbins,...
        totalTime.Em.(Ech_nm_init).RTbins,...
        efficacy_with2first.Em.(Ech_nm_init).RTbins,...
        efficacy_pureNback.Em.(Ech_nm_init).RTbins,...
        RT.Em.(Ech_nm_init).RTbins] = deal(NaN(nBins, NS));
    for iMoney = absMoney_levels
        money_nm_init = ['M',num2str(iMoney)];
        [perf.Ep.(Ech_nm_init).R.(money_nm_init),...
            AUC.Ep.(Ech_nm_init).R.(money_nm_init),...
            peakF.Ep.(Ech_nm_init).R.(money_nm_init),...
            latency.Ep.(Ech_nm_init).R.(money_nm_init),...
            AUC_overshoot.Ep.(Ech_nm_init).R.(money_nm_init),...
            RT.Ep.(Ech_nm_init).R.(money_nm_init),...
            perf.Em.(Ech_nm_init).R.(money_nm_init),...
            nCorrect.Em.(Ech_nm_init).R.(money_nm_init),...
            nErrors.Em.(Ech_nm_init).R.(money_nm_init),...
            totalTime.Em.(Ech_nm_init).R.(money_nm_init),...
            efficacy_with2first.Em.(Ech_nm_init).R.(money_nm_init),...
            efficacy_pureNback.Em.(Ech_nm_init).R.(money_nm_init),...
            RT.Em.(Ech_nm_init).R.(money_nm_init),...
            perf.Ep.(Ech_nm_init).P.(money_nm_init),...
            AUC.Ep.(Ech_nm_init).P.(money_nm_init),...
            peakF.Ep.(Ech_nm_init).P.(money_nm_init),...
            latency.Ep.(Ech_nm_init).P.(money_nm_init),...
            AUC_overshoot.Ep.(Ech_nm_init).P.(money_nm_init),...
            RT.Ep.(Ech_nm_init).P.(money_nm_init),...
            perf.Em.(Ech_nm_init).P.(money_nm_init),...
            nCorrect.Em.(Ech_nm_init).P.(money_nm_init),...
            nErrors.Em.(Ech_nm_init).P.(money_nm_init),...
            totalTime.Em.(Ech_nm_init).P.(money_nm_init),...
            efficacy_with2first.Em.(Ech_nm_init).P.(money_nm_init),...
            efficacy_pureNback.Em.(Ech_nm_init).P.(money_nm_init),...
            RT.Em.(Ech_nm_init).P.(money_nm_init)] = deal(NaN(nBins, NS));
    end % money
end % effort level

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];
    % extract information of runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    nTotalRuns = length(runsStruct.tasks);

    for iTask = 1:nTasks
        task_id = task_names{iTask};
        switch task_id
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end
        runs.(task_id) = strcmp(runsStruct.tasks, task_id);
        jRun = 0;
        for iRun = 1:nTotalRuns
            runToInclude = 0;
            if runs.(task_id)(iRun) == 1
                jRun = jRun + 1;
                runToInclude = 1;
            end
            run_nm = num2str(runsStruct.runsToKeep(iRun));

            if runToInclude == 1
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);

                %% extract relevant variables depending on task
                switch task_id
                    case 'Ep'
                         [latency_tmp,...
                             AUC_tmp,...
                             peakF_tmp,...
                             AUC_overshoot_tmp] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
                         latency.Ep.all(runTrials_idx,iS) = latency_tmp.allTrials;
                         AUC.Ep.all(runTrials_idx,iS) = AUC_tmp.allTrials;
                         peakF.Ep.all(runTrials_idx,iS) = peakF_tmp.allTrials;
                         AUC_overshoot.Ep.all(runTrials_idx,iS) = AUC_overshoot_tmp.allTrials;
                    case 'Em'
                        [~, nCorrect_tmp,...
                            nErrors_tmp, ~,...
                            totalTime_tmp,...
                            ~,...
                            efficacy_with2first_tmp,...
                            efficacy_pureNback_tmp] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
                        nCorrect.Em.all(runTrials_idx,iS) = nCorrect_tmp.allTrials;
                        nErrors.Em.all(runTrials_idx,iS) = nErrors_tmp.allTrials;
                        totalTime.Em.all(runTrials_idx,iS) = totalTime_tmp.allTrials;
                        efficacy_with2first.Em.all(runTrials_idx,iS) = efficacy_with2first_tmp.allTrials;
                        efficacy_pureNback.Em.all(runTrials_idx,iS) = efficacy_pureNback_tmp.allTrials;
                end
                %% load data
                behaviorStruct_tmp = load([subBehaviorFolder,...
                    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
                    '_task.mat']);
                % load relevant data
                switch task_id
                    case 'Ep'
                        onsets_tmp = behaviorStruct_tmp.physicalPerf.onsets;
                        choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                    case 'Em'
                        onsets_tmp = behaviorStruct_tmp.mentalE_perf.onsets;
                        choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
                end
                %% extract RT
                RT_tmp = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
                RT_tmp(choice_LR_tmp == 0) = NaN;
                RT.(task_id).all(runTrials_idx,iS) = RT_tmp;
                %% extract performance
                for iTrial = 1:nTrialsPerRun
                    jTrial = iTrial +  nTrialsPerRun*(jRun - 1);
                    switch task_id
                        case 'Ep'
                            perf.Ep.all(jTrial,iS) = behaviorStruct_tmp.physicalPerf.perfSummary{1,iTrial}.performance;
                        case 'Em'
                            perf.Em.all(jTrial,iS) = behaviorStruct_tmp.mentalE_perf.perfSummary{1,iTrial}.performance;
                    end
                end
                %% extract effort chosen
                [E_chosen_tmp] = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                Ech.(task_id).all(runTrials_idx,iS) = E_chosen_tmp;
                %% extract R/P and money levels
                [~, money_level_chosen] = extract_money_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                switch task_id
                    case 'Ep'
                        R_trials = strcmp(behaviorStruct_tmp.physicalPerf.choiceOptions.R_or_P,'R');
                    case 'Em'
                        R_trials = strcmp(behaviorStruct_tmp.mentalE_perf.choiceOptions.R_or_P,'R');
                end
                RP.(task_id).all(runTrials_idx,iS) = R_trials;
                money_level_ch.(task_id).all(runTrials_idx,iS) = money_level_chosen.*R_trials +...
                    -money_level_chosen.*(R_trials == 0);
            end % run to include?
        end % run loop

        %% extract bins
        switch task_id
            case 'Ep'
                [perf.Ep.RTbins(:,iS), RT.Ep.RTbins(:,iS)] = do_bin2(perf.Ep.all(:,iS), RT.Ep.all(:,iS), nBins,0);
                AUC.Ep.RTbins(:,iS) = do_bin2(AUC.Ep.all(:,iS), RT.Ep.all(:,iS), nBins,0);
                peakF.Ep.RTbins(:,iS) = do_bin2(peakF.Ep.all(:,iS), RT.Ep.all(:,iS), nBins,0);
                latency.Ep.RTbins(:,iS) = do_bin2(latency.Ep.all(:,iS), RT.Ep.all(:,iS), nBins,0);
                AUC_overshoot.Ep.RTbins(:,iS) = do_bin2(AUC_overshoot.Ep.all(:,iS), RT.Ep.all(:,iS), nBins,0);
            case 'Em'
                [perf.Em.RTbins(:,iS), RT.Em.RTbins(:,iS)] = do_bin2(perf.Em.all(:,iS), RT.Em.all(:,iS), nBins,0);
                nCorrect.Em.RTbins(:,iS) = do_bin2(nCorrect.Em.all(:,iS), RT.Em.all(:,iS), nBins,0);
                nErrors.Em.RTbins(:,iS) = do_bin2(nErrors.Em.all(:,iS), RT.Em.all(:,iS), nBins,0);
                totalTime.Em.RTbins(:,iS) = do_bin2(totalTime.Em.all(:,iS), RT.Em.all(:,iS), nBins,0);
                efficacy_with2first.Em.RTbins(:,iS) = do_bin2(efficacy_with2first.Em.all(:,iS), RT.Em.all(:,iS), nBins,0);
                efficacy_pureNback.Em.RTbins(:,iS) = do_bin2(efficacy_pureNback.Em.all(:,iS), RT.Em.all(:,iS), nBins,0);
        end
        %% extract info split by effort level
        for iEch = Ech_levels
            Ech_nm = ['Ech',num2str(iEch)];
            curr_Ech_trials = Ech.(task_id).all(:,iS) == iEch;
            switch task_id
                case 'Ep'
                    [perf.Ep.(Ech_nm).RTbins(:,iS), RT.Ep.(Ech_nm).RTbins(:,iS)] = do_bin2(perf.Ep.all(curr_Ech_trials,iS), RT.Ep.all(curr_Ech_trials,iS), nBins,0);
                    AUC.Ep.(Ech_nm).RTbins(:,iS) = do_bin2(AUC.Ep.all(curr_Ech_trials,iS), RT.Ep.all(curr_Ech_trials,iS), nBins,0);
                    peakF.Ep.(Ech_nm).RTbins(:,iS) = do_bin2(peakF.Ep.all(curr_Ech_trials,iS), RT.Ep.all(curr_Ech_trials,iS), nBins,0);
                    latency.Ep.(Ech_nm).RTbins(:,iS) = do_bin2(latency.Ep.all(curr_Ech_trials,iS), RT.Ep.all(curr_Ech_trials,iS), nBins,0);
                    AUC_overshoot.Ep.(Ech_nm).RTbins(:,iS) = do_bin2(AUC_overshoot.Ep.all(curr_Ech_trials,iS), RT.Ep.all(curr_Ech_trials,iS), nBins,0);
                case 'Em'
                    [perf.Em.(Ech_nm).RTbins(:,iS), RT.Em.(Ech_nm).RTbins(:,iS)] = do_bin2(perf.Em.all(curr_Ech_trials,iS), RT.Em.all(curr_Ech_trials,iS), nBins,0);
                    nCorrect.Em.(Ech_nm).RTbins(:,iS) = do_bin2(nCorrect.Em.all(curr_Ech_trials,iS), RT.Em.all(curr_Ech_trials,iS), nBins,0);
                    nErrors.Em.(Ech_nm).RTbins(:,iS) = do_bin2(nErrors.Em.all(curr_Ech_trials,iS), RT.Em.all(curr_Ech_trials,iS), nBins,0);
                    totalTime.Em.(Ech_nm).RTbins(:,iS) = do_bin2(totalTime.Em.all(curr_Ech_trials,iS), RT.Em.all(curr_Ech_trials,iS), nBins,0);
                    efficacy_with2first.Em.(Ech_nm).RTbins(:,iS) = do_bin2(efficacy_with2first.Em.all(curr_Ech_trials,iS), RT.Em.all(curr_Ech_trials,iS), nBins,0);
                    efficacy_pureNback.Em.(Ech_nm).RTbins(:,iS) = do_bin2(efficacy_pureNback.Em.all(curr_Ech_trials,iS), RT.Em.all(curr_Ech_trials,iS), nBins,0);
            end

            %% extract info split by money
            for iMoney = absMoney_levels
                money_nm_init = ['M',num2str(iMoney)];
                curr_money_level = money_level_ch.(task_id).all(:,iS) == iMoney;
                % reward trials
                curr_R_trials = RP.(task_id).all(:,iS) == 1;
                R_lvl_idx = (curr_R_trials.*curr_money_level.*curr_Ech_trials) == 1;

                % punishment trials
                curr_P_trials = RP.(task_id).all(:,iS) == 0;
                P_lvl_idx = (curr_P_trials.*curr_money_level.*curr_Ech_trials) == 1;
            end
        end % effort level
    end % physical/mental loop
    disp(['Subject ',num2str(iS),'/',num2str(NS),' - done']);
end % subject loop

%% average
[avg.Ep.perf,...
    sem.Ep.perf] = mean_sem_sd(perf.Ep.RTbins,2);
[avg.Ep.AUC,...
    sem.Ep.AUC] = mean_sem_sd(AUC.Ep.RTbins,2);
[avg.Ep.AUCovershoot,...
    sem.Ep.AUCovershoot] = mean_sem_sd(AUC_overshoot.Ep.RTbins,2);
[avg.Ep.peakF,...
    sem.Ep.peakF] = mean_sem_sd(peakF.Ep.RTbins,2);
[avg.Ep.latency,...
    sem.Ep.latency] = mean_sem_sd(latency.Ep.RTbins,2);
[avg.Ep_RT,...
    sem.Ep_RT] = mean_sem_sd(RT.Ep.RTbins,2);
[avg.Em.perf,...
    sem.Em.perf] = mean_sem_sd(perf.Em.RTbins,2);
[avg.Em.nCorrect,...
    sem.Em.nCorrect] = mean_sem_sd(nCorrect.Em.RTbins,2);
[avg.Em.nErrors,...
    sem.Em.nErrors] = mean_sem_sd(nErrors.Em.RTbins,2);
[avg.Em.totalTime,...
    sem.Em.totalTime] = mean_sem_sd(totalTime.Em.RTbins,2);
[avg.Em.efficacy_with2first,...
    sem.Em.efficacy_with2first] = mean_sem_sd(efficacy_with2first.Em.RTbins,2);
[avg.Em.efficacy_pureNback,...
    sem.Em.efficacy_pureNback] = mean_sem_sd(efficacy_pureNback.Em.RTbins,2);
[avg.Em_RT,...
    sem.Em_RT] = mean_sem_sd(RT.Em.RTbins,2);
% average per effort chosen
for iEch = Ech_levels
    Ech_nm = ['Ech',num2str(iEch)];
    curr_Ech_trials = Ech.(task_id).all(:,iS) == iEch;
    [avg.(['Ep_',Ech_nm]).perf,...
        sem.(['Ep_',Ech_nm]).perf] = mean_sem_sd(perf.Ep.(Ech_nm).RTbins,2);
    [avg.(['Ep_',Ech_nm]).AUC,...
        sem.(['Ep_',Ech_nm]).AUC] = mean_sem_sd(AUC.Ep.(Ech_nm).RTbins,2);
    [avg.(['Ep_',Ech_nm]).AUCovershoot,...
        sem.(['Ep_',Ech_nm]).AUCovershoot] = mean_sem_sd(AUC_overshoot.Ep.(Ech_nm).RTbins,2);
    [avg.(['Ep_',Ech_nm]).peakF,...
        sem.(['Ep_',Ech_nm]).peakF] = mean_sem_sd(peakF.Ep.(Ech_nm).RTbins,2);
    [avg.(['Ep_',Ech_nm]).latency,...
        sem.(['Ep_',Ech_nm]).latency] = mean_sem_sd(latency.Ep.(Ech_nm).RTbins,2);
    [avg.(['Ep_',Ech_nm,'_RT']),...
        sem.(['Ep_',Ech_nm,'_RT'])] = mean_sem_sd(RT.Ep.(Ech_nm).RTbins,2);
    [avg.(['Em_',Ech_nm]).perf,...
        sem.(['Em_',Ech_nm]).perf] = mean_sem_sd(perf.Em.(Ech_nm).RTbins,2);
    [avg.(['Em_',Ech_nm]).nCorrect,...
        sem.(['Em_',Ech_nm]).nCorrect] = mean_sem_sd(nCorrect.Em.(Ech_nm).RTbins,2);
    [avg.(['Em_',Ech_nm]).nErrors,...
        sem.(['Em_',Ech_nm]).nErrors] = mean_sem_sd(nErrors.Em.(Ech_nm).RTbins,2);
    [avg.(['Em_',Ech_nm]).totalTime,...
        sem.(['Em_',Ech_nm]).totalTime] = mean_sem_sd(totalTime.Em.(Ech_nm).RTbins,2);
    [avg.(['Em_',Ech_nm]).efficacy_with2first,...
        sem.(['Em_',Ech_nm]).efficacy_with2first] = mean_sem_sd(efficacy_with2first.Em.(Ech_nm).RTbins,2);
    [avg.(['Em_',Ech_nm]).efficacy_pureNback,...
        sem.(['Em_',Ech_nm]).efficacy_pureNback] = mean_sem_sd(efficacy_pureNback.Em.(Ech_nm).RTbins,2);
    [avg.(['Em_',Ech_nm,'_RT']),...
        sem.(['Em_',Ech_nm,'_RT'])] = mean_sem_sd(RT.Em.(Ech_nm).RTbins,2);
end

%% figure
lWidth = 3;
pSize = 40;

%% physical task - different correlations
physical_var_names = fieldnames(avg.Ep);
nEp_vars = length(physical_var_names);
fig;

for iEpFig = 1:nEp_vars
    Ep_var_nm = physical_var_names{iEpFig};
    subplot(2,3,iEpFig);
    hold on;
    er_Ep_hdl = errorbar(avg.Ep_RT,...
        avg.Ep.(Ep_var_nm),...
        sem.Ep.(Ep_var_nm));
    er_Ep_hdl.LineWidth = lWidth;
    er_Ep_hdl.Color = 'k';
    xlabel('choice RT (s)');
    switch Ep_var_nm
        case 'perf'
            ylabel('Performance (%)');
        case 'AUC'
            ylabel('AUC');
        case 'AUCovershoot'
            ylabel('AUC overshoot');
        case 'peakF'
            ylabel('Peak performance (%)');
        case 'latency'
            ylabel('Latency (s)');
    end
    legend_size(pSize);
end % physical variable loop

%% mental task - different correlations
mental_var_names = fieldnames(avg.Em);
nEm_vars = length(mental_var_names);
fig;

for iEmFig = 1:nEm_vars
    Em_var_nm = mental_var_names{iEmFig};
    subplot(2,3,iEmFig);
    hold on;
    er_Em_hdl = errorbar(avg.Em_RT,...
        avg.Em.(Em_var_nm),...
        sem.Em.(Em_var_nm));
    er_Em_hdl.LineWidth = lWidth;
    er_Em_hdl.Color = 'k';
    xlabel('choice RT (s)');
    switch Em_var_nm
        case 'perf'
            ylabel('Performance (%)');
        case 'nCorrect'
            ylabel({'Number correct';'answers'});
        case 'nErrors'
            ylabel('Number errors');
        case 'totalTime'
            ylabel('total time (s)');
        case 'efficacy_with2first'
            ylabel({'Efficacy';'(including 2 first replies)'});
        case 'efficacy_pureNback'
            ylabel('Efficacy');
    end
    legend_size(pSize);
end % mental variable loop

%% Effort chosen
%% physical task - different correlations
fig;

for iEpFig = 1:nEp_vars
    Ep_var_nm = physical_var_names{iEpFig};
    subplot(2,3,iEpFig);
    hold on;
    legend_hdl = NaN(1,length(Ech_levels));
    legend_nm = cell(1,length(Ech_levels));
    for iEch = Ech_levels
        Ech_nm = ['Ech',num2str(iEch)];
        legend_nm{iEch+1} = Ech_nm;
        curr_Ech_trials = Ech.(task_id).all(:,iS) == iEch;
        er_Ep_hdl = errorbar(avg.(['Ep_',Ech_nm,'_RT']),...
            avg.(['Ep_',Ech_nm]).(Ep_var_nm),...
            sem.(['Ep_',Ech_nm]).(Ep_var_nm));
        legend_hdl(iEch+1) = er_Ep_hdl;
        er_Ep_hdl.LineWidth = lWidth;
        er_Ep_hdl.Color = [0.2*iEch 1-0.2*iEch 0.2*iEch];
    end
    legend(legend_hdl, legend_nm);
    legend('boxoff');
    xlabel('choice RT (s)');
    switch Ep_var_nm
        case 'perf'
            ylabel('Performance (%)');
        case 'AUC'
            ylabel('AUC');
        case 'AUCovershoot'
            ylabel('AUC overshoot');
        case 'peakF'
            ylabel('Peak performance (%)');
        case 'latency'
            ylabel('Latency (s)');
    end
    legend_size(pSize);
end % physical variable loop

%% mental task - different correlations
fig;

for iEmFig = 1:nEm_vars
    Em_var_nm = mental_var_names{iEmFig};
    subplot(2,3,iEmFig);
    hold on;
    legend_hdl = NaN(1,length(Ech_levels));
    legend_nm = cell(1,length(Ech_levels));
    for iEch = Ech_levels
        Ech_nm = ['Ech',num2str(iEch)];
        legend_nm{iEch+1} = Ech_nm;
        er_Em_hdl = errorbar(avg.(['Em_',Ech_nm,'_RT']),...
            avg.(['Em_',Ech_nm]).(Em_var_nm),...
            sem.(['Em_',Ech_nm]).(Em_var_nm));
        legend_hdl(iEch+1) = er_Em_hdl;
        er_Em_hdl.LineWidth = lWidth;
        er_Em_hdl.Color = [0.2*iEch 1-0.2*iEch 0.2*iEch];
    end
    legend(legend_hdl, legend_nm);
    legend('boxoff');
    xlabel('choice RT (s)');
    switch Em_var_nm
        case 'perf'
            ylabel('Performance (%)');
        case 'nCorrect'
            ylabel({'Number correct';'answers'});
        case 'nErrors'
            ylabel('Number errors');
        case 'totalTime'
            ylabel('total time (s)');
        case 'efficacy_with2first'
            ylabel({'Efficacy';'(including 2 first replies)'});
        case 'efficacy_pureNback'
            ylabel('Efficacy');
    end
    legend_size(pSize);
end % mental variable loop