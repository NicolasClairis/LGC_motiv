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
    timePeriod_nm] = extract_ROI_betas_onsets_only_bis(computerRoot,...
    study_nm, subject_id, condition);

%% general parameters
nRuns = 4;
nTrialsPerRun = 54;
nTrials = nTrialsPerRun*nRuns;
nRunsPerTask = 2;
nTrialsPerTask = nTrialsPerRun.*nRunsPerTask;
task_names = {'physical','mental'};
task_ids = {'Ep','Em'};
nTasks = length(task_names);
[ROI_activity.allTrials,...
    perf.latency.allTrials,...
    perf.AUC.allTrials,...
    perf.forcePeak.allTrials,...
    perf.AUC_overshoot.allTrials,...
    perf.forcePeak_N.allTrials,...
    perf.n_errors.allTrials,...
    perf.successSpeed.allTrials,...
    choice_hE,...
    hE_level,...
    E_chosen,...
    task_cstt.physical, task_cstt.mental,...
    r1_cstt.physical, r2_cstt.physical,...
    r1_cstt.mental, r2_cstt.mental,...
    perf_fit.AUC.allTrials,...
    perf_fit.forcePeak.allTrials,...
    perf_fit.AUC_overshoot.allTrials,...
    perf_fit.forcePeak_N.allTrials,...
    perf_fit.n_errors.allTrials,...
    perf_fit.successSpeed.allTrials,...
    perf_fit.latency.allTrials,...
    perf_fit.latency.allTrials,...
    perf_fit.perf.allTrials,...
    perf_fit.perf.allTrials] = deal(NaN(nTrials, NS));
n_hE_levels = 3;
hE_levels = 1:3;
n_Ech = 4;
Ech_levels = 0:3;
Ep_fields = {'latency','AUC','forcePeak','forcePeak_N','AUC_overshoot','perf'};
n_Ep_fields = length(Ep_fields);
Em_fields = {'n_errors','successSpeed','latency','perf'};
n_Em_fields = length(Em_fields);
for iTask = 1:nTasks
    task_nm = task_names{iTask};
    task_id = task_ids{iTask};
    ROI_activity.(task_nm).perEch= deal(NaN(n_Ech,NS));
    [ROI_activity.(task_nm).choice_highE.perHighElevel,...
        ROI_activity.(task_nm).choice_lowE.perHighElevel] = deal(NaN(n_hE_levels,NS));
    switch task_id
        case 'Ep'
            for iEp_field1 = 1:n_Ep_fields
                Ep_field_nm = Ep_fields{iEp_field1};
                % data split for E chosen*ROI median split
                [perf.(task_nm).(Ep_field_nm).perEch.(['low_',ROI_short_nm]),...
                    perf.(task_nm).(Ep_field_nm).perEch.(['high_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Ep_field_nm).perEch.(['low_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Ep_field_nm).perEch.(['high_',ROI_short_nm])] = deal(NaN(n_Ech,NS));
                % high E*ROI median split*choice
                [perf.(task_nm).(Ep_field_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
                    perf.(task_nm).(Ep_field_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Ep_field_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Ep_field_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
                    perf.(task_nm).(Ep_field_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
                    perf.(task_nm).(Ep_field_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Ep_field_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Ep_field_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm])] = deal(NaN(n_hE_levels,NS));
                % betas for tests
                [betas.(task_nm).(Ep_field_nm).bROI.allS,...
                    betas.(task_nm).(Ep_field_nm).bE.allS,...
                    betas.(task_nm).(Ep_field_nm).bROI_x_E.allS] = deal(NaN(1,NS));
            end % field
        case 'Em'
            for iEm_field1 = 1:n_Em_fields
                Em_field_nm = Em_fields{iEm_field1};
                % data split for E chosen*ROI median split
                [perf.(task_nm).(Em_field_nm).perEch.(['low_',ROI_short_nm]),...
                    perf.(task_nm).(Em_field_nm).perEch.(['high_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Em_field_nm).perEch.(['low_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Em_field_nm).perEch.(['high_',ROI_short_nm])] = deal(NaN(n_Ech,NS));
                % high E*ROI median split*choice
                [perf.(task_nm).(Em_field_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
                    perf.(task_nm).(Em_field_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Em_field_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Em_field_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
                    perf.(task_nm).(Em_field_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
                    perf.(task_nm).(Em_field_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Em_field_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
                    perf_fit.(task_nm).(Em_field_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm])] = deal(NaN(n_hE_levels,NS));
                % betas for tests
                [betas.(task_nm).(Em_field_nm).bROI.allS,...
                    betas.(task_nm).(Em_field_nm).bE.allS,...
                    betas.(task_nm).(Em_field_nm).bROI_x_E.allS] = deal(NaN(1,NS));
            end % field
    end % task
end % task loop

%% load performance and extract corresponding ROI activity
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
    % extract runs
    [runsStruct] = runs_definition(study_nm, sub_nm, condition);
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    for iRun = 1:length(okRuns)
        kRun = okRuns(iRun);
        run_nm = num2str(kRun);
        task_nm_tmp = taskNames{iRun};
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
                r1_cstt.(task_fullName)(runTrials_idx, iS) = 1;
                r2_cstt.(task_fullName)(runTrials_idx, iS) = 0;
            case {3,4}
                run_nm_bis = ['run',num2str(2)];
                r1_cstt.(task_fullName)(runTrials_idx, iS) = 0;
                r2_cstt.(task_fullName)(runTrials_idx, iS) = 1;
        end
        
        %% load choice informations
        choice_hE(runTrials_idx, iS) = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        hE_level(runTrials_idx, iS) = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        E_chosen(runTrials_idx, iS) = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        perf.perf.allTrials(runTrials_idx,iS) = extract_perf(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        task_cstt.(task_fullName)(runTrials_idx,iS) = 1;
        %% load the performance data
        switch task_nm_tmp
            case 'Ep'
                [latency_Ep, AUC, forcePeak, AUC_overshoot, ~, forcePeak_N, ~] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
                perf.latency.allTrials(runTrials_idx,iS) = latency_Ep.allTrials;
                perf.AUC.allTrials(runTrials_idx,iS) = AUC.allTrials;
                perf.forcePeak.allTrials(runTrials_idx,iS) = forcePeak.allTrials;
                perf.forcePeak_N.allTrials(runTrials_idx,iS) = forcePeak_N.allTrials;
                perf.AUC_overshoot.allTrials(runTrials_idx,iS) = AUC_overshoot.allTrials;
                
            case 'Em'
                [successSpeed, ~, n_errors, RT_avg,...
                    ~,~,~,~,~,~,latency_Em] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
                perf.successSpeed.allTrials(runTrials_idx,iS) = successSpeed.allTrials;
                perf.n_errors.allTrials(runTrials_idx,iS) = n_errors.allTrials;
                perf.RT_avg.allTrials(runTrials_idx,iS) = RT_avg.allTrials;
                perf.latency.allTrials(runTrials_idx,iS) = latency_Em.allTrials;
        end
        
        %% extract fMRI ROI mediator
        ROI_activity.allTrials(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
    end % run loop
    
    %% loop over tasks
    for iTask = 1:nTasks
        task_nm = task_names{iTask};
        task_idx = task_cstt.(task_nm)(:,iS);
        
        %% statistical test and fit
        % general parameters
        task_trials = task_cstt.(task_nm)(:,iS) == 1;
        r1_task_cstt = r1_cstt.(task_nm)(task_trials, iS);
        r2_task_cstt = r2_cstt.(task_nm)(task_trials, iS);
        switch task_nm
            case 'physical'
                % general parameters
                E_chosen_Ep = E_chosen(task_trials,iS);
                ROI_Ep = ROI_activity.allTrials(task_trials,iS);
                ROI_x_Ech_Ep = ROI_Ep.*E_chosen_Ep;
                x_Ep = [E_chosen_Ep, ROI_Ep, ROI_x_Ech_Ep];
                x_reg_Ep_names = {'bE','bROI','bROI_x_E'};
                
                % Ep: AUC
                [betas_Ep_AUC_tmp, perf_fit.AUC.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Ep, perf.AUC.allTrials(task_trials,iS), x_reg_Ep_names);
                for iR = 1:length(x_reg_Ep_names)
                    b_nm = x_reg_Ep_names{iR};
                    betas.Ep.AUC.(b_nm).allS(iS) = betas_Ep_AUC_tmp.(b_nm);
                end
                
                % Ep: force peak
                [betas_Ep_forcePeak_tmp, perf_fit.forcePeak.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Ep, perf.forcePeak.allTrials(task_trials,iS), x_reg_Ep_names);
                for iR = 1:length(x_reg_Ep_names)
                    b_nm = x_reg_Ep_names{iR};
                    betas.Ep.forcePeak.(b_nm).allS(iS) = betas_Ep_forcePeak_tmp.(b_nm);
                end
                
                % Ep: AUC overshoot
                [betas_Ep_AUC_overshoot_tmp, perf_fit.AUC_overshoot.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Ep, perf.AUC_overshoot.allTrials(task_trials,iS), x_reg_Ep_names);
                for iR = 1:length(x_reg_Ep_names)
                    b_nm = x_reg_Ep_names{iR};
                    betas.Ep.AUC_overshoot.(b_nm).allS(iS) = betas_Ep_AUC_overshoot_tmp.(b_nm);
                end
                
                % Ep: force peak in Newtons
                [betas_Ep_forcePeak_N_tmp, perf_fit.forcePeak_N.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Ep, perf.forcePeak_N.allTrials(task_trials,iS), x_reg_Ep_names);
                for iR = 1:length(x_reg_Ep_names)
                    b_nm = x_reg_Ep_names{iR};
                    betas.Ep.forcePeak_N.(b_nm).allS(iS) = betas_Ep_forcePeak_N_tmp.(b_nm);
                end
                
                % Ep latency
                [betas_Ep_latency_tmp, perf_fit.latency.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Ep, perf.latency.allTrials(task_trials,iS), x_reg_Ep_names);
                for iR = 1:length(x_reg_Ep_names)
                    b_nm = x_reg_Ep_names{iR};
                    betas.Ep.latency.(b_nm).allS(iS) = betas_Ep_latency_tmp.(b_nm);
                end
                
                % Ep performance
                [betas_Ep_perf_tmp, perf_fit.perf.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Ep, perf.perf.allTrials(task_trials,iS), x_reg_Ep_names);
                for iR = 1:length(x_reg_Ep_names)
                    b_nm = x_reg_Ep_names{iR};
                    betas.Ep.perf.(b_nm).allS(iS) = betas_Ep_perf_tmp.(b_nm);
                end
            case 'mental'
                % general parameters
                E_chosen_Em = E_chosen(task_trials,iS);
                ROI_Em = ROI_activity.allTrials(task_trials,iS);
                ROI_x_Ech_Em = ROI_Em.*E_chosen_Em;
                x_Em = [E_chosen_Em, ROI_Em, ROI_x_Ech_Em];
                x_reg_Em_names = {'bE','bROI','bROI_x_E'};
                
                % Em: number of errors
                [betas_Em_nErrors_tmp, perf_fit.n_errors.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Em, perf.n_errors.allTrials(task_trials,iS), x_reg_Em_names);
                for iR = 1:length(x_reg_Em_names)
                    b_nm = x_reg_Em_names{iR};
                    betas.Em.n_errors.(b_nm).allS(iS) = betas_Em_nErrors_tmp.(b_nm);
                end
                
                % Em: success speed
                [betas_Em_successSpeed_tmp, perf_fit.successSpeed.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Em, perf.successSpeed.allTrials(task_trials,iS), x_reg_Em_names);
                for iR = 1:length(x_reg_Em_names)
                    b_nm = x_reg_Em_names{iR};
                    betas.Em.successSpeed.(b_nm).allS(iS) = betas_Em_successSpeed_tmp.(b_nm);
                end
                
                % Em latency
                [betas_Em_latency_tmp, perf_fit.latency.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Em, perf.latency.allTrials(task_trials,iS), x_reg_Em_names);
                for iR = 1:length(x_reg_Em_names)
                    b_nm = x_reg_Em_names{iR};
                    betas.Em.latency.(b_nm).allS(iS) = betas_Em_latency_tmp.(b_nm);
                end
                
                % Em performance
                [betas_Em_perf_tmp, perf_fit.perf.allTrials(task_trials,iS)] = glmfit_adapted_for_runConstants(r1_task_cstt, r2_task_cstt,...
                    x_Em, perf.perf.allTrials(task_trials,iS), x_reg_Em_names);
                for iR = 1:length(x_reg_Em_names)
                    b_nm = x_reg_Em_names{iR};
                    betas.Em.perf.(b_nm).allS(iS) = betas_Em_perf_tmp.(b_nm);
                end
        end
        
        %% split data by effort levels and effort chosen and ROI activity
        for iEch = 1:n_Ech
            jEch_idx = (E_chosen(:,iS) == (iEch - 1)).*(task_idx == 1) == 1;
            ROI_activity.(task_nm).perEch(iEch, iS)       = mean(ROI_activity.allTrials(jEch_idx,iS),1,'omitnan');
            median_ROI_Ech = median(ROI_activity.allTrials(jEch_idx, iS),1,'omitnan');
            low_ROI_Ech = ((ROI_activity.allTrials(:,iS) <= median_ROI_Ech).*jEch_idx) == 1;
            high_ROI_Ech = ((ROI_activity.allTrials(:,iS) > median_ROI_Ech).*jEch_idx) == 1;
            switch task_nm
                case 'physical'
                    for iEp_field2 = 1:n_Ep_fields
                        Ep_field_nm2 = Ep_fields{iEp_field2};
                        % low ROI level
                        perf.(task_nm).(Ep_field_nm2).perEch.(['low_',ROI_short_nm])(iEch, iS)       = mean(perf.(Ep_field_nm2).allTrials(low_ROI_Ech,iS),1,'omitnan');
                        perf_fit.(task_nm).(Ep_field_nm2).perEch.(['low_',ROI_short_nm])(iEch, iS)   = mean(perf_fit.(Ep_field_nm2).allTrials(low_ROI_Ech,iS),1,'omitnan');
                        % high ROI level
                        perf.(task_nm).(Ep_field_nm2).perEch.(['high_',ROI_short_nm])(iEch, iS)       = mean(perf.(Ep_field_nm2).allTrials(high_ROI_Ech,iS),1,'omitnan');
                        perf_fit.(task_nm).(Ep_field_nm2).perEch.(['high_',ROI_short_nm])(iEch, iS)   = mean(perf_fit.(Ep_field_nm2).allTrials(high_ROI_Ech,iS),1,'omitnan');
                    end
                case 'mental'
                    for iEm_field2 = 1:n_Em_fields
                        Em_field_nm2 = Em_fields{iEm_field2};
                        % low ROI level
                        perf.(task_nm).(Em_field_nm2).perEch.(['low_',ROI_short_nm])(iEch, iS)      = mean(perf.(Em_field_nm2).allTrials(low_ROI_Ech,iS),1,'omitnan');
                        perf_fit.(task_nm).(Em_field_nm2).perEch.(['low_',ROI_short_nm])(iEch, iS)  = mean(perf_fit.(Em_field_nm2).allTrials(low_ROI_Ech,iS),1,'omitnan');
                        % high ROI level
                        perf.(task_nm).(Em_field_nm2).perEch.(['high_',ROI_short_nm])(iEch, iS)     = mean(perf.(Em_field_nm2).allTrials(high_ROI_Ech,iS),1,'omitnan');
                        perf_fit.(task_nm).(Em_field_nm2).perEch.(['high_',ROI_short_nm])(iEch, iS) = mean(perf_fit.(Em_field_nm2).allTrials(high_ROI_Ech,iS),1,'omitnan');
                    end
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
                    for iEp_field3 = 1:n_Ep_fields
                        Ep_field_nm3 = Ep_fields{iEp_field3};
                        % low ROI - low effort chosen
                        perf.(task_nm).(Ep_field_nm3).choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)       = mean(perf.(Ep_field_nm3).allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                        perf_fit.(task_nm).(Ep_field_nm3).choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)       = mean(perf_fit.(Ep_field_nm3).allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                        % high ROI - low effort chosen
                        perf.(task_nm).(Ep_field_nm3).choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)       = mean(perf.(Ep_field_nm3).allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                        perf_fit.(task_nm).(Ep_field_nm3).choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)       = mean(perf_fit.(Ep_field_nm3).allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                        % low ROI - high effort chosen
                        perf.(task_nm).(Ep_field_nm3).choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)      = mean(perf.(Ep_field_nm3).allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                        perf_fit.(task_nm).(Ep_field_nm3).choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)      = mean(perf_fit.(Ep_field_nm3).allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                        % high ROI - high effort chosen
                        perf.(task_nm).(Ep_field_nm3).choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)      = mean(perf.(Ep_field_nm3).allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                        perf_fit.(task_nm).(Ep_field_nm3).choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)      = mean(perf_fit.(Ep_field_nm3).allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                    end
                case 'mental'
                    for iEm_field3 = 1:n_Em_fields
                        Em_field_nm3 = Em_fields{iEm_field3};
                        % low ROI - low effort chosen
                        perf.(task_nm).(Em_field_nm3).choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)         = mean(perf.(Em_field_nm3).allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                        perf_fit.(task_nm).(Em_field_nm3).choice_lowE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)         = mean(perf_fit.(Em_field_nm3).allTrials(low_ROI_hE_lowCh,iS),1,'omitnan');
                        % high ROI - low effort chosen
                        perf.(task_nm).(Em_field_nm3).choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)        = mean(perf.(Em_field_nm3).allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                        perf_fit.(task_nm).(Em_field_nm3).choice_lowE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)        = mean(perf_fit.(Em_field_nm3).allTrials(high_ROI_hE_lowCh,iS),1,'omitnan');
                        % low ROI - high effort chosen
                        perf.(task_nm).(Em_field_nm3).choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)        = mean(perf.(Em_field_nm3).allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                        perf_fit.(task_nm).(Em_field_nm3).choice_highE.perHighElevel.(['low_',ROI_short_nm])(iEff, iS)        = mean(perf_fit.(Em_field_nm3).allTrials(low_ROI_hE_highCh,iS),1,'omitnan');
                        % high ROI - high effort chosen
                        perf.(task_nm).(Em_field_nm3).choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)       = mean(perf.(Em_field_nm3).allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                        perf_fit.(task_nm).(Em_field_nm3).choice_highE.perHighElevel.(['high_',ROI_short_nm])(iEff, iS)       = mean(perf_fit.(Em_field_nm3).allTrials(high_ROI_hE_highCh,iS),1,'omitnan');
                    end
            end
            
        end % effort level
        
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
     
     %% check if betas are significant at the group level
     % physical variables
     Ep_fields = fieldnames(betas.Ep);
     for iField_Ep = 1:length(Ep_fields)
         Ep_field_nm = Ep_fields{iField_Ep};
         % bE
         [betas.Ep.(Ep_field_nm).bE.mean,...
             betas.Ep.(Ep_field_nm).bE.sem] = mean_sem_sd(betas.Ep.(Ep_field_nm).bE.allS,2);
         [~,pval.Ep.(Ep_field_nm).bE] = ttest(betas.Ep.(Ep_field_nm).bE.allS);
         
         % bROI
         [betas.Ep.(Ep_field_nm).bROI.mean,...
             betas.Ep.(Ep_field_nm).bROI.sem] = mean_sem_sd(betas.Ep.(Ep_field_nm).bROI.allS,2);
         [~,pval.Ep.(Ep_field_nm).bROI] = ttest(betas.Ep.(Ep_field_nm).bROI.allS);
         
         % bROI*E
         [betas.Ep.(Ep_field_nm).bROI_x_E.mean,...
             betas.Ep.(Ep_field_nm).bROI_x_E.sem] = mean_sem_sd(betas.Ep.(Ep_field_nm).bROI_x_E.allS,2);
         [~,pval.Ep.(Ep_field_nm).bROI_x_E] = ttest(betas.Ep.(Ep_field_nm).bROI_x_E.allS);
     end % Ep field loop
     
     % mental variables
     Em_fields = fieldnames(betas.Em);
     for iField_Em = 1:length(Em_fields)
         Em_field_nm = Em_fields{iField_Em};
         % bE
         [betas.Em.(Em_field_nm).bE.mean,...
             betas.Em.(Em_field_nm).bE.sem] = mean_sem_sd(betas.Em.(Em_field_nm).bE.allS,2);
         [~,pval.Em.(Em_field_nm).bE] = ttest(betas.Em.(Em_field_nm).bE.allS);
         
         % bROI
         [betas.Em.(Em_field_nm).bROI.mean,...
             betas.Em.(Em_field_nm).bROI.sem] = mean_sem_sd(betas.Em.(Em_field_nm).bROI.allS,2);
         [~,pval.Em.(Em_field_nm).bROI] = ttest(betas.Em.(Em_field_nm).bROI.allS);
         
         % bROI*E
         [betas.Em.(Em_field_nm).bROI_x_E.mean,...
             betas.Em.(Em_field_nm).bROI_x_E.sem] = mean_sem_sd(betas.Em.(Em_field_nm).bROI_x_E.allS,2);
         [~,pval.Em.(Em_field_nm).bROI_x_E] = ttest(betas.Em.(Em_field_nm).bROI_x_E.allS);
     end % Ep field loop
end % task loop
     
%% loop through levels of activity of the ROI
ROI_level = {'low','high'};
for iROI_lvl = 1:length(ROI_level)
    ROI_info_nm = [ROI_level{iROI_lvl},'_',ROI_short_nm];
    
    %% physical
    for iEp_field4 = 1:n_Ep_fields
        Ep_field_nm4 = Ep_fields{iEp_field4};
        % E chosen
        [m_perf.physical.(Ep_field_nm4).perEch.(ROI_info_nm),...
            sem_perf.physical.(Ep_field_nm4).perEch.(ROI_info_nm)] = mean_sem_sd(perf.physical.(Ep_field_nm4).perEch.(ROI_info_nm), 2);
        [m_perf_fit.physical.(Ep_field_nm4).perEch.(ROI_info_nm),...
            sem_perf_fit.physical.(Ep_field_nm4).perEch.(ROI_info_nm)] = mean_sem_sd(perf_fit.physical.(Ep_field_nm4).perEch.(ROI_info_nm), 2);
        % low effort chosen
        [m_perf.physical.(Ep_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm),...
            sem_perf.physical.(Ep_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.(Ep_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm), 2);
        [m_perf_fit.physical.(Ep_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm),...
            sem_perf_fit.physical.(Ep_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf_fit.physical.(Ep_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm), 2);
        % high effort chosen
        [m_perf.physical.(Ep_field_nm4).choice_highE.perHighElevel.(ROI_info_nm),...
            sem_perf.physical.(Ep_field_nm4).choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.physical.(Ep_field_nm4).choice_highE.perHighElevel.(ROI_info_nm), 2);
        [m_perf_fit.physical.(Ep_field_nm4).choice_highE.perHighElevel.(ROI_info_nm),...
            sem_perf_fit.physical.(Ep_field_nm4).choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf_fit.physical.(Ep_field_nm4).choice_highE.perHighElevel.(ROI_info_nm), 2);
    end
    %% mental
    for iEm_field4 = 1:n_Em_fields
        Em_field_nm4 = Em_fields{iEm_field4};
        % E chosen
        [m_perf.mental.(Em_field_nm4).perEch.(ROI_info_nm),...
            sem_perf.mental.(Em_field_nm4).perEch.(ROI_info_nm)] = mean_sem_sd(perf.mental.(Em_field_nm4).perEch.(ROI_info_nm), 2);
        [m_perf_fit.mental.(Em_field_nm4).perEch.(ROI_info_nm),...
            sem_perf_fit.mental.(Em_field_nm4).perEch.(ROI_info_nm)] = mean_sem_sd(perf_fit.mental.(Em_field_nm4).perEch.(ROI_info_nm), 2);
        % low effort chosen
        [m_perf.mental.(Em_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm),...
            sem_perf.mental.(Em_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.mental.(Em_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm), 2);
        [m_perf_fit.mental.(Em_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm),...
            sem_perf_fit.mental.(Em_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf_fit.mental.(Em_field_nm4).choice_lowE.perHighElevel.(ROI_info_nm), 2);
        % high effort chosen
        [m_perf.mental.(Em_field_nm4).choice_highE.perHighElevel.(ROI_info_nm),...
            sem_perf.mental.(Em_field_nm4).choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf.mental.(Em_field_nm4).choice_highE.perHighElevel.(ROI_info_nm), 2);
        [m_perf_fit.mental.(Em_field_nm4).choice_highE.perHighElevel.(ROI_info_nm),...
            sem_perf_fit.mental.(Em_field_nm4).choice_highE.perHighElevel.(ROI_info_nm)] = mean_sem_sd(perf_fit.mental.(Em_field_nm4).choice_highE.perHighElevel.(ROI_info_nm), 2);
    end
end % ROI activity

%% display figures
pSize = 30;
lWidth = 1;
blue = [0 143 255]./255;
purple = [128 0 128]./255;
%% E chosen - grip
fig;

Ep_feature = {'latency','forcePeak','forcePeak_N',...
    'AUC','AUC_overshoot'};
for iEp = 1:length(Ep_feature)
    Ep_nm = Ep_feature{iEp};
    
    subplot(2,3,iEp);
    hold on;
    Ech_low_ROI_hdl = errorbar(Ech_levels,...
        m_perf.physical.(Ep_nm).perEch.(['low_',ROI_short_nm]),...
        sem_perf.physical.(Ep_nm).perEch.(['low_',ROI_short_nm]));
    Ech_high_ROI_hdl = errorbar(Ech_levels,...
        m_perf.physical.(Ep_nm).perEch.(['high_',ROI_short_nm]),...
        sem_perf.physical.(Ep_nm).perEch.(['high_',ROI_short_nm]));
    % add the fit
    Ech_low_ROI_fit_hdl = plot(Ech_levels,...
        m_perf_fit.physical.(Ep_nm).perEch.(['low_',ROI_short_nm]));
    Ech_high_ROI_fit_hdl = plot(Ech_levels,...
        m_perf_fit.physical.(Ep_nm).perEch.(['high_',ROI_short_nm]));
    % modify curve properties
    Ech_low_ROI_hdl.Marker = 'o';
    Ech_low_ROI_hdl.LineStyle = 'none';
    Ech_low_ROI_hdl.LineWidth = lWidth;
    Ech_low_ROI_hdl.Color = blue;
    Ech_high_ROI_hdl.Marker = 'o';
    Ech_high_ROI_hdl.LineStyle = 'none';
    Ech_high_ROI_hdl.LineWidth = lWidth;
    Ech_high_ROI_hdl.Color = purple;
    % fit handle
    Ech_low_ROI_fit_hdl.LineStyle = '--';
    Ech_low_ROI_fit_hdl.LineWidth = lWidth;
    Ech_low_ROI_fit_hdl.Color = blue;
    Ech_high_ROI_fit_hdl.LineStyle = '--';
    Ech_high_ROI_fit_hdl.LineWidth = lWidth;
    Ech_high_ROI_fit_hdl.Color = purple;
    legend([Ech_high_ROI_hdl, Ech_low_ROI_hdl],...
        {['high ',ROI_short_nm],['low ',ROI_short_nm]});
    legend('boxoff');
    legend('Location','NorthWest');
    xticks(Ech_levels);
    xlabel('E chosen');
    switch Ep_nm
        case 'latency'
            % latency = f(ROI)
            ylabel('Latency to squeeze (s)');
        case 'forcePeak'
            % force peak = f(ROI)
            ylabel('Peak force (%)');
        case 'forcePeak_N'
            % force peak (Newtons) = f(ROI)
            ylabel('Peak force (N)');
        case 'AUC'
            % AUC force = f(ROI)
            ylabel('AUC force');
        case 'AUC_overshoot'
            % AUC overshoot force = f(ROI)
            ylabel('AUC force overshoot');
    end
    legend_size(pSize);
end % Ep feature loop

%% E chosen - 2-back
fig;

Em_feature = {'successSpeed','n_errors','latency'};
for iEm = 1:length(Em_feature)
    Em_nm = Em_feature{iEm};
    
    subplot(2,2,iEm);
    hold on;
    Ech_low_ROI_hdl = errorbar(Ech_levels,...
        m_perf.mental.(Em_nm).perEch.(['low_',ROI_short_nm]),...
        sem_perf.mental.(Em_nm).perEch.(['low_',ROI_short_nm]));
    Ech_high_ROI_hdl = errorbar(Ech_levels,...
        m_perf.mental.(Em_nm).perEch.(['high_',ROI_short_nm]),...
        sem_perf.mental.(Em_nm).perEch.(['high_',ROI_short_nm]));
    % add the fit
    Ech_low_ROI_fit_hdl = plot(Ech_levels,...
        m_perf_fit.mental.(Em_nm).perEch.(['low_',ROI_short_nm]));
    Ech_high_ROI_fit_hdl = plot(Ech_levels,...
        m_perf_fit.mental.(Em_nm).perEch.(['high_',ROI_short_nm]));
    % modify curve properties
    Ech_low_ROI_hdl.Marker = 'o';
    Ech_low_ROI_hdl.LineStyle = 'none';
    Ech_low_ROI_hdl.LineWidth = lWidth;
    Ech_low_ROI_hdl.Color = blue;
    Ech_high_ROI_hdl.Marker = 'o';
    Ech_high_ROI_hdl.LineStyle = 'none';
    Ech_high_ROI_hdl.LineWidth = lWidth;
    Ech_high_ROI_hdl.Color = purple;
    % fit handle
    Ech_low_ROI_fit_hdl.LineStyle = '--';
    Ech_low_ROI_fit_hdl.LineWidth = lWidth;
    Ech_low_ROI_fit_hdl.Color = blue;
    Ech_high_ROI_fit_hdl.LineStyle = '--';
    Ech_high_ROI_fit_hdl.LineWidth = lWidth;
    Ech_high_ROI_fit_hdl.Color = purple;
    legend([Ech_high_ROI_hdl, Ech_low_ROI_hdl],...
        {['high ',ROI_short_nm],['low ',ROI_short_nm]});
    legend('boxoff');
    xticks(Ech_levels);
    xlabel('E chosen');
    switch Em_nm
        case 'successSpeed'
            % success speed = f(ROI)
            ylabel('success speed (s)');
        case 'n_errors'
            % number of errors made = f(ROI)
            ylabel('Number of errors');
        case 'latency'
            % latency = f(ROI)
            ylabel('latency (s)');
    end
    legend_size(pSize);
end % Em feature loop

%% effort level and choice - grip
fig;

for iEp = 1:length(Ep_feature)
    Ep_nm = Ep_feature{iEp};
    subplot(2,3,iEp);
    hold on;
    hE_lowEchoice_lowROI_hdl = errorbar(hE_levels,...
        m_perf.physical.(Ep_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
        sem_perf.physical.(Ep_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
    hE_lowEchoice_highROI_hdl = errorbar(hE_levels,...
        m_perf.physical.(Ep_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
        sem_perf.physical.(Ep_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
    hE_highEchoice_lowROI_hdl = errorbar(hE_levels,...
        m_perf.physical.(Ep_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
        sem_perf.physical.(Ep_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]));
    hE_highEchoice_highROI_hdl = errorbar(hE_levels,...
        m_perf.physical.(Ep_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
        sem_perf.physical.(Ep_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]));
    % add fit
    hE_lowEchoice_lowROI_fit_hdl = plot(hE_levels,...
        m_perf_fit.physical.(Ep_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
    hE_lowEchoice_highROI_fit_hdl = plot(hE_levels,...
        m_perf_fit.physical.(Ep_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
    hE_highEchoice_lowROI_fit_hdl = plot(hE_levels,...
        m_perf_fit.physical.(Ep_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]));
    hE_highEchoice_highROI_fit_hdl = plot(hE_levels,...
        m_perf_fit.physical.(Ep_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]));
    % modify curve properties
    hE_lowEchoice_lowROI_hdl.Marker = 'o';
    hE_lowEchoice_lowROI_hdl.LineStyle = 'none';
    hE_lowEchoice_lowROI_hdl.LineWidth = lWidth;
    hE_lowEchoice_lowROI_hdl.Color = blue;
    hE_lowEchoice_highROI_hdl.Marker = 'o';
    hE_lowEchoice_highROI_hdl.LineStyle = 'none';
    hE_lowEchoice_highROI_hdl.LineWidth = lWidth;
    hE_lowEchoice_highROI_hdl.Color = purple;
    hE_highEchoice_lowROI_hdl.Marker = 'o';
    hE_highEchoice_lowROI_hdl.LineStyle = 'none';
    hE_highEchoice_lowROI_hdl.LineWidth = lWidth;
    hE_highEchoice_lowROI_hdl.Color = blue;
    hE_highEchoice_highROI_hdl.Marker = 'o';
    hE_highEchoice_highROI_hdl.LineStyle = 'none';
    hE_highEchoice_highROI_hdl.LineWidth = lWidth;
    hE_highEchoice_highROI_hdl.Color = purple;
    % curve for fit
    % modify curve properties
    hE_lowEchoice_lowROI_fit_hdl.LineStyle = '--';
    hE_lowEchoice_lowROI_fit_hdl.LineWidth = lWidth;
    hE_lowEchoice_lowROI_fit_hdl.Color = blue;
    hE_lowEchoice_highROI_fit_hdl.LineStyle = '--';
    hE_lowEchoice_highROI_fit_hdl.LineWidth = lWidth;
    hE_lowEchoice_highROI_fit_hdl.Color = purple;
    hE_highEchoice_lowROI_fit_hdl.LineStyle = '-';
    hE_highEchoice_lowROI_fit_hdl.LineWidth = lWidth;
    hE_highEchoice_lowROI_fit_hdl.Color = blue;
    hE_highEchoice_highROI_fit_hdl.LineStyle = '-';
    hE_highEchoice_highROI_fit_hdl.LineWidth = lWidth;
    hE_highEchoice_highROI_fit_hdl.Color = purple;
%     legend([hE_highEchoice_highROI_hdl,...
%         hE_highEchoice_lowROI_hdl,...
%         hE_lowEchoice_highROI_hdl,...
%         hE_lowEchoice_lowROI_hdl],...
%         {['hEch - high ',ROI_short_nm],...
%         ['hEch - low ',ROI_short_nm],...
%         ['lowEch - high ',ROI_short_nm],...
%         ['lowEch - low ',ROI_short_nm]});
%     legend('boxoff');
    xticks(hE_levels);
    xlabel('high effort level');
    switch Ep_nm
        case 'latency'
            % latency = f(ROI)
            ylabel('Latency to squeeze (s)');
        case 'forcePeak'
            % force peak = f(ROI)
            ylabel('Peak force (%)');
        case 'forcePeak_N'
            % force peak (Newtons) = f(ROI)
            ylabel('Peak force (N)');
        case 'AUC'
            % AUC force = f(ROI)
            ylabel('AUC force');
        case 'AUC_overshoot'
            % AUC overshoot force = f(ROI)
            ylabel('AUC force overshoot');
    end
    legend_size(pSize);
end % Ep feature loop

%% effort level and choice - 2-back
fig;

for iEm = 1:length(Em_feature)
    Em_nm = Em_feature{iEm};
    
    subplot(2,2,iEm);
    hold on;
    hE_lowEchoice_lowROI_hdl = errorbar(hE_levels,...
        m_perf.mental.(Em_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]),...
        sem_perf.mental.(Em_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
    hE_lowEchoice_highROI_hdl = errorbar(hE_levels,...
        m_perf.mental.(Em_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm]),...
        sem_perf.mental.(Em_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
    hE_highEchoice_lowROI_hdl = errorbar(hE_levels,...
        m_perf.mental.(Em_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]),...
        sem_perf.mental.(Em_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]));
    hE_highEchoice_highROI_hdl = errorbar(hE_levels,...
        m_perf.mental.(Em_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]),...
        sem_perf.mental.(Em_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]));
    % add fit
    hE_lowEchoice_lowROI_fit_hdl = plot(hE_levels,...
        m_perf_fit.mental.(Em_nm).choice_lowE.perHighElevel.(['low_',ROI_short_nm]));
    hE_lowEchoice_highROI_fit_hdl = plot(hE_levels,...
        m_perf_fit.mental.(Em_nm).choice_lowE.perHighElevel.(['high_',ROI_short_nm]));
    hE_highEchoice_lowROI_fit_hdl = plot(hE_levels,...
        m_perf_fit.mental.(Em_nm).choice_highE.perHighElevel.(['low_',ROI_short_nm]));
    hE_highEchoice_highROI_fit_hdl = plot(hE_levels,...
        m_perf_fit.mental.(Em_nm).choice_highE.perHighElevel.(['high_',ROI_short_nm]));
    % modify curve properties
    hE_lowEchoice_lowROI_hdl.Marker = 'o';
    hE_lowEchoice_lowROI_hdl.LineStyle = 'none';
    hE_lowEchoice_lowROI_hdl.LineWidth = lWidth;
    hE_lowEchoice_lowROI_hdl.Color = blue;
    hE_lowEchoice_highROI_hdl.Marker = 'o';
    hE_lowEchoice_highROI_hdl.LineStyle = 'none';
    hE_lowEchoice_highROI_hdl.LineWidth = lWidth;
    hE_lowEchoice_highROI_hdl.Color = purple;
    hE_highEchoice_lowROI_hdl.Marker = 'o';
    hE_highEchoice_lowROI_hdl.LineStyle = 'none';
    hE_highEchoice_lowROI_hdl.LineWidth = lWidth;
    hE_highEchoice_lowROI_hdl.Color = blue;
    hE_highEchoice_highROI_hdl.Marker = 'o';
    hE_highEchoice_highROI_hdl.LineStyle = 'none';
    hE_highEchoice_highROI_hdl.LineWidth = lWidth;
    hE_highEchoice_highROI_hdl.Color = purple;
    % modify curve properties for fit
    hE_lowEchoice_lowROI_fit_hdl.LineStyle = '--';
    hE_lowEchoice_lowROI_fit_hdl.LineWidth = lWidth;
    hE_lowEchoice_lowROI_fit_hdl.Color = blue;
    hE_lowEchoice_highROI_fit_hdl.LineStyle = '--';
    hE_lowEchoice_highROI_fit_hdl.LineWidth = lWidth;
    hE_lowEchoice_highROI_fit_hdl.Color = purple;
    hE_highEchoice_lowROI_fit_hdl.LineStyle = '-';
    hE_highEchoice_lowROI_fit_hdl.LineWidth = lWidth;
    hE_highEchoice_lowROI_fit_hdl.Color = blue;
    hE_highEchoice_highROI_fit_hdl.LineStyle = '-';
    hE_highEchoice_highROI_fit_hdl.LineWidth = lWidth;
    hE_highEchoice_highROI_fit_hdl.Color = purple;
%     legend([hE_highEchoice_highROI_hdl,...
%         hE_highEchoice_lowROI_hdl,...
%         hE_lowEchoice_highROI_hdl,...
%         hE_lowEchoice_lowROI_hdl],...
%         {['hEch - high ',ROI_short_nm],...
%         ['hEch - low ',ROI_short_nm],...
%         ['lowEch - high ',ROI_short_nm],...
%         ['lowEch - low ',ROI_short_nm]});
%     legend('boxoff');
    xticks(hE_levels);
    xlabel('high effort level');
    switch Em_nm
        case 'successSpeed'
            % success speed = f(ROI)
            ylabel('success speed (s)');
        case 'n_errors'
            % number of errors made = f(ROI)
            ylabel('Number of errors');
        case 'latency'
            % latency = f(ROI)
            ylabel('latency (s)');
    end
    legend_size(pSize);
end % Em feature loop

%% Performance (%) for each task
fig;

% Performance (physical) = f(ROI)
subplot(1,2,1);
hold on;
perf_low_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.perf.perEch.(['low_',ROI_short_nm]),...
    sem_perf.physical.perf.perEch.(['low_',ROI_short_nm]));
perf_high_ROI_hdl = errorbar(Ech_levels,...
    m_perf.physical.perf.perEch.(['high_',ROI_short_nm]),...
    sem_perf.physical.perf.perEch.(['high_',ROI_short_nm]));
% add fit
perf_low_ROI_fit_hdl = plot(Ech_levels,...
    m_perf_fit.physical.perf.perEch.(['low_',ROI_short_nm]));
perf_high_ROI_fit_hdl = plot(Ech_levels,...
    m_perf_fit.physical.perf.perEch.(['high_',ROI_short_nm]));
% modify curve properties
perf_low_ROI_hdl.Marker = 'o';
perf_low_ROI_hdl.LineStyle = 'none';
perf_low_ROI_hdl.LineWidth = lWidth;
perf_low_ROI_hdl.Color = blue;
perf_high_ROI_hdl.Marker = 'o';
perf_high_ROI_hdl.LineStyle = 'none';
perf_high_ROI_hdl.LineWidth = lWidth;
perf_high_ROI_hdl.Color = purple;
% add fit
perf_low_ROI_fit_hdl.LineStyle = '--';
perf_low_ROI_fit_hdl.LineWidth = lWidth;
perf_low_ROI_fit_hdl.Color = blue;
perf_high_ROI_fit_hdl.LineStyle = '-';
perf_high_ROI_fit_hdl.LineWidth = lWidth;
perf_high_ROI_fit_hdl.Color = purple;
legend([perf_high_ROI_hdl, perf_low_ROI_hdl],...
    {['high ',ROI_short_nm],['low ',ROI_short_nm]});
legend('boxoff');
legend('Location','SouthWest');
xticks(Ech_levels);
xlabel('E chosen');
ylabel('Performance (%) - Physical task');
legend_size(pSize);

% Performance (mental) = f(ROI)
subplot(1,2,2);
hold on;
perf_low_ROI_hdl = errorbar(Ech_levels,...
    m_perf.mental.perf.perEch.(['low_',ROI_short_nm]),...
    sem_perf.mental.perf.perEch.(['low_',ROI_short_nm]));
perf_high_ROI_hdl = errorbar(Ech_levels,...
    m_perf.mental.perf.perEch.(['high_',ROI_short_nm]),...
    sem_perf.mental.perf.perEch.(['high_',ROI_short_nm]));
% add fit
perf_low_ROI_fit_hdl = plot(Ech_levels,...
    m_perf_fit.mental.perf.perEch.(['low_',ROI_short_nm]));
perf_high_ROI_fit_hdl = plot(Ech_levels,...
    m_perf_fit.mental.perf.perEch.(['high_',ROI_short_nm]));
% modify curve properties
perf_low_ROI_hdl.Marker = 'o';
perf_low_ROI_hdl.LineStyle = 'none';
perf_low_ROI_hdl.LineWidth = lWidth;
perf_low_ROI_hdl.Color = blue;
perf_high_ROI_hdl.Marker = 'o';
perf_high_ROI_hdl.LineStyle = 'none';
perf_high_ROI_hdl.LineWidth = lWidth;
perf_high_ROI_hdl.Color = purple;
% add fit
perf_low_ROI_fit_hdl.LineStyle = '--';
perf_low_ROI_fit_hdl.LineWidth = lWidth;
perf_low_ROI_fit_hdl.Color = blue;
perf_high_ROI_fit_hdl.LineStyle = '-';
perf_high_ROI_fit_hdl.LineWidth = lWidth;
perf_high_ROI_fit_hdl.Color = purple;
legend([perf_high_ROI_hdl, perf_low_ROI_hdl],...
    {['high ',ROI_short_nm],['low ',ROI_short_nm]});
legend('boxoff');
legend('Location','SouthWest');
xticks(Ech_levels);
xlabel('E chosen');
ylabel('Performance (%) - Mental task');
legend_size(pSize);