%% extract the level of activity of a given ROI (extracted trial by trial 
% with extract_ROI_betas_onsets_only_bis.m) in function of the varying
% reward, punishment and the level of effort chosen

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% working directory
computerRoot = LGCM_root_paths;
studyFolder = [computerRoot, study_nm, filesep];

%% load ROI
[condition_for_fMRI_extraction] = condition_for_fMRI(condition);
% [ROI_trial_b_trial, ROI_subList,...
%     ROI_nm, ROI_short_nm,...
%     ~, timePeriod_nm] = extract_ROI_betas_onsets_only_bis(computerRoot,...
%     study_nm, subject_id, condition_for_fMRI_extraction);

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
    task_nm1 = tasks{iT};
    % data for all trials
    [fMRI_ROI.(task_nm1).allTrials,...
        fMRI_ROI_fit.(task_nm1).f_dR.allTrials,...
        fMRI_ROI_fit.(task_nm1).f_dP.allTrials,...
        fMRI_ROI_fit.(task_nm1).f_dE.allTrials,...
        fMRI_ROI_fit.(task_nm1).f_dR_lowEch.allTrials,...
        fMRI_ROI_fit.(task_nm1).f_dR_highEch.allTrials,...
        fMRI_ROI_fit.(task_nm1).f_dP_lowEch.allTrials,...
        fMRI_ROI_fit.(task_nm1).f_dP_highEch.allTrials,...
        fMRI_ROI_fit.(task_nm1).f_dE_lowEch.allTrials,...
        fMRI_ROI_fit.(task_nm1).f_dE_highEch.allTrials,...
        choice_hE.(task_nm1).allTrials,...
        RP.(task_nm1).allTrials,...
        hR_level.(task_nm1).allTrials, Rchosen.(task_nm1).allTrials,...
        hP_level.(task_nm1).allTrials, Pchosen.(task_nm1).allTrials,...
        hE_level.(task_nm1).allTrials, Echosen.(task_nm1).allTrials,...
        run1_cstt.(task_nm1).allTrials, run2_cstt.(task_nm1).allTrials] = deal(NaN(nTrialsPerTask, NS));
    % data split by level of incentive proposed for high effort option
    [fMRI_ROI.(task_nm1).hR_level,...
        fMRI_ROI_fit.(task_nm1).hR_level,...
        choice_hE.(task_nm1).hR_level] = deal(NaN(n_hR_levels, NS));
    [fMRI_ROI.(task_nm1).hP_level,...
        fMRI_ROI_fit.(task_nm1).hP_level,...
        choice_hE.(task_nm1).hP_level] = deal(NaN(n_hP_levels, NS));
    
    % data split by level of effort proposed for high effort option
    [fMRI_ROI.(task_nm1).hE_level,...
        fMRI_ROI_fit.(task_nm1).hE_level,...
        choice_hE.(task_nm1).hE_level] = deal(NaN(n_hE_levels, NS));
    
    % data split by incentive*choice
    [fMRI_ROI.(task_nm1).lEch.hR_level,...
        fMRI_ROI_fit.(task_nm1).lEch.hR_level,...
        fMRI_ROI.(task_nm1).hEch.hR_level,...
        fMRI_ROI_fit.(task_nm1).hEch.hR_level] = deal(NaN(n_hR_levels, NS));
    [fMRI_ROI.(task_nm1).lEch.hP_level,...
        fMRI_ROI_fit.(task_nm1).lEch.hP_level,...
        fMRI_ROI.(task_nm1).hEch.hP_level,...
        fMRI_ROI_fit.(task_nm1).hEch.hP_level] = deal(NaN(n_hP_levels, NS));
    % data split by effort*choice
    [fMRI_ROI.(task_nm1).lEch.hE_level,...
        fMRI_ROI_fit.(task_nm1).lEch.hE_level,...
        fMRI_ROI.(task_nm1).hEch.hE_level,...
        fMRI_ROI_fit.(task_nm1).hEch.hE_level] = deal(NaN(n_hP_levels, NS));
    
    % fit for correlation with dR/dP/dE and dR/dP/dE*choice
    [betas.(task_nm1).ROI_f_dR.allTrials.b_r1,...
        betas.(task_nm1).ROI_f_dR.allTrials.b_r2,...
        betas.(task_nm1).ROI_f_dR.allTrials.bR,...
        betas.(task_nm1).ROI_f_dP.allTrials.b_r1,...
        betas.(task_nm1).ROI_f_dP.allTrials.b_r2,...
        betas.(task_nm1).ROI_f_dP.allTrials.bP,...
        betas.(task_nm1).ROI_f_dE.allTrials.b_r1,...
        betas.(task_nm1).ROI_f_dE.allTrials.b_r2,...
        betas.(task_nm1).ROI_f_dE.allTrials.bE,...
        betas.(task_nm1).ROI_f_dR.lowEch.b_r1,...
        betas.(task_nm1).ROI_f_dR.lowEch.b_r2,...
        betas.(task_nm1).ROI_f_dR.lowEch.bR,...
        betas.(task_nm1).ROI_f_dP.lowEch.b_r1,...
        betas.(task_nm1).ROI_f_dP.lowEch.b_r2,...
        betas.(task_nm1).ROI_f_dP.lowEch.bP,...
        betas.(task_nm1).ROI_f_dE.lowEch.b_r1,...
        betas.(task_nm1).ROI_f_dE.lowEch.b_r2,...
        betas.(task_nm1).ROI_f_dE.lowEch.bE,...
        betas.(task_nm1).ROI_f_dR.highEch.b_r1,...
        betas.(task_nm1).ROI_f_dR.highEch.b_r2,...
        betas.(task_nm1).ROI_f_dR.highEch.bR,...
        betas.(task_nm1).ROI_f_dP.highEch.b_r1,...
        betas.(task_nm1).ROI_f_dP.highEch.b_r2,...
        betas.(task_nm1).ROI_f_dP.highEch.bP,...
        betas.(task_nm1).ROI_f_dE.highEch.b_r1,...
        betas.(task_nm1).ROI_f_dE.highEch.b_r2,...
        betas.(task_nm1).ROI_f_dE.highEch.bE] = deal(NaN(1,NS));
end % task loop

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyFolder, 'CID',sub_nm, filesep, 'behavior',filesep];

    % extract runs
    [runsStruct, n_runs] = runs_definition(study_nm, sub_nm, condition);
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    for iRun = 1:n_runs
        task_nm2 = taskNames{iRun};
        kRun = okRuns(iRun);
        run_nm = num2str(kRun);
        switch kRun
            case {1,2}
                jRun = 1;
            case {3,4}
                jRun = 2;
        end
        run_nm_bis = ['run',num2str(jRun)];
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun - 1);
        [task_fullName] = task_fullName_extraction(task_nm2);
        
        %% extract relevant variables
        switch jRun
            case 1
                run1_cstt.(task_nm2).allTrials(runTrials_idx,iS) = 1;
                run2_cstt.(task_nm2).allTrials(runTrials_idx,iS) = 0;
            case 2
                run1_cstt.(task_nm2).allTrials(runTrials_idx,iS) = 0;
                run2_cstt.(task_nm2).allTrials(runTrials_idx,iS) = 1;
        end
        % R/P trials
        RP_var_tmp = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        R_trials = RP_var_tmp == 1;
        P_trials = RP_var_tmp == -1;
        RP.(task_nm2).allTrials(runTrials_idx,iS) = RP_var_tmp;
        % choice
        choice_hE.(task_nm2).allTrials(runTrials_idx,iS) = extract_choice_hE(subBehaviorFolder,sub_nm,run_nm,task_fullName);
        % extract R, P and E level for high effort option and for choice
        hR_level.(task_nm2).allTrials(runTrials_idx,iS) = extract_hR_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        hP_level.(task_nm2).allTrials(runTrials_idx,iS) = extract_hP_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        hE_level.(task_nm2).allTrials(runTrials_idx,iS) = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        %% extract fMRI ROI mediator
        fMRI_ROI.(task_nm2).allTrials(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm2).(run_nm_bis).(timePeriod_nm)(:, iS);
    end % run loop
    
    %% now perform the fit and average the data according to R/P/E and choice
    for iT = 1:nTasks
        task_nm3 = tasks{iT};
        
        %% perform the fits
        % reward
        R_okTrials = ~isnan(run1_cstt.(task_nm3).allTrials(:,iS).*hR_level.(task_nm3).allTrials(:,iS).*fMRI_ROI.(task_nm3).allTrials(:, iS));
        x_R_regs = [run1_cstt.(task_nm3).allTrials(R_okTrials,iS),...
            run2_cstt.(task_nm3).allTrials(R_okTrials,iS),...
            hR_level.(task_nm3).allTrials(R_okTrials,iS)];
        betas_R_tmp = glmfit(x_R_regs,...
            fMRI_ROI.(task_nm3).allTrials(R_okTrials, iS),...
            'normal','Constant','off');
        betas.(task_nm3).ROI_f_dR.allTrials.b_r1(iS) = betas_R_tmp(1);
        betas.(task_nm3).ROI_f_dR.allTrials.b_r2(iS) = betas_R_tmp(2);
        betas.(task_nm3).ROI_f_dR.allTrials.bR(iS) = betas_R_tmp(3);
        fMRI_ROI_fit.(task_nm3).f_dR.allTrials(R_okTrials, iS) = glmval(betas_R_tmp, x_R_regs, 'identity','Constant','off');
        
        % punishment
        P_okTrials = ~isnan(run1_cstt.(task_nm3).allTrials(:,iS).*hP_level.(task_nm3).allTrials(:,iS).*fMRI_ROI.(task_nm3).allTrials(:, iS));
        x_P_regs = [run1_cstt.(task_nm3).allTrials(P_okTrials,iS),...
            run2_cstt.(task_nm3).allTrials(P_okTrials,iS),...
            hP_level.(task_nm3).allTrials(P_okTrials,iS)];
        betas_P_tmp = glmfit(x_P_regs,...
            fMRI_ROI.(task_nm3).allTrials(P_okTrials, iS),...
            'normal','Constant','off');
        betas.(task_nm3).ROI_f_dP.allTrials.b_r1(iS) = betas_P_tmp(1);
        betas.(task_nm3).ROI_f_dP.allTrials.b_r2(iS) = betas_P_tmp(2);
        betas.(task_nm3).ROI_f_dP.allTrials.bP(iS) = betas_P_tmp(3);
        fMRI_ROI_fit.(task_nm3).f_dP.allTrials(P_okTrials, iS) = glmval(betas_P_tmp, x_P_regs, 'identity','Constant','off');
        
        % effort
        E_okTrials = ~isnan(run1_cstt.(task_nm3).allTrials(:,iS).*hE_level.(task_nm3).allTrials(:,iS).*fMRI_ROI.(task_nm3).allTrials(:, iS));
        x_E_regs = [run1_cstt.(task_nm3).allTrials(E_okTrials,iS),...
            run2_cstt.(task_nm3).allTrials(E_okTrials,iS),...
            hE_level.(task_nm3).allTrials(E_okTrials,iS)];
        betas_E_tmp = glmfit(x_E_regs,...
            fMRI_ROI.(task_nm3).allTrials(E_okTrials, iS),...
            'normal','Constant','off');
        betas.(task_nm3).ROI_f_dE.allTrials.b_r1(iS) = betas_E_tmp(1);
        betas.(task_nm3).ROI_f_dE.allTrials.b_r2(iS) = betas_E_tmp(2);
        betas.(task_nm3).ROI_f_dE.allTrials.bE(iS) = betas_E_tmp(3);
        fMRI_ROI_fit.(task_nm3).f_dE.allTrials(E_okTrials, iS) = glmval(betas_E_tmp, x_E_regs, 'identity','Constant','off');        
        
        % fits for dR/dP/dE*choice
        choice_hE_tmp = choice_hE.(task_nm3).allTrials(:, iS) == 1;
        choice_lE_tmp = choice_hE.(task_nm3).allTrials(:, iS) == 0;
        % reward
        % high E chosen
        R_hEch_okTrials = (R_okTrials.*choice_hE_tmp) == 1;
        x_R_hEch_regs = [run1_cstt.(task_nm3).allTrials(R_hEch_okTrials,iS),...
            run2_cstt.(task_nm3).allTrials(R_hEch_okTrials,iS),...
            hR_level.(task_nm3).allTrials(R_hEch_okTrials,iS)];
        betas_R_hEch_tmp = glmfit(x_R_hEch_regs,...
            fMRI_ROI.(task_nm3).allTrials(R_hEch_okTrials, iS),...
            'normal','Constant','off');
        betas.(task_nm3).ROI_f_dR.highEch.b_r1(iS) = betas_R_hEch_tmp(1);
        betas.(task_nm3).ROI_f_dR.highEch.b_r2(iS) = betas_R_hEch_tmp(2);
        betas.(task_nm3).ROI_f_dR.highEch.bR(iS) = betas_R_hEch_tmp(3);
        fMRI_ROI_fit.(task_nm3).f_dR_highEch.allTrials(R_hEch_okTrials, iS) = glmval(betas_R_hEch_tmp, x_R_hEch_regs, 'identity','Constant','off');
        % low E chosen
        R_lEch_okTrials = (R_okTrials.*choice_lE_tmp) == 1;
        x_R_lEch_regs = [run1_cstt.(task_nm3).allTrials(R_lEch_okTrials,iS),...
            run2_cstt.(task_nm3).allTrials(R_lEch_okTrials,iS),...
            hR_level.(task_nm3).allTrials(R_lEch_okTrials,iS)];
        betas_R_lEch_tmp = glmfit(x_R_lEch_regs,...
            fMRI_ROI.(task_nm3).allTrials(R_lEch_okTrials, iS),...
            'normal','Constant','off');
        betas.(task_nm3).ROI_f_dR.lowEch.b_r1(iS) = betas_R_lEch_tmp(1);
        betas.(task_nm3).ROI_f_dR.lowEch.b_r2(iS) = betas_R_lEch_tmp(2);
        betas.(task_nm3).ROI_f_dR.lowEch.bR(iS) = betas_R_lEch_tmp(3);
        fMRI_ROI_fit.(task_nm3).f_dR_lowEch.allTrials(R_lEch_okTrials, iS) = glmval(betas_R_lEch_tmp, x_R_lEch_regs, 'identity','Constant','off');
        
        % punishment
        % high E chosen
        P_hEch_okTrials = (P_okTrials.*choice_hE_tmp) == 1;
        x_P_hEch_regs = [run1_cstt.(task_nm3).allTrials(P_hEch_okTrials,iS),...
            run2_cstt.(task_nm3).allTrials(P_hEch_okTrials,iS),...
            hP_level.(task_nm3).allTrials(P_hEch_okTrials,iS)];
        betas_P_hEch_tmp = glmfit(x_P_hEch_regs,...
            fMRI_ROI.(task_nm3).allTrials(P_hEch_okTrials, iS),...
            'normal','Constant','off');
        betas.(task_nm3).ROI_f_dP.highEch.b_r1(iS) = betas_P_hEch_tmp(1);
        betas.(task_nm3).ROI_f_dP.highEch.b_r2(iS) = betas_P_hEch_tmp(2);
        betas.(task_nm3).ROI_f_dP.highEch.bP(iS) = betas_P_hEch_tmp(3);
        fMRI_ROI_fit.(task_nm3).f_dP_highEch.allTrials(P_hEch_okTrials, iS) = glmval(betas_P_hEch_tmp, x_P_hEch_regs, 'identity','Constant','off');
        % low E chosen
        P_lEch_okTrials = (P_okTrials.*choice_lE_tmp) == 1;
        x_P_lEch_regs = [run1_cstt.(task_nm3).allTrials(P_lEch_okTrials,iS),...
            run2_cstt.(task_nm3).allTrials(P_lEch_okTrials,iS),...
            hP_level.(task_nm3).allTrials(P_lEch_okTrials,iS)];
        betas_P_lEch_tmp = glmfit(x_P_lEch_regs,...
            fMRI_ROI.(task_nm3).allTrials(P_lEch_okTrials, iS),...
            'normal','Constant','off');
        betas.(task_nm3).ROI_f_dP.lowEch.b_r1(iS) = betas_P_lEch_tmp(1);
        betas.(task_nm3).ROI_f_dP.lowEch.b_r2(iS) = betas_P_lEch_tmp(2);
        betas.(task_nm3).ROI_f_dP.lowEch.bP(iS) = betas_P_lEch_tmp(3);
        fMRI_ROI_fit.(task_nm3).f_dP_lowEch.allTrials(P_lEch_okTrials, iS) = glmval(betas_P_lEch_tmp, x_P_lEch_regs, 'identity','Constant','off');
        
        
        % effort
        % high E chosen
        E_hEch_okTrials = (E_okTrials.*choice_hE_tmp) == 1;
        x_E_hEch_regs = [run1_cstt.(task_nm3).allTrials(E_hEch_okTrials,iS),...
            run2_cstt.(task_nm3).allTrials(E_hEch_okTrials,iS),...
            hE_level.(task_nm3).allTrials(E_hEch_okTrials,iS)];
        betas_E_hEch_tmp = glmfit(x_E_hEch_regs,...
            fMRI_ROI.(task_nm3).allTrials(E_hEch_okTrials, iS),...
            'normal','Constant','off');
        betas.(task_nm3).ROI_f_dE.highEch.b_r1(iS) = betas_E_hEch_tmp(1);
        betas.(task_nm3).ROI_f_dE.highEch.b_r2(iS) = betas_E_hEch_tmp(2);
        betas.(task_nm3).ROI_f_dE.highEch.bE(iS) = betas_E_hEch_tmp(3);
        fMRI_ROI_fit.(task_nm3).f_dE_highEch.allTrials(E_hEch_okTrials, iS) = glmval(betas_E_hEch_tmp, x_E_hEch_regs, 'identity','Constant','off');
        % low E chosen
        E_lEch_okTrials = (E_okTrials.*choice_lE_tmp) == 1;
        x_E_lEch_regs = [run1_cstt.(task_nm3).allTrials(E_lEch_okTrials,iS),...
            run2_cstt.(task_nm3).allTrials(E_lEch_okTrials,iS),...
            hE_level.(task_nm3).allTrials(E_lEch_okTrials,iS)];
        betas_E_lEch_tmp = glmfit(x_E_lEch_regs,...
            fMRI_ROI.(task_nm3).allTrials(E_lEch_okTrials, iS),...
            'normal','Constant','off');
        betas.(task_nm3).ROI_f_dE.lowEch.b_r1(iS) = betas_E_lEch_tmp(1);
        betas.(task_nm3).ROI_f_dE.lowEch.b_r2(iS) = betas_E_lEch_tmp(2);
        betas.(task_nm3).ROI_f_dE.lowEch.bE(iS) = betas_E_lEch_tmp(3);
        fMRI_ROI_fit.(task_nm3).f_dE_lowEch.allTrials(E_lEch_okTrials, iS) = glmval(betas_E_lEch_tmp, x_E_lEch_regs, 'identity','Constant','off');
        
        %% extract data per high effort option features
        % extract data per reward level
        for iR = 1:n_hR_levels
            R_idx = hR_level.(task_nm3).allTrials(:,iS) == iR;
            fMRI_ROI.(task_nm3).hR_level(iR,iS) = mean(fMRI_ROI.(task_nm3).allTrials(R_idx, iS),1,'omitnan');
            choice_hE.(task_nm3).hR_level(iR,iS) = mean(choice_hE.(task_nm3).allTrials(R_idx, iS),1,'omitnan');
            fMRI_ROI_fit.(task_nm3).hR_level(iR,iS) = mean(fMRI_ROI_fit.(task_nm3).f_dR.allTrials(R_idx, iS),1,'omitnan');
        end % reward

        % extract data per punishment level
        for iP = 1:n_hP_levels
            P_idx = hP_level.(task_nm3).allTrials(:,iS) == iP;
            fMRI_ROI.(task_nm3).hP_level(iP,iS) = mean(fMRI_ROI.(task_nm3).allTrials(P_idx, iS),1,'omitnan');
            choice_hE.(task_nm3).hP_level(iP,iS) = mean(choice_hE.(task_nm3).allTrials(P_idx, iS),1,'omitnan');
            fMRI_ROI_fit.(task_nm3).hP_level(iP,iS) = mean(fMRI_ROI_fit.(task_nm3).f_dP.allTrials(P_idx, iS),1,'omitnan');
        end % punishment

        % extract data per effort level
        for iE = 1:n_hE_levels
            E_idx = hE_level.(task_nm3).allTrials(:,iS) == iE;
            fMRI_ROI.(task_nm3).hE_level(iE,iS) = mean(fMRI_ROI.(task_nm3).allTrials(E_idx, iS),1,'omitnan');
            choice_hE.(task_nm3).hE_level(iE,iS) = mean(choice_hE.(task_nm3).allTrials(E_idx, iS),1,'omitnan');
            fMRI_ROI_fit.(task_nm3).hE_level(iE,iS) = mean(fMRI_ROI_fit.(task_nm3).f_dE.allTrials(E_idx, iS),1,'omitnan');
        end % effort

        %% extract data depending on choice (high/low) and high effort option features
        choice_highE_idx_tmp    = choice_hE.(task_nm3).allTrials(:, iS) == 1;

        % reward
        for iR = 1:n_hR_levels
            R_hE_idx = (hR_level.(task_nm3).allTrials(:,iS) == iR).*(choice_highE_idx_tmp == 1) == 1;
            R_lE_idx = (hR_level.(task_nm3).allTrials(:,iS) == iR).*(choice_highE_idx_tmp == 0) == 1;
            fMRI_ROI.(task_nm3).hEch.hR_level(iR,iS) = mean(fMRI_ROI.(task_nm3).allTrials(R_hE_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm3).lEch.hR_level(iR,iS) = mean(fMRI_ROI.(task_nm3).allTrials(R_lE_idx, iS),1,'omitnan');
            % fit
            fMRI_ROI_fit.(task_nm3).hEch.hR_level(iR,iS) = mean(fMRI_ROI_fit.(task_nm3).f_dR_highEch.allTrials(R_hE_idx, iS),1,'omitnan');
            fMRI_ROI_fit.(task_nm3).lEch.hR_level(iR,iS) = mean(fMRI_ROI_fit.(task_nm3).f_dR_lowEch.allTrials(R_lE_idx, iS),1,'omitnan');
        end

        % punishment
        for iP = 1:n_hP_levels
            P_hE_idx = (hP_level.(task_nm3).allTrials(:,iS) == iP).*(choice_highE_idx_tmp == 1) == 1;
            P_lE_idx = (hP_level.(task_nm3).allTrials(:,iS) == iP).*(choice_highE_idx_tmp == 0) == 1;
            fMRI_ROI.(task_nm3).hEch.hP_level(iP,iS) = mean(fMRI_ROI.(task_nm3).allTrials(P_hE_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm3).lEch.hP_level(iP,iS) = mean(fMRI_ROI.(task_nm3).allTrials(P_lE_idx, iS),1,'omitnan');
            % fit
            fMRI_ROI_fit.(task_nm3).hEch.hP_level(iP,iS) = mean(fMRI_ROI_fit.(task_nm3).f_dP_highEch.allTrials(P_hE_idx, iS),1,'omitnan');
            fMRI_ROI_fit.(task_nm3).lEch.hP_level(iP,iS) = mean(fMRI_ROI_fit.(task_nm3).f_dP_lowEch.allTrials(P_lE_idx, iS),1,'omitnan');
        end

        % effort
        for iE = 1:n_hE_levels
            E_hE_idx = (hE_level.(task_nm3).allTrials(:,iS) == iE).*(choice_highE_idx_tmp == 1) == 1;
            E_lE_idx = (hE_level.(task_nm3).allTrials(:,iS) == iE).*(choice_highE_idx_tmp == 0) == 1;
            fMRI_ROI.(task_nm3).hEch.hE_level(iE,iS) = mean(fMRI_ROI.(task_nm3).allTrials(E_hE_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm3).lEch.hE_level(iE,iS) = mean(fMRI_ROI.(task_nm3).allTrials(E_lE_idx, iS),1,'omitnan');
            % fit
            fMRI_ROI_fit.(task_nm3).hEch.hE_level(iE,iS) = mean(fMRI_ROI_fit.(task_nm3).f_dE_highEch.allTrials(E_hE_idx, iS),1,'omitnan');
            fMRI_ROI_fit.(task_nm3).lEch.hE_level(iE,iS) = mean(fMRI_ROI_fit.(task_nm3).f_dE_lowEch.allTrials(E_lE_idx, iS),1,'omitnan');
        end
    end % task loop
end % subject loop

%% figure parameter
[pSize, lWidth, col, mSize] = general_fig_prm;
figHighE = fig;
figChoice_x_highE = fig;
xRange = [0.8 3.2];
y_ROI_range = [-4.2 6.5];

%% loop through tasks for average + figure
for iT = 1:nTasks
    task_nm4 = tasks{iT};

    %% average across subjects
    % data split by level of incentive proposed for high effort option
    [mean_fMRI_ROI.(task_nm4).hR_level,...
        sem_fMRI_ROI.(task_nm4).hR_level] = mean_sem_sd(fMRI_ROI.(task_nm4).hR_level, 2);
    [mean_fMRI_ROI_fit.(task_nm4).hR_level,...
        sem_fMRI_ROI_fit.(task_nm4).hR_level] = mean_sem_sd(fMRI_ROI_fit.(task_nm4).hR_level, 2);
    [mean_choice_hE.(task_nm4).hR_level,...
        sem_choice_hE.(task_nm4).hR_level] = mean_sem_sd(choice_hE.(task_nm4).hR_level, 2);
    [mean_fMRI_ROI.(task_nm4).hP_level,...
        sem_fMRI_ROI.(task_nm4).hP_level] = mean_sem_sd(fMRI_ROI.(task_nm4).hP_level, 2);
    [mean_fMRI_ROI_fit.(task_nm4).hP_level,...
        sem_fMRI_ROI_fit.(task_nm4).hP_level] = mean_sem_sd(fMRI_ROI_fit.(task_nm4).hP_level, 2);
    [mean_choice_hE.(task_nm4).hP_level,...
        sem_choice_hE.(task_nm4).hP_level] = mean_sem_sd(choice_hE.(task_nm4).hP_level, 2);
    
    % data split by level of effort proposed for high effort option
    [mean_fMRI_ROI.(task_nm4).hE_level,...
        sem_fMRI_ROI.(task_nm4).hE_level] = mean_sem_sd(fMRI_ROI.(task_nm4).hE_level, 2);
    [mean_fMRI_ROI_fit.(task_nm4).hE_level,...
        sem_fMRI_ROI_fit.(task_nm4).hE_level] = mean_sem_sd(fMRI_ROI_fit.(task_nm4).hE_level, 2);
    [mean_choice_hE.(task_nm4).hE_level,...
        sem_choice_hE.(task_nm4).hE_level] = mean_sem_sd(choice_hE.(task_nm4).hE_level, 2);
    
    % data split by incentive*choice
    [mean_fMRI_ROI.(task_nm4).lEch.hR_level,...
        sem_fMRI_ROI.(task_nm4).lEch.hR_level] = mean_sem_sd(fMRI_ROI.(task_nm4).lEch.hR_level, 2);
    [mean_fMRI_ROI.(task_nm4).hEch.hR_level,...
        sem_fMRI_ROI.(task_nm4).hEch.hR_level] = mean_sem_sd(fMRI_ROI.(task_nm4).hEch.hR_level, 2);
    [mean_fMRI_ROI.(task_nm4).lEch.hP_level,...
        sem_fMRI_ROI.(task_nm4).lEch.hP_level] = mean_sem_sd(fMRI_ROI.(task_nm4).lEch.hP_level, 2);
    [mean_fMRI_ROI.(task_nm4).hEch.hP_level,...
        sem_fMRI_ROI.(task_nm4).hEch.hP_level] = mean_sem_sd(fMRI_ROI.(task_nm4).hEch.hP_level, 2);
    % same for fit
    [mean_fMRI_ROI_fit.(task_nm4).lEch.hR_level,...
        sem_fMRI_ROI_fit.(task_nm4).lEch.hR_level] = mean_sem_sd(fMRI_ROI_fit.(task_nm4).lEch.hR_level, 2);
    [mean_fMRI_ROI_fit.(task_nm4).hEch.hR_level,...
        sem_fMRI_ROI_fit.(task_nm4).hEch.hR_level] = mean_sem_sd(fMRI_ROI_fit.(task_nm4).hEch.hR_level, 2);
    [mean_fMRI_ROI_fit.(task_nm4).lEch.hP_level,...
        sem_fMRI_ROI_fit.(task_nm4).lEch.hP_level] = mean_sem_sd(fMRI_ROI_fit.(task_nm4).lEch.hP_level, 2);
    [mean_fMRI_ROI_fit.(task_nm4).hEch.hP_level,...
        sem_fMRI_ROI_fit.(task_nm4).hEch.hP_level] = mean_sem_sd(fMRI_ROI_fit.(task_nm4).hEch.hP_level, 2);

    % data split by effort*choice
    [mean_fMRI_ROI.(task_nm4).lEch.hE_level,...
        sem_fMRI_ROI.(task_nm4).lEch.hE_level] = mean_sem_sd(fMRI_ROI.(task_nm4).lEch.hE_level, 2);
    [mean_fMRI_ROI.(task_nm4).hEch.hE_level,...
        sem_fMRI_ROI.(task_nm4).hEch.hE_level] = mean_sem_sd(fMRI_ROI.(task_nm4).hEch.hE_level, 2);
    % same for fit
    [mean_fMRI_ROI_fit.(task_nm4).lEch.hE_level,...
        sem_fMRI_ROI_fit.(task_nm4).lEch.hE_level] = mean_sem_sd(fMRI_ROI_fit.(task_nm4).lEch.hE_level, 2);
    [mean_fMRI_ROI_fit.(task_nm4).hEch.hE_level,...
        sem_fMRI_ROI_fit.(task_nm4).hEch.hE_level] = mean_sem_sd(fMRI_ROI_fit.(task_nm4).hEch.hE_level, 2);
    
    %% average betas
    % across all trials
    [betas.(task_nm4).ROI_f_dR.allTrials.m_bR,...
        betas.(task_nm4).ROI_f_dR.allTrials.sem_bR] = mean_sem_sd(betas.(task_nm4).ROI_f_dR.allTrials.bR, 2);
    [betas.(task_nm4).ROI_f_dP.allTrials.m_bP,...
        betas.(task_nm4).ROI_f_dP.allTrials.sem_bP] = mean_sem_sd(betas.(task_nm4).ROI_f_dP.allTrials.bP, 2);
    [betas.(task_nm4).ROI_f_dE.allTrials.m_bE,...
        betas.(task_nm4).ROI_f_dE.allTrials.sem_bE] = mean_sem_sd(betas.(task_nm4).ROI_f_dE.allTrials.bE, 2);
    
    % split per type of effort chosen
    % reward
    [betas.(task_nm4).ROI_f_dR.highEch.m_bR,...
        betas.(task_nm4).ROI_f_dR.highEch.sem_bR] = mean_sem_sd(betas.(task_nm4).ROI_f_dR.highEch.bR, 2);
    [betas.(task_nm4).ROI_f_dR.lowEch.m_bR,...
        betas.(task_nm4).ROI_f_dR.lowEch.sem_bR] = mean_sem_sd(betas.(task_nm4).ROI_f_dR.lowEch.bR, 2);
    % punishment
    [betas.(task_nm4).ROI_f_dP.highEch.m_bP,...
        betas.(task_nm4).ROI_f_dP.highEch.sem_bP] = mean_sem_sd(betas.(task_nm4).ROI_f_dP.highEch.bP, 2);
    [betas.(task_nm4).ROI_f_dP.lowEch.m_bP,...
        betas.(task_nm4).ROI_f_dP.lowEch.sem_bP] = mean_sem_sd(betas.(task_nm4).ROI_f_dP.lowEch.bP, 2);
    % effort
    [betas.(task_nm4).ROI_f_dE.highEch.m_bE,...
        betas.(task_nm4).ROI_f_dE.highEch.sem_bE] = mean_sem_sd(betas.(task_nm4).ROI_f_dE.highEch.bE, 2);
    [betas.(task_nm4).ROI_f_dE.lowEch.m_bE,...
        betas.(task_nm4).ROI_f_dE.lowEch.sem_bE] = mean_sem_sd(betas.(task_nm4).ROI_f_dE.lowEch.bE, 2);
    
    %% t.test on betas + comparing slopes
    
    % are betas of the slopes significant?
    % pool across all trials
    [~,pval.(task_nm4).ROI_f_dR.allTrials.bR] = ttest(betas.(task_nm4).ROI_f_dR.allTrials.bR);
    [~,pval.(task_nm4).ROI_f_dP.allTrials.bP] = ttest(betas.(task_nm4).ROI_f_dP.allTrials.bP);
    [~,pval.(task_nm4).ROI_f_dE.allTrials.bE] = ttest(betas.(task_nm4).ROI_f_dE.allTrials.bE);
    
    % only high effort chosen slopes
    [~,pval.(task_nm4).ROI_f_dR.highEch.bR] = ttest(betas.(task_nm4).ROI_f_dR.highEch.bR);
    [~,pval.(task_nm4).ROI_f_dP.highEch.bP] = ttest(betas.(task_nm4).ROI_f_dP.highEch.bP);
    [~,pval.(task_nm4).ROI_f_dE.highEch.bE] = ttest(betas.(task_nm4).ROI_f_dE.highEch.bE);
    
    % only low effort chosen slopes
    [~,pval.(task_nm4).ROI_f_dR.lowEch.bR] = ttest(betas.(task_nm4).ROI_f_dR.lowEch.bR);
    [~,pval.(task_nm4).ROI_f_dP.lowEch.bP] = ttest(betas.(task_nm4).ROI_f_dP.lowEch.bP);
    [~,pval.(task_nm4).ROI_f_dE.lowEch.bE] = ttest(betas.(task_nm4).ROI_f_dE.lowEch.bE);
    
    % comparing slopes (with paired t.test)
    [~,pval.(task_nm4).ROI_f_dR.high_vs_low_Ech.bR] = ttest(betas.(task_nm4).ROI_f_dR.highEch.bR, betas.(task_nm4).ROI_f_dR.lowEch.bR);
    [~,pval.(task_nm4).ROI_f_dP.high_vs_low_Ech.bP] = ttest(betas.(task_nm4).ROI_f_dP.highEch.bP, betas.(task_nm4).ROI_f_dP.lowEch.bP);
    [~,pval.(task_nm4).ROI_f_dE.high_vs_low_Ech.bE] = ttest(betas.(task_nm4).ROI_f_dE.highEch.bE, betas.(task_nm4).ROI_f_dE.lowEch.bE);
    
    %% figures

    %% figure focused on varying option
    figure(figHighE);

    % reward
    subplot(nTasks,3,1 + 3*(iT-1)); hold on;
    R_er_hdl = errorbar(1:n_hR_levels,...
        mean_fMRI_ROI.(task_nm4).hR_level,...
        sem_fMRI_ROI.(task_nm4).hR_level);
    errorbar_hdl_upgrade(R_er_hdl);
    R_fit_hdl = plot(1:n_hR_levels,...
        mean_fMRI_ROI_fit.(task_nm4).hR_level);
    fit_hdl_upgrade(R_fit_hdl);
    ylim(y_ROI_range);
    xlim(xRange);
    xlabel('Reward');
    ylabel(ROI_short_nm);
    legend_size(pSize);

    % punishment
    subplot(nTasks,3,2 + 3*(iT-1)); hold on;
    P_er_hdl = errorbar(1:n_hP_levels,...
        mean_fMRI_ROI.(task_nm4).hP_level,...
        sem_fMRI_ROI.(task_nm4).hP_level);
    errorbar_hdl_upgrade(P_er_hdl);
    P_fit_hdl = plot(1:n_hP_levels,...
        mean_fMRI_ROI_fit.(task_nm4).hP_level);
    fit_hdl_upgrade(P_fit_hdl);
    ylim(y_ROI_range);
    xlim(xRange);
    xlabel('Punishment');
    ylabel(ROI_short_nm);
    legend_size(pSize);

    % effort
    subplot(nTasks,3,3 + 3*(iT-1)); hold on;
    E_er_hdl = errorbar(1:n_hE_levels,...
        mean_fMRI_ROI.(task_nm4).hE_level,...
        sem_fMRI_ROI.(task_nm4).hE_level);
    errorbar_hdl_upgrade(E_er_hdl);
    E_fit_hdl = plot(1:n_hE_levels,...
        mean_fMRI_ROI_fit.(task_nm4).hE_level);
    fit_hdl_upgrade(E_fit_hdl);
    ylim(y_ROI_range);
    xlim(xRange);
    xlabel('Effort');
    ylabel(ROI_short_nm);
    legend_size(pSize);


    %% figure for choice*high effort option
    figure(figChoice_x_highE);

    % reward
    subplot(nTasks,3,1 + 3*(iT-1)); hold on;
    R_hE_er_hdl = errorbar(1:n_hR_levels,...
        mean_fMRI_ROI.(task_nm4).hEch.hR_level,...
        sem_fMRI_ROI.(task_nm4).hEch.hR_level);
    R_lE_er_hdl = errorbar(1:n_hR_levels,...
        mean_fMRI_ROI.(task_nm4).lEch.hR_level,...
        sem_fMRI_ROI.(task_nm4).lEch.hR_level);
    R_hE_fit_hdl = plot(1:n_hR_levels,...
        mean_fMRI_ROI_fit.(task_nm4).hEch.hR_level);
    R_lE_fit_hdl = plot(1:n_hR_levels,...
        mean_fMRI_ROI_fit.(task_nm4).lEch.hR_level);
    errorbar_hdl_upgrade(R_hE_er_hdl, col.orange);
    errorbar_hdl_upgrade(R_lE_er_hdl, col.blue);
    fit_hdl_upgrade(R_hE_fit_hdl, col.orange_light);
    fit_hdl_upgrade(R_lE_fit_hdl, col.blue_light);
    ylim(y_ROI_range);
    xlim(xRange);
    xlabel('Reward');
    ylabel(ROI_short_nm);
    legend_size(pSize);

    % punishment
    subplot(nTasks,3,2 + 3*(iT-1)); hold on;
    P_hE_er_hdl = errorbar(1:n_hP_levels,...
        mean_fMRI_ROI.(task_nm4).hEch.hP_level,...
        sem_fMRI_ROI.(task_nm4).hEch.hP_level);
    P_lE_er_hdl = errorbar(1:n_hP_levels,...
        mean_fMRI_ROI.(task_nm4).lEch.hP_level,...
        sem_fMRI_ROI.(task_nm4).lEch.hP_level);
    P_hE_fit_hdl = plot(1:n_hR_levels,...
        mean_fMRI_ROI_fit.(task_nm4).hEch.hP_level);
    P_lE_fit_hdl = plot(1:n_hR_levels,...
        mean_fMRI_ROI_fit.(task_nm4).lEch.hP_level);
    errorbar_hdl_upgrade(P_hE_er_hdl, col.orange);
    errorbar_hdl_upgrade(P_lE_er_hdl, col.blue);
    fit_hdl_upgrade(P_hE_fit_hdl, col.orange_light);
    fit_hdl_upgrade(P_lE_fit_hdl, col.blue_light);
    ylim(y_ROI_range);
    xlim(xRange);
    xlabel('Punishment');
    ylabel(ROI_short_nm);
    legend_size(pSize);

    % effort
    subplot(nTasks,3,3 + 3*(iT-1)); hold on;
    E_hE_er_hdl = errorbar(1:n_hE_levels,...
        mean_fMRI_ROI.(task_nm4).hEch.hE_level,...
        sem_fMRI_ROI.(task_nm4).hEch.hE_level);
    E_lE_er_hdl = errorbar(1:n_hE_levels,...
        mean_fMRI_ROI.(task_nm4).lEch.hE_level,...
        sem_fMRI_ROI.(task_nm4).lEch.hE_level);
    E_hE_fit_hdl = plot(1:n_hE_levels,...
        mean_fMRI_ROI_fit.(task_nm4).hEch.hE_level);
    E_lE_fit_hdl = plot(1:n_hE_levels,...
        mean_fMRI_ROI_fit.(task_nm4).lEch.hE_level);
    errorbar_hdl_upgrade(E_hE_er_hdl, col.orange);
    errorbar_hdl_upgrade(E_lE_er_hdl, col.blue);
    fit_hdl_upgrade(E_hE_fit_hdl, col.orange_light);
    fit_hdl_upgrade(E_lE_fit_hdl, col.blue_light);
    ylim(y_ROI_range);
    xlim(xRange);
    xlabel('Effort');
    ylabel(ROI_short_nm);
    legend_size(pSize);
end % task loop