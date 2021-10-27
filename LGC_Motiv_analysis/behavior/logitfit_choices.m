%% check kR, kF and kE for pilots

resultsFolder = 'C:\Users\Loco\Documents\GitHub\LGC_motiv\LGC_Motiv_results\fMRI_pilots';

subid = 'pilot_s1';
% subid = 'pilot_s2';

switch subid
    case 'pilot_s1'
        EpData_r1 = load([resultsFolder,filesep,'AG_s1_session1_physical_task_behavioral_tmp.mat']);
        EmData_r1 = load([resultsFolder,filesep,'AG_s1_session2_mental_task_behavioral_tmp.mat']);
    case 'pilot_s2'
        EmData_r1 = load([resultsFolder,filesep,'LC_s2_session1_mental_task_behavioral_tmp.mat']);
end

n_bins = 6;
pSize = 40;

% extract kR and kEp
switch subid
    case 'pilot_s1'
        R_or_P_Ep = strcmp(summary.choiceOptions.R_or_P,'R');
        deltaR_Ep = EpData_r1.summary.choiceOptions.R.left(R_or_P_Ep == 1) - EpData_r1.summary.choiceOptions.R.right(R_or_P_Ep == 1);
        deltaEp = EpData_r1.summary.choiceOptions.E.left(R_or_P_Ep == 1) - EpData_r1.summary.choiceOptions.E.right(R_or_P_Ep == 1);
        choice_LR_Ep = (EpData_r1.summary.choice(R_or_P_Ep == 1) < 0);
        allTrials = 1:48;
        trialN = allTrials(R_or_P_Ep == 1);
        
        fatigue_Ep = deltaEp.*trialN;
        x = [deltaR_Ep', deltaEp', fatigue_Ep'];
        [betas_Ep] = glmfit(x, choice_LR_Ep', 'binomial','link','logit','constant','off');
        choice_LR_Ep_fit = glmval(betas_Ep, x,'logit','constant','off');
        
        % perform bins
        [choice_LR_Ep_deltaR_binned, deltaR_Ep_binned] = do_bin2(choice_LR_Ep, deltaR_Ep, n_bins, 0);
        [choice_LR_Ep_deltaE_binned, deltaE_Ep_binned] = do_bin2(choice_LR_Ep, deltaEp, n_bins, 0);
        [choice_LR_Ep_fit_deltaR_binned, deltaR_Ep_binned] = do_bin2(choice_LR_Ep_fit, deltaR_Ep, n_bins, 0);
        [choice_LR_Ep_fit_deltaE_binned, deltaE_Ep_binned] = do_bin2(choice_LR_Ep_fit, deltaEp, n_bins, 0);
        
        fig;
        subplot(1,2,1);
        scatter(deltaR_Ep_binned, choice_LR_Ep_deltaR_binned);
        hold on;
        plot(deltaR_Ep_binned, choice_LR_Ep_fit_deltaR_binned, 'k');
        ylim([-0.2 1.2])
        xlabel('R_l_e_f_t - R_r_i_g_h_t');
        ylabel('choice = left');
        legend_size(pSize);
        
        subplot(1,2,2);
        scatter(deltaE_Ep_binned, choice_LR_Ep_deltaE_binned);
        hold on;
        plot(deltaE_Ep_binned, choice_LR_Ep_fit_deltaE_binned);
        ylim([-0.2 1.2])
        xlabel('E_l_e_f_t - E_r_i_g_h_t');
        ylabel('choice = left');
        legend_size(pSize);
end

% extract kR and kEm
switch subid
    case {'pilot_s1','pilot_s2'}
        R_or_P_Em = strcmp(summary.choiceOptions.R_or_P,'R');
        deltaR_Em = EmData_r1.summary.choiceOptions.R.left(R_or_P_Em == 1) - EmData_r1.summary.choiceOptions.R.right(R_or_P_Em == 1);
        deltaEm = EmData_r1.summary.choiceOptions.E.left(R_or_P_Em == 1) - EmData_r1.summary.choiceOptions.E.right(R_or_P_Em == 1);
        choice_LR_Em = (EmData_r1.summary.choice(R_or_P_Em == 1) < 0);
        allTrials = 1:48;
        trialN = allTrials(R_or_P_Em == 1);
        
        fatigue_Em = deltaEm.*trialN;
        x = [deltaR_Em', deltaEm', fatigue_Em'];
        [betas_Em] = glmfit(x, choice_LR_Em', 'binomial','link','logit','constant','off');
        choice_LR_Em_fit = glmval(betas_Em, x,'logit','constant','off');
        
        % perform bins
        [choice_LR_Em_deltaR_binned, deltaR_Em_binned] = do_bin2(choice_LR_Em, deltaR_Em, n_bins, 0);
        [choice_LR_Em_deltaE_binned, deltaE_Em_binned] = do_bin2(choice_LR_Em, deltaEm, n_bins, 0);
        [choice_LR_Em_fit_deltaR_binned, deltaR_Em_binned] = do_bin2(choice_LR_Em_fit, deltaR_Em, n_bins, 0);
        [choice_LR_Em_fit_deltaE_binned, deltaE_Em_binned] = do_bin2(choice_LR_Em_fit, deltaEm, n_bins, 0);
        
        fig;
        subplot(1,2,1);
        scatter(deltaR_Em_binned, choice_LR_Em_deltaR_binned);
        hold on;
        plot(deltaR_Em_binned, choice_LR_Em_fit_deltaR_binned, 'k');
        ylim([-0.2 1.2])
        xlabel('R_l_e_f_t - R_r_i_g_h_t');
        ylabel('choice = left');
        legend_size(pSize);
        
        subplot(1,2,2);
        scatter(deltaE_Em_binned, choice_LR_Em_deltaE_binned);
        hold on;
        plot(deltaE_Em_binned, choice_LR_Em_fit_deltaE_binned);
        ylim([-0.2 1.2])
        xlabel('E_l_e_f_t - E_r_i_g_h_t');
        ylabel('choice = left');
        legend_size(pSize);
end
        