% choice_f_ROI__metabolites_mSplit will extract the choices depending on
% ROI level with choice_f_ROI.m. Then it will perform a median split on the
% data based on the level of metabolites (metabolite to be selected
% manually with medSplit_metabolites.m). You will then be able to see
% whether the proportion of high effort choices varies depending on the
% level of metabolites.
%
% See also choice_f_ROI.m and medSplit_metabolites.m

%% subject selection
study_nm = 'study1';
nBins = 6;
condition = subject_condition;
subject_id = LGCM_subject_selection(study_nm, condition);

%% general parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);

%% load data
figDisp = 0;
[b_choice_f_fMRI, pval, fMRI_bins,...
    choice_hE_bins, choice_hE_fit_bins] = choice_f_ROI(nBins, study_nm, subject_id, condition, figDisp);

%% extract level of metabolites
[low_met_subs, high_met_subs, metabolite_nm, MRS_ROI_nm,...
    metabolite_allSubs] = medSplit_metabolites(study_nm, subject_id);
low_mb_nm = ['low_',MRS_ROI_nm,'_',metabolite_nm];
high_mb_nm = ['high_',MRS_ROI_nm,'_',metabolite_nm];
low_mb_nm_bis = ['low ',MRS_ROI_nm,' ',metabolite_nm];
high_mb_nm_bis = ['high ',MRS_ROI_nm,' ',metabolite_nm];
%% loop on tasks
for iT = 1:nTasks
    task_nm_tmp = tasks{iT};
    
    %% average and split data depending on metabolite levels
    for iM = 1:2
        switch iM
            case 1
                mb_nm = low_mb_nm;
                mb_idx = low_met_subs;
            case 2
                mb_nm = high_mb_nm;
                mb_idx = high_met_subs;
        end
        
        % all trials
        [m_fMRI_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm),...
            sem_fMRI_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm)] = mean_sem_sd(fMRI_bins.(task_nm_tmp).allTrials(:,mb_idx),2);
        [m_choice_hE_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm),...
            sem_choice_hE_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm)] = mean_sem_sd(choice_hE_bins.(task_nm_tmp).allTrials(:,mb_idx),2);
        [m_choice_hE_fit_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm),...
            sem_choice_hE_fit_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm)] = mean_sem_sd(choice_hE_fit_bins.(task_nm_tmp).allTrials(:,mb_idx),2);
        b_choice_f_fMRI_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).allTrials.(mb_nm) = b_choice_f_fMRI.(task_nm_tmp).allTrials(2,mb_idx);
        [~,pval_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).choice_f_fMRI.allTrials.(mb_nm)] = ttest(b_choice_f_fMRI_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).allTrials.(mb_nm));
        for iE = 1:n_E_levels
            b_choice_f_fMRI_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).perElevel.(mb_nm).(['E',num2str(iE)]) = b_choice_f_fMRI.(task_nm_tmp).perElevel(2,mb_idx,iE);
            [~,pval_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).choice_f_fMRI.perE.(mb_nm).(['E',num2str(iE)])] = ttest(b_choice_f_fMRI_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).perElevel.(mb_nm).(['E',num2str(iE)]));
        end
        % split per E level
        [m_fMRI_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm),...
            sem_fMRI_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm)] = mean_sem_sd(fMRI_bins.(task_nm_tmp).perElevel(:,:,mb_idx),3);
        [m_choice_hE_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm),...
            sem_choice_hE_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm)] = mean_sem_sd(choice_hE_bins.(task_nm_tmp).perElevel(:,:,mb_idx),3);
        [m_choice_hE_fit_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm),...
            sem_choice_hE_fit_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm)] = mean_sem_sd(choice_hE_fit_bins.(task_nm_tmp).perElevel(:,:,mb_idx),3);
    end % loop on median split
    
    % compare betas of low vs high metabolite group
    [~,pval_mSplit.(task_nm_tmp).(['low_vs_high_',MRS_ROI_nm,'_',metabolite_nm]).allTrials] = ttest2(b_choice_f_fMRI_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).allTrials.(low_mb_nm),...
        b_choice_f_fMRI_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).allTrials.(high_mb_nm));
    for iE = 1:n_E_levels
        [~,pval_mSplit.(task_nm_tmp).(['low_vs_high_',MRS_ROI_nm,'_',metabolite_nm]).perElevel.(['E',num2str(iE)])] = ttest2(b_choice_f_fMRI_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).perElevel.(low_mb_nm).(['E',num2str(iE)]),...
        b_choice_f_fMRI_mSplit.(task_nm_tmp).([metabolite_nm,'_split']).perElevel.(high_mb_nm).(['E',num2str(iE)]));
    end % effort level
    
    %% figure
    pSize = 30;
    lWidth = 3;
    low_mb_col = [103 169 207]./255;
    high_mb_col = [239 138 98]./255;
    E1_col_low = [222 235 247]./255;
    E2_col_low = [158 202 225]./255;
    E3_col_low = [49 130 189]./255;
    E1_col_high = [254 224 210]./255;
    E2_col_high = [252 146 114]./255;
    E3_col_high = [222 45 38]./255;
    %% choice = f(ROI)
    fig;
    for iM = 1:2
        switch iM
            case 1
                mb_nm = low_mb_nm;
                curve_col = low_mb_col;
            case 2
                mb_nm = high_mb_nm;
                curve_col = high_mb_col;
        end
        choice_f_fMRI_hdl_mSplit.(mb_nm) = errorbar(m_fMRI_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm),...
            m_choice_hE_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm).*100,...
            sem_choice_hE_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm).*100);
        choice_f_fMRI_hdl_mSplit.(mb_nm).LineWidth = lWidth;
        choice_f_fMRI_hdl_mSplit.(mb_nm).Marker = 'o';
        choice_f_fMRI_hdl_mSplit.(mb_nm).Color = curve_col;
        choice_f_fMRI_hdl_mSplit.(mb_nm).LineStyle = 'none';
        choice_fit_f_fMRI_hdl_mSplit = plot(m_fMRI_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm),...
            m_choice_hE_fit_bins_mSplit.(task_nm_tmp).allTrials.(mb_nm).*100);
        choice_fit_f_fMRI_hdl_mSplit.LineWidth = lWidth;
        choice_fit_f_fMRI_hdl_mSplit.Color = curve_col;
        choice_fit_f_fMRI_hdl_mSplit.LineStyle = '--';
    end
    line([0 0],[0 100],'Color','k','LineStyle','-','LineWidth',lWidth);
    legend([choice_f_fMRI_hdl_mSplit.(low_mb_nm),...
        choice_f_fMRI_hdl_mSplit.(high_mb_nm)],...
        {low_mb_nm_bis, high_mb_nm_bis});
    legend('boxoff');
    ylim([0 100]);
    xlabel([fMRI_ROI_short_nm,' BOLD during ',timePeriod_nm,' - ',task_nm_tmp]);
    ylabel(['Choice (%) - ',task_nm_tmp]);
    legend_size(pSize);
    
    
    %% choice = f(ROI*E_level)
    fig;
    [E_hdl_mSplit, E_fit_hdl_mSplit] = deal(cell(n_E_levels*2, 1));
    for iM = 1:2
        switch iM
            case 1
                mb_nm = low_mb_nm;
                curve_col = low_mb_col;
            case 2
                mb_nm = high_mb_nm;
                curve_col = high_mb_col;
        end
        for iE = 1:n_E_levels
            jE = iE + n_E_levels*(iM-1);
            E_hdl_mSplit{jE} = errorbar(m_fMRI_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm)(:,iE),...
                m_choice_hE_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm)(:,iE).*100,...
                sem_choice_hE_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm)(:,iE).*100);
            E_fit_hdl_mSplit{jE} = plot(m_fMRI_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm)(:,iE),...
                m_choice_hE_fit_bins_mSplit.(task_nm_tmp).perElevel.(mb_nm)(:,iE).*100);
            E_hdl_mSplit{jE}.Marker = 'o';
            E_hdl_mSplit{jE}.LineWidth = lWidth;
            switch iM
                case 1
                    switch iE
                        case 1
                            E_hdl_mSplit{jE}.Color = E1_col_low;
                            E_fit_hdl_mSplit{jE}.Color = E1_col_low;
                        case 2
                            E_hdl_mSplit{jE}.Color = E2_col_low;
                            E_fit_hdl_mSplit{jE}.Color = E2_col_low;
                        case 3
                            E_hdl_mSplit{jE}.Color = E3_col_low;
                            E_fit_hdl_mSplit{jE}.Color = E3_col_low;
                    end % effort
                case 2
                    switch iE
                        case 1
                            E_hdl_mSplit{jE}.Color = E1_col_high;
                            E_fit_hdl_mSplit{jE}.Color = E1_col_high;
                        case 2
                            E_hdl_mSplit{jE}.Color = E2_col_high;
                            E_fit_hdl_mSplit{jE}.Color = E2_col_high;
                        case 3
                            E_hdl_mSplit{jE}.Color = E3_col_high;
                            E_fit_hdl_mSplit{jE}.Color = E3_col_high;
                    end % effort
            end % metabolite
            E_hdl_mSplit{jE}.LineStyle = 'none';
            E_fit_hdl_mSplit{jE}.LineWidth = lWidth;
            E_fit_hdl_mSplit{jE}.LineStyle = '--';
        end % effort level
    end % metabolite level
    line([0 0],[0 100],'Color','k','LineStyle','-','LineWidth',lWidth);
    ylim([0 100]);
    legend([E_hdl_mSplit{1},E_hdl_mSplit{2},E_hdl_mSplit{3},E_hdl_mSplit{4},E_hdl_mSplit{5},E_hdl_mSplit{6}],...
        {['E1 - low ',metabolite_nm],...
        ['E2 - low ',metabolite_nm],...
        ['E3 - low ',metabolite_nm],...
        ['E1 - high ',metabolite_nm],...
        ['E2 - high ',metabolite_nm],...
        ['E3 - high ',metabolite_nm]});
    legend('Location','SouthEast');
    legend('boxoff');
    xlabel([fMRI_ROI_short_nm,' BOLD during ',timePeriod_nm,' - ',task_nm_tmp]);
    ylabel(['Choice (%) - ',task_nm_tmp]);
    legend_size(pSize);

end % task loop