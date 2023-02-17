function[mediation_paths, pval] = RT_f_NV_uncertainty__mSplit_bhvPrm(study_nm, subject_id, condition, figDisp)
% [mediation_paths, pval] = RT_f_NV_uncertainty__mSplit_bhvPrm(study_nm, subject_id, condition, figDisp)
% RT_f_NV_uncertainty__mSplit_bhvPrm will look at how reaction times (RT) during choice
% vary with net value depending on the uncertainty rating and uncertainty
% derived from the model, looking at all tasks together and each task
% separately. Subjects will be split depending on their level of the
% selected parameter to see if there are differences based on the selected
% parameter.
%
% INPUTS
% study_nm: study name
%
% subject_id: list of subjects (can be left empty)
%
% condition: condition used to be defined
% (suggestion: you should rather use a non-saturation run filter to avoid
% runs where there is no variation in confidence ratings)
%
% figDisp: display figure (1) or not (0)
%
% OUTPUTS
% RT: structure with RT
%
% mediation_paths: structure with average (mean) and SEM for each path of
% the mediation for each task
%
% pval: structure with p.value for the t.test against zero of each path of
% the mediation for each task

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm) ||...
        ~ismember(study_nm,{'study1','study2'})
    study_nm = 'study1';
end
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    subject_id = LGCM_subject_selection(study_nm, condition);
end
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
end

%% general parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
prm_mSplit = {'low','high'};
nMsplit = length(prm_mSplit);

%% extract RT = f(NV)*uncertainty
[RT_f_pChoice_bin, pChoice_f_pChoice_bin, mediation_paths_all, pval_all,...
    RT_orth_ConfMdl_f_pChoice_bin, RT_orth_ConfMdl_f_RT_orth_ConfMdl_bin,...
    conf_rtg_f_RT_orth_ConfMdl_bin,...
    path_a_confMdl_RT, path_b_RT_confRtg,...
    path_c_confMdl_confRtg, path_cprime_confMdl_confRtg,...
    path_a_lowConfMdl_RT, path_b_lowConfMdl_RT_confRtg,...
    path_c_lowConfMdl_confRtg, path_cprime_lowConfMdl_confRtg] = RT_f_NV_uncertainty(study_nm, subject_id, condition, 0);

%% extract behavioural parameters
[low_prm_CID, high_prm_CID] = bhvPrm_mSplit(study_nm, subject_id);

%% average data across subjects, split per level of the parameter, split with median split
conditions = {'high_confRtg','low_confRtg','high_confMdl','low_confMdl'};
nCond = length(conditions);
for iMS = 1:nMsplit
    mSplit_nm = prm_mSplit{iMS};
    switch mSplit_nm
        case 'low'
            mSplit_idx = ismember(subject_id, low_prm_CID);
        case 'high'
            mSplit_idx = ismember(subject_id, high_prm_CID);
    end
    
    for iT = 1:nTasks
        task_nm = tasks{iT};
        for iC = 1:nCond
            cond_nm = conditions{iC};
            [m_RT_f_pChoice_bin.(task_nm).(cond_nm).(mSplit_nm),...
                sem_RT_f_pChoice_bin.(task_nm).(cond_nm).(mSplit_nm)] = mean_sem_sd(RT_f_pChoice_bin.(task_nm).(cond_nm)(:,mSplit_idx),2);
            [m_pChoice_f_pChoice_bin.(task_nm).(cond_nm).(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).(cond_nm).(mSplit_nm)] = mean_sem_sd(pChoice_f_pChoice_bin.(task_nm).(cond_nm)(:,mSplit_idx),2);
            switch cond_nm
                case {'high_confRtg','low_confRtg'}
                    [m_RT_orth_ConfMdl_f_pChoice_bin.(task_nm).(cond_nm).(mSplit_nm),...
                        sem_RT_orth_ConfMdl_f_pChoice_bin.(task_nm).(cond_nm).(mSplit_nm)] = mean_sem_sd(RT_orth_ConfMdl_f_pChoice_bin.(task_nm).(cond_nm)(:,mSplit_idx),2);
            end
        end % condition loop
        
        % average bin
        [m_conf_rtg_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm),...
            sem_conf_rtg_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm)] = mean_sem_sd(conf_rtg_f_RT_orth_ConfMdl_bin.(task_nm)(:,mSplit_idx),2);
        [m_RT_orth_ConfMdl_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm),...
            sem_RT_orth_ConfMdl_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm)] = mean_sem_sd(RT_orth_ConfMdl_f_RT_orth_ConfMdl_bin.(task_nm)(:,mSplit_idx),2);
        
        % average and test mediation paths
        [mediation_paths.(task_nm).mean.allTrials.a_confMdl_RT.(mSplit_nm),...
            mediation_paths.(task_nm).sem.allTrials.a_confMdl_RT.(mSplit_nm)] = mean_sem_sd(path_a_confMdl_RT.(task_nm)(:,mSplit_idx),2);
        [mediation_paths.(task_nm).mean.allTrials.b_RT_confRtg.(mSplit_nm),...
            mediation_paths.(task_nm).sem.allTrials.b_RT_confRtg.(mSplit_nm)] = mean_sem_sd(path_b_RT_confRtg.(task_nm)(:,mSplit_idx),2);
        [mediation_paths.(task_nm).mean.allTrials.c_confMdl_confRtg.(mSplit_nm),...
            mediation_paths.(task_nm).sem.allTrials.c_confMdl_confRtg.(mSplit_nm)] = mean_sem_sd(path_c_confMdl_confRtg.(task_nm)(:,mSplit_idx),2);
        [mediation_paths.(task_nm).mean.allTrials.cprime_confMdl_confRtg.(mSplit_nm),...
            mediation_paths.(task_nm).sem.allTrials.cprime_confMdl_confRtg.(mSplit_nm)] = mean_sem_sd(path_cprime_confMdl_confRtg.(task_nm)(:,mSplit_idx),2);
        [~,pval.(task_nm).allTrials.a_confMdl_RT.(mSplit_nm)] = ttest(path_a_confMdl_RT.(task_nm)(:,mSplit_idx));
        [~,pval.(task_nm).allTrials.b_RT_confRtg.(mSplit_nm)] = ttest(path_b_RT_confRtg.(task_nm)(:,mSplit_idx));
        [~,pval.(task_nm).allTrials.c_confMdl_confRtg.(mSplit_nm)] = ttest(path_c_confMdl_confRtg.(task_nm)(:,mSplit_idx));
        [~,pval.(task_nm).allTrials.cprime_confMdl_confRtg.(mSplit_nm)] = ttest(path_cprime_confMdl_confRtg.(task_nm)(:,mSplit_idx));
        
        % average and test mediation paths when focusing on low initial
        % confidence trials
        [mediation_paths.(task_nm).mean.lowConfMdl_trials.a_confMdl_RT.(mSplit_nm),...
            mediation_paths.(task_nm).sem.lowConfMdl_trials.a_confMdl_RT.(mSplit_nm)] = mean_sem_sd(path_a_lowConfMdl_RT.(task_nm)(:,mSplit_idx),2);
        [mediation_paths.(task_nm).mean.lowConfMdl_trials.b_RT_confRtg.(mSplit_nm),...
            mediation_paths.(task_nm).sem.lowConfMdl_trials.b_RT_confRtg.(mSplit_nm)] = mean_sem_sd(path_b_lowConfMdl_RT_confRtg.(task_nm)(:,mSplit_idx),2);
        [mediation_paths.(task_nm).mean.lowConfMdl_trials.c_confMdl_confRtg.(mSplit_nm),...
            mediation_paths.(task_nm).sem.lowConfMdl_trials.c_confMdl_confRtg.(mSplit_nm)] = mean_sem_sd(path_c_lowConfMdl_confRtg.(task_nm)(:,mSplit_idx),2);
        [mediation_paths.(task_nm).mean.lowConfMdl_trials.cprime_confMdl_confRtg.(mSplit_nm),...
            mediation_paths.(task_nm).sem.lowConfMdl_trials.cprime_confMdl_confRtg.(mSplit_nm)] = mean_sem_sd(path_cprime_lowConfMdl_confRtg.(task_nm)(:,mSplit_idx),2);
        [~,pval.(task_nm).lowConfMdl_trials.a_confMdl_RT.(mSplit_nm)] = ttest(path_a_lowConfMdl_RT.(task_nm)(:,mSplit_idx));
        [~,pval.(task_nm).lowConfMdl_trials.b_RT_confRtg.(mSplit_nm)] = ttest(path_b_lowConfMdl_RT_confRtg.(task_nm)(:,mSplit_idx));
        [~,pval.(task_nm).lowConfMdl_trials.c_confMdl_confRtg.(mSplit_nm)] = ttest(path_c_lowConfMdl_confRtg.(task_nm)(:,mSplit_idx));
        [~,pval.(task_nm).lowConfMdl_trials.cprime_confMdl_confRtg.(mSplit_nm)] = ttest(path_cprime_lowConfMdl_confRtg.(task_nm)(:,mSplit_idx));
    end % task loop
end % median split

%% display figure
if figDisp == 1
    lWidth = 3;
    pSize = 30;
    lowConf_col = [5 113 176]./255;
    highConf_col = [202 0 32]./255;

    for iT = 1:nTasks
        task_nm = tasks{iT};
        
        %% RT = f(confidence)
        fig;
        %% RT = f(NV) split per confidence based on the model
        subplot(1,3,1);
        hold on;
        
        for iMS = 1:nMsplit
            mSplit_nm = prm_mSplit{iMS};
            
            % low confidence
            lowConfMdl_hdl.(mSplit_nm) = errorbar(m_pChoice_f_pChoice_bin.(task_nm).low_confMdl.(mSplit_nm),...
                m_RT_f_pChoice_bin.(task_nm).low_confMdl.(mSplit_nm),...
                sem_RT_f_pChoice_bin.(task_nm).low_confMdl.(mSplit_nm),...
                sem_RT_f_pChoice_bin.(task_nm).low_confMdl.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).low_confMdl.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).low_confMdl.(mSplit_nm));
            lowConfMdl_hdl.(mSplit_nm).LineWidth = lWidth;
            lowConfMdl_hdl.(mSplit_nm).Color = lowConf_col;
            
            % high confidence
            highConfMdl_hdl.(mSplit_nm) = errorbar(m_pChoice_f_pChoice_bin.(task_nm).high_confMdl.(mSplit_nm),...
                m_RT_f_pChoice_bin.(task_nm).high_confMdl.(mSplit_nm),...
                sem_RT_f_pChoice_bin.(task_nm).high_confMdl.(mSplit_nm),...
                sem_RT_f_pChoice_bin.(task_nm).high_confMdl.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).high_confMdl.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).high_confMdl.(mSplit_nm));
            highConfMdl_hdl.(mSplit_nm).LineWidth = lWidth;
            highConfMdl_hdl.(mSplit_nm).Color = highConf_col;
            
            switch mSplit_nm
                case 'low'
                    lowConfMdl_hdl.low.LineStyle = '--';
                    highConfMdl_hdl.low.LineStyle = '--';
                case 'high'
                    lowConfMdl_hdl.high.LineStyle = '-';
                    highConfMdl_hdl.high.LineStyle = '-';
            end
        end
        legend([highConfMdl_hdl.high, lowConfMdl_hdl.high,...
            highConfMdl_hdl.low, lowConfMdl_hdl.low],...
            {'high model Conf', 'low model Conf',...
            'high model Conf', 'low model Conf'});
        legend('boxoff');
        legend('Location','South');
        xlabel(['p(high effort) (',task_nm,')']);
        ylabel('RT (s)');
        legend_size(pSize);
        
        %% RT = f(NV) split per confidence rating
        subplot(1,3,2);
        hold on;
        
        for iMS = 1:nMsplit
            mSplit_nm = prm_mSplit{iMS};
            % low confidence
            lowConfRtg_hdl.(mSplit_nm) = errorbar(m_pChoice_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                m_RT_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                sem_RT_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                sem_RT_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm));
            lowConfRtg_hdl.(mSplit_nm).LineWidth = lWidth;
            lowConfRtg_hdl.(mSplit_nm).Color = lowConf_col;
            
            % high confidence
            highConfRtg_hdl.(mSplit_nm) = errorbar(m_pChoice_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                m_RT_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                sem_RT_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                sem_RT_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm));
            highConfRtg_hdl.(mSplit_nm).LineWidth = lWidth;
            highConfRtg_hdl.(mSplit_nm).Color = highConf_col;
            
            switch mSplit_nm
                case 'low'
                    lowConfRtg_hdl.low.LineStyle = '--';
                    highConfRtg_hdl.low.LineStyle = '--';
                case 'high'
                    lowConfRtg_hdl.high.LineStyle = '-';
                    highConfRtg_hdl.high.LineStyle = '-';
            end
        end
        legend([highConfRtg_hdl.low, lowConfRtg_hdl.low,...
            highConfRtg_hdl.high, lowConfRtg_hdl.high],...
            {'high rating Conf', 'low rating Conf',...
            'high rating Conf', 'low rating Conf'});
        legend('boxoff');
        legend('Location','South');
        xlabel(['p(high effort) - ',task_nm]);
        ylabel('RT (s)');
        legend_size(pSize);
        
        %% RT (orthogonalized to model confidence) = f(NV) split per confidence rating
        subplot(1,3,3);
        hold on;
        
        for iMS = 1:nMsplit
            mSplit_nm = prm_mSplit{iMS};
            % low confidence
            lowConfRtg_RTorth_hdl.(mSplit_nm) = errorbar(m_pChoice_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                m_RT_orth_ConfMdl_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                sem_RT_orth_ConfMdl_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                sem_RT_orth_ConfMdl_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).low_confRtg.(mSplit_nm));
            lowConfRtg_RTorth_hdl.(mSplit_nm).LineWidth = lWidth;
            lowConfRtg_RTorth_hdl.(mSplit_nm).Color = lowConf_col;
            
            % high confidence
            highConfRtg_RTorth_hdl.(mSplit_nm) = errorbar(m_pChoice_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                m_RT_orth_ConfMdl_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                sem_RT_orth_ConfMdl_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                sem_RT_orth_ConfMdl_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm),...
                sem_pChoice_f_pChoice_bin.(task_nm).high_confRtg.(mSplit_nm));
            highConfRtg_RTorth_hdl.(mSplit_nm).LineWidth = lWidth;
            highConfRtg_RTorth_hdl.(mSplit_nm).Color = highConf_col;
            
            switch mSplit_nm
                case 'low'
                    lowConfRtg_RTorth_hdl.low.LineStyle = '--';
                    highConfRtg_RTorth_hdl.low.LineStyle = '--';
                case 'high'
                    lowConfRtg_RTorth_hdl.high.LineStyle = '-';
                    highConfRtg_RTorth_hdl.high.LineStyle = '-';
            end
        end % median split loop
        legend([highConfRtg_RTorth_hdl.low, lowConfRtg_RTorth_hdl.low,...
            highConfRtg_RTorth_hdl.high, lowConfRtg_RTorth_hdl.high],...
            {'high rating Conf', 'low rating Conf',...
            'high rating Conf', 'low rating Conf'});
        legend('boxoff');
        legend('Location','South');
        xlabel(['p(high effort) - ',task_nm]);
        ylabel('RT (s) (orthogonalized to model conf.)');
        legend_size(pSize);

        
        %% confidence rating = f(RT orthogonalized to confidence model)
        fig;
        for iMS = 1:nMsplit
            mSplit_nm = prm_mSplit{iMS};
            conf_rtg_RTorth_hdl.(mSplit_nm) = errorbar(m_RT_orth_ConfMdl_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm),...
                m_conf_rtg_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm),...
                sem_conf_rtg_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm),...
                sem_conf_rtg_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm),...
                sem_RT_orth_ConfMdl_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm),...
                sem_RT_orth_ConfMdl_f_RT_orth_ConfMdl_bin.(task_nm).(mSplit_nm));
            conf_rtg_RTorth_hdl.(mSplit_nm).LineWidth = lWidth;
            conf_rtg_RTorth_hdl.(mSplit_nm).Color = 'k';
            switch mSplit_nm
                case 'low'
                    conf_rtg_RTorth_hdl.low.LineStyle = '--';
                case 'high'
                    conf_rtg_RTorth_hdl.high.LineStyle = '-';
            end
        end % median split loop
        xlabel('RT (s) (orthogonalized to model conf.)');
        ylabel('Confidence rating');
        legend_size(pSize);
    end % task loop
end % figure display
end % function