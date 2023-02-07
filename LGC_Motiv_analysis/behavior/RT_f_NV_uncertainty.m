function[RT, b_RT, pval] = RT_f_NV_uncertainty(study_nm, subject_id, condition, figDisp)
% [RT, b_RT, pval] = RT_f_NV_uncertainty(study_nm, subject_id, condition, figDisp)
% RT_f_NV_uncertainty will look at how reaction times (RT) during choice
% vary with net value depending on the uncertainty rating and uncertainty
% derived from the model, looking at all tasks together and each task
% separately.
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
% b_RT: beta for linear regressions
%
% pval: p.value for different linear regressions

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm) ||...
        ~ismember(study_nm,{'study1','study2'})
    study_nm = 'study1';
end
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
end

%% working directories
computerRoot = LGCM_root_paths;
dataRoot = [computerRoot, filesep, study_nm, filesep];
gitResultsFolder = [fullfile('C:','Users','clairis','Desktop',...
    'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];

%% main parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
nRunsPerTask = 2;
mdl_nm = 'mdl_3';
nBins = 6;
[RT_f_NV_bin.Ep.high_confRtg, NV_f_NV_bin.Ep.high_confRtg,...
    RT_f_NV_bin.Ep.low_confRtg, NV_f_NV_bin.Ep.low_confRtg,...
    RT_f_NV_bin.Em.high_confRtg, NV_f_NV_bin.Em.high_confRtg,...
    RT_f_NV_bin.Em.low_confRtg, NV_f_NV_bin.Em.low_confRtg,...
    RT_f_NV_bin.Ep.high_confMdl, NV_f_NV_bin.Ep.high_confMdl,...
    RT_f_NV_bin.Ep.low_confMdl, NV_f_NV_bin.Ep.low_confMdl,...
    RT_f_NV_bin.Em.high_confMdl, NV_f_NV_bin.Em.high_confMdl,...
    RT_f_NV_bin.Em.low_confMdl, NV_f_NV_bin.Em.low_confMdl] = deal(NaN(nBins,NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [dataRoot, filesep, 'CID',sub_nm, filesep, 'behavior',filesep];
    runs = runs_definition(study_nm, sub_nm, condition);

    for iTask = 1:nTasks
        task_nm = tasks{iTask};
        switch task_nm
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end

        % re-initialize data for each subject
        [RT_f_NV_bin_highConfRtg_tmp, RT_f_NV_bin_lowConfRtg_tmp,...
            NV_f_NV_bin_highConfRtg_tmp, NV_f_NV_bin_lowConfRtg_tmp,...
            RT_f_NV_bin_highConfMdl_tmp, RT_f_NV_bin_lowConfMdl_tmp,...
            NV_f_NV_bin_highConfMdl_tmp, NV_f_NV_bin_lowConfMdl_tmp] = deal(NaN(nBins, nRunsPerTask));

        for iRun = 1:runs.nb_runs.(task_nm)
            jRun = runs.(task_nm).runsToKeep(iRun);
            run_nm = num2str(jRun);
            % extract RT, conf from the model, conf rating and NV
            [~, deltaNV_tmp, modelConf_tmp] = extract_bayesian_mdl(gitResultsFolder, subBehaviorFolder,...
                sub_nm, run_nm, task_fullName, mdl_nm);
            % RT
            [RT_tmp] = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            % confidence rating
            [conf_rtg_tmp] = extract_confidence_rating(subBehaviorFolder, sub_nm, run_nm, task_fullName);

            % split data per confidence  levels
            high_confRtg = conf_rtg_tmp == 1;
            low_confRtg = conf_rtg_tmp == 0;
            high_confMdl = modelConf_tmp > 0.5;
            low_confMdl = modelConf_tmp <= 0.5;

            % then do bins
            % bins for rated confidence
            if sum(high_confRtg) > 0
                [RT_f_NV_bin_highConfRtg_tmp(:,iRun), NV_f_NV_bin_highConfRtg_tmp(:,iRun)] = do_bin2(RT_tmp(high_confRtg), deltaNV_tmp(high_confRtg), nBins, 0);
            end
            if sum(low_confRtg) > 0
                [RT_f_NV_bin_lowConfRtg_tmp(:,iRun), NV_f_NV_bin_lowConfRtg_tmp(:,iRun)] = do_bin2(RT_tmp(low_confRtg), deltaNV_tmp(low_confRtg), nBins, 0);
            end
            % bins for confidence inferred by the model
            if sum(high_confMdl) > 0
                [RT_f_NV_bin_highConfMdl_tmp(:,iRun), NV_f_NV_bin_highConfMdl_tmp(:,iRun)] = do_bin2(RT_tmp(high_confMdl), deltaNV_tmp(high_confMdl), nBins, 0);
            end
            if sum(low_confMdl) > 0
                [RT_f_NV_bin_lowConfMdl_tmp(:,iRun), NV_f_NV_bin_lowConfMdl_tmp(:,iRun)] = do_bin2(RT_tmp(low_confMdl), deltaNV_tmp(low_confMdl), nBins, 0);
            end
        end % run loop

        %% average bins RT=f(NV) for each group (low/high uncertainty) across runs
        % same for confidence derived from the ratings
        RT_f_NV_bin.(task_nm).high_confRtg(:,iS) = mean(RT_f_NV_bin_highConfRtg_tmp,2,'omitnan');
        NV_f_NV_bin.(task_nm).high_confRtg(:,iS) = mean(NV_f_NV_bin_highConfRtg_tmp,2,'omitnan');
        RT_f_NV_bin.(task_nm).low_confRtg(:,iS) = mean(RT_f_NV_bin_lowConfRtg_tmp,2,'omitnan');
        NV_f_NV_bin.(task_nm).low_confRtg(:,iS) = mean(NV_f_NV_bin_lowConfRtg_tmp,2,'omitnan');
        % same for confidence derived from the model
        RT_f_NV_bin.(task_nm).high_confMdl(:,iS) = mean(RT_f_NV_bin_highConfMdl_tmp,2,'omitnan');
        NV_f_NV_bin.(task_nm).high_confMdl(:,iS) = mean(NV_f_NV_bin_highConfMdl_tmp,2,'omitnan');
        RT_f_NV_bin.(task_nm).low_confMdl(:,iS) = mean(RT_f_NV_bin_lowConfMdl_tmp,2,'omitnan');
        NV_f_NV_bin.(task_nm).low_confMdl(:,iS) = mean(NV_f_NV_bin_lowConfMdl_tmp,2,'omitnan');
    end % task loop
end % subject loop
%% average data across subjects
conditions = {'high_confRtg','low_confRtg','high_confMdl','low_confMdl'};
nCond = length(conditions);
for iT = 1:nTasks
    task_nm = tasks{iT};
    for iC = 1:nCond
        cond_nm = conditions{iC};
        [m_RT_f_NV_bin.(task_nm).(cond_nm),...
            sem_RT_f_NV_bin.(task_nm).(cond_nm)] = mean_sem_sd(RT_f_NV_bin.(task_nm).(cond_nm),2,'omitnan');
        [m_NV_f_NV_bin.(task_nm).(cond_nm),...
            sem_NV_f_NV_bin.(task_nm).(cond_nm)] = mean_sem_sd(NV_f_NV_bin.(task_nm).(cond_nm),2,'omitnan');
    end % condition loop
end % task loop

%% display figure
if figDisp == 1
    lWidth = 3;
    pSize = 40;
    lowConf_col = [5 113 176]./255;
    highConf_col = [202 0 32]./255;

    for iT = 1:nTasks
        task_nm = tasks{iT};

        fig;
        % RT = f(NV) split per confidence rating
        subplot(1,2,1);
        % low confidence
        lowConfRtg_hdl = errorbar(m_NV_f_NV_bin.(task_nm).low_confRtg,...
            m_RT_f_NV_bin.(task_nm).low_confRtg,...
            sem_RT_f_NV_bin.(task_nm).low_confRtg,...
            sem_RT_f_NV_bin.(task_nm).low_confRtg,...
            sem_NV_f_NV_bin.(task_nm).low_confRtg,...
            sem_NV_f_NV_bin.(task_nm).low_confRtg);
        lowConfRtg_hdl.LineWidth = lWidth;
        lowConfRtg_hdl.Color = lowConf_col;
        % high confidence
        highConfRtg_hdl = errorbar(m_NV_f_NV_bin.(task_nm).high_confRtg,...
            m_RT_f_NV_bin.(task_nm).high_confRtg,...
            sem_RT_f_NV_bin.(task_nm).high_confRtg,...
            sem_RT_f_NV_bin.(task_nm).high_confRtg,...
            sem_NV_f_NV_bin.(task_nm).high_confRtg,...
            sem_NV_f_NV_bin.(task_nm).high_confRtg);
        highConfRtg_hdl.LineWidth = lWidth;
        highConfRtg_hdl.Color = highConf_col;
        legend([highConfRtg_hdl, lowConfRtg_hdl],...
            {'high confidence rating', 'low confidence rating'});
        legend('boxoff');
        xlabel(['NV - ',task_nm]);
        ylabel('RT (s)');
        legend_size(pSize);

        % RT = f(NV) split per confidence based on the model
        subplot(1,2,2);
        % low confidence
        lowConfMdl_hdl = errorbar(m_NV_f_NV_bin.(task_nm).low_confMdl,...
            m_RT_f_NV_bin.(task_nm).low_confMdl,...
            sem_RT_f_NV_bin.(task_nm).low_confMdl,...
            sem_RT_f_NV_bin.(task_nm).low_confMdl,...
            sem_NV_f_NV_bin.(task_nm).low_confMdl,...
            sem_NV_f_NV_bin.(task_nm).low_confMdl);
        lowConfMdl_hdl.LineWidth = lWidth;
        lowConfMdl_hdl.Color = lowConf_col;
        % high confidence
        highConfMdl_hdl = errorbar(m_NV_f_NV_bin.(task_nm).high_confMdl,...
            m_RT_f_NV_bin.(task_nm).high_confMdl,...
            sem_RT_f_NV_bin.(task_nm).high_confMdl,...
            sem_RT_f_NV_bin.(task_nm).high_confMdl,...
            sem_NV_f_NV_bin.(task_nm).high_confMdl,...
            sem_NV_f_NV_bin.(task_nm).high_confMdl);
        highConfMdl_hdl.LineWidth = lWidth;
        highConfMdl_hdl.Color = highConf_col;
        legend([highConfRtg_hdl, lowConfRtg_hdl],...
            {'high model confidence', 'low model confidence'});
        legend('boxoff');
        xlabel(['NV - ',task_nm]);
        ylabel('RT (s)');
        legend_size(pSize);

    end % task loop
end % figure display
end % function