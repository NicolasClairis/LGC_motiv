function[b_choice_f_fMRI, pval, fMRI_bins,...
    choice_hE_bins, choice_hE_fit_bins,...
    fMRI_ROI_short_nm, timePeriod_nm] = choice_f_ROI(nBins,...
    study_nm, subject_id, condition,...
    RT_orth,...
    figDisp)
% [b_choice_f_fMRI, pval, fMRI_bins,...
%     choice_hE_bins, choice_hE_fit_bins,...
%     fMRI_ROI_short_nm, timePeriod_nm] = choice_f_ROI(nBins,...
%     study_nm, subject_id, condition,...
%     RT_orth,...
%     figDisp)
% choice_f_ROI will look at percentage of choices depending on ROI level of
% activity. Then it will also split the data depending on the effort level
% proposed and try to look whether the results vary accordingly.
%
% INPUTS
% nBins : number of bins
%
% study_nm: study name ('study1'/'study2')
%
% subject_id: list of subjects to consider
%
% condition: condition to use for the extraction (which subjects and runs
% to include)
%
% RT_orth: do you want to orthogonalize the fMRI activity to reaction times
% (RT_orth=1) or not (RT_orth=0)? by default, no orthogonalization is
% applied
%
% figDisp: display figure (1) or not (0)?
%
% OUTPUTS
% b_choice_f_fMRI: structure with beta for the slope of choice = f(ROI)
%
% pval: p.value for t.test against zero of b_choice_f_fMRI
%
% fMRI_bins: structure with ROI activity binned
%
% choice_hE_bins: structure with choices binned depending on ROI level of
% activity and effort level
%
% choice_hE_fit_bins: structure with fit of the choices=f(ROI)
%
% fMRI_ROI_short_nm: fMRI ROI short name
%
% timePeriod_nm: time period of fMRI to look at

%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% subject selection
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end

%% general parameters
if ~exist('nBins','var') || isempty(nBins)
    nBins = 6;
end
if ~exist('RT_orth','var') || isempty(RT_orth)
    RT_orth = 0;
end
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
end
tasks = {'Ep','Em'};
nTasks = length(tasks);
nTrialsPerRun = 54;
nRunsPerTask = 2;
nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
nTrials = nTrialsPerRun*nRunsPerTask*nTasks;
n_E_levels = 3;
[fMRI_allTrials.Ep, choice_hE_allTrials.Ep, E_level.Ep, choice_hE_fit_allTrials.Ep, RT_allTrials.Ep,...
    fMRI_allTrials.Em, choice_hE_allTrials.Em, E_level.Em, choice_hE_fit_allTrials.Em, RT_allTrials.Em] = deal(NaN(nTrialsPerTask, NS));
[fMRI_allTrials.EpEmPool, choice_hE_allTrials.EpEmPool, E_level.EpEmPool,...
    choice_hE_fit_allTrials.EpEmPool, RT_allTrials.EpEmPool] = deal(NaN(nTrials, NS));
[choice_hE_fit_perElevel.Ep, choice_hE_fit_perElevel.Em] = deal(NaN(nTrialsPerTask/n_E_levels, n_E_levels, NS));
choice_hE_fit_perElevel.EpEmPool = NaN(nTrials/n_E_levels, n_E_levels, NS);
[fMRI_bins.Ep.allTrials, choice_hE_bins.Ep.allTrials, choice_hE_fit_bins.Ep.allTrials,...
    fMRI_bins.Em.allTrials, choice_hE_bins.Em.allTrials, choice_hE_fit_bins.Em.allTrials,...
    fMRI_bins.EpEmPool.allTrials, choice_hE_bins.EpEmPool.allTrials, choice_hE_fit_bins.EpEmPool.allTrials] = deal(NaN(nBins, NS));
[b_choice_f_fMRI.Ep.allTrials, b_choice_f_fMRI.Em.allTrials, b_choice_f_fMRI.EpEmPool.allTrials] = deal(NaN(2,NS));
[fMRI_bins.Ep.perElevel, choice_hE_bins.Ep.perElevel, choice_hE_fit_bins.Ep.perElevel,...
    fMRI_bins.Em.perElevel, choice_hE_bins.Em.perElevel, choice_hE_fit_bins.Em.perElevel,...
    fMRI_bins.EpEmPool.perElevel, choice_hE_bins.EpEmPool.perElevel, choice_hE_fit_bins.EpEmPool.perElevel] = deal(NaN(nBins, n_E_levels, NS));
[b_choice_f_fMRI.Ep.perElevel, b_choice_f_fMRI.Em.perElevel, b_choice_f_fMRI.EpEmPool.perElevel] = deal(NaN(2, NS, n_E_levels));

%% extract ROI activity for all subjects
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);
% define which ROI, and which time period is of interest to you
% define ROI
[fMRI_ROI_nm, fMRI_ROI_short_nm,...
    ~,...
    timePeriod_nm] = extract_ROI_betas_onsets_only_questInfos(ROI_trial_b_trial);

%% fit parameters
fitType = 'binomial'; % binomial or normal
switch fitType
    case 'binomial'
        fitType_glmval = 'logit';
    case 'normal'
        fitType_glmval = 'identity';
end

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
%     if ~ismember(sub_nm,{'054','061'})
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
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
        fMRI_allTrials.(task_nm_tmp)(runTrials_idx, iS) = ROI_trial_b_trial.(fMRI_ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
        
        %% same but for Ep+Em pool
        choice_hE_allTrials.EpEmPool(runTrials_EpEmPool_idx, iS) = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        E_level.EpEmPool(runTrials_EpEmPool_idx, iS) = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        RT_allTrials.EpEmPool(runTrials_EpEmPool_idx, iS) = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        fMRI_allTrials.EpEmPool(runTrials_EpEmPool_idx, iS) = ROI_trial_b_trial.(fMRI_ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
    end % run loop

    %% remove RT effect from fMRI ROI
    if RT_orth == 1
        % estimate how ROI is impacted by RT
        beta_Ep_tmp = glmfit(RT_allTrials.Ep(:,iS), fMRI_allTrials.Ep(:,iS),'normal');
        beta_Em_tmp = glmfit(RT_allTrials.Em(:,iS), fMRI_allTrials.Em(:,iS),'normal');
        beta_EpEmPool_tmp = glmfit(RT_allTrials.EpEmPool(:,iS), fMRI_allTrials.EpEmPool(:,iS),'normal');
        % then orthogonalize ROI activity to RT
        fMRI_allTrials.Ep(:,iS) = fMRI_allTrials.Ep(:,iS) - beta_Ep_tmp(2).*RT_allTrials.Ep(:,iS);
        fMRI_allTrials.Em(:,iS) = fMRI_allTrials.Em(:,iS) - beta_Em_tmp(2).*RT_allTrials.Em(:,iS);
        fMRI_allTrials.EpEmPool(:,iS) = fMRI_allTrials.EpEmPool(:,iS) - beta_EpEmPool_tmp(2).*RT_allTrials.EpEmPool(:,iS);
    end 
    %% fit
    for iT = 1:(nTasks+1)
        switch iT
            case {1,2}
                task_nm_tmp = tasks{iT};
            case 3
                task_nm_tmp = 'EpEmPool';
        end
        % choice = f(fMRI) logit
        b_choice_f_fMRI.(task_nm_tmp).allTrials(:,iS) = glmfit(fMRI_allTrials.(task_nm_tmp)(:, iS),...
            choice_hE_allTrials.(task_nm_tmp)(:, iS),fitType);
        choice_hE_fit_allTrials.(task_nm_tmp)(:, iS) = glmval(b_choice_f_fMRI.(task_nm_tmp).allTrials(:,iS),...
            fMRI_allTrials.(task_nm_tmp)(:, iS), fitType_glmval);
        
        % choice = f(fMRI) logit split per effort level
        for iE = 1:n_E_levels
            E_lvl_idx = E_level.(task_nm_tmp)(:, iS) == iE;
            b_choice_f_fMRI.(task_nm_tmp).perElevel(:,iS, iE) = glmfit(fMRI_allTrials.(task_nm_tmp)(E_lvl_idx, iS),...
                choice_hE_allTrials.(task_nm_tmp)(E_lvl_idx, iS),fitType);
            if (sum(E_lvl_idx) == nTrialsPerTask/n_E_levels && ismember(task_nm_tmp,{'Ep','Em'})) ||...
                    (sum(E_lvl_idx) == nTrials/n_E_levels && strcmp(task_nm_tmp,'EpEmPool'))
                choice_hE_fit_perElevel.(task_nm_tmp)(:, iE, iS) = glmval(b_choice_f_fMRI.(task_nm_tmp).perElevel(:,iS,iE),...
                    fMRI_allTrials.(task_nm_tmp)(E_lvl_idx, iS), fitType_glmval);
            else % prevent bug from subjects like CID040 where some runs were not performed in fMRI
                choice_hE_fit_perElevel.(task_nm_tmp)(1:sum(E_lvl_idx), iE, iS) = glmval(b_choice_f_fMRI.(task_nm_tmp).perElevel(:, iS, iE),...
                    fMRI_allTrials.(task_nm_tmp)(E_lvl_idx, iS), fitType_glmval);
            end
        end % effort loop
        
        %% extract the bins
        % all trials
        [choice_hE_bins.(task_nm_tmp).allTrials(:,iS),...
            fMRI_bins.(task_nm_tmp).allTrials(:,iS)] = do_bin2(choice_hE_allTrials.(task_nm_tmp)(:, iS),...
            fMRI_allTrials.(task_nm_tmp)(:, iS), nBins, 0);
        [choice_hE_fit_bins.(task_nm_tmp).allTrials(:,iS),...
            fMRI_bins.(task_nm_tmp).allTrials(:,iS)] = do_bin2(choice_hE_fit_allTrials.(task_nm_tmp)(:, iS),...
            fMRI_allTrials.(task_nm_tmp)(:, iS), nBins, 0);
        
        % split per effort level
        for iE = 1:n_E_levels
            E_lvl_idx = E_level.(task_nm_tmp)(:, iS) == iE;
            if sum(E_lvl_idx) > 0
                [choice_hE_bins.(task_nm_tmp).perElevel(:,iE,iS),...
                    fMRI_bins.(task_nm_tmp).perElevel(:,iE,iS)] = do_bin2(choice_hE_allTrials.(task_nm_tmp)(E_lvl_idx, iS),...
                    fMRI_allTrials.(task_nm_tmp)(E_lvl_idx, iS), nBins, 0);
                [choice_hE_fit_bins.(task_nm_tmp).perElevel(:,iE,iS),...
                    fMRI_bins.(task_nm_tmp).perElevel(:,iE,iS)] = do_bin2(choice_hE_fit_perElevel.(task_nm_tmp)(1:sum(E_lvl_idx), iE, iS),...
                    fMRI_allTrials.(task_nm_tmp)(E_lvl_idx, iS), nBins, 0);
            end
        end
    end % task loop
%     end % bad subject Em filter
end % subject loop

if figDisp == 1
    pSize = 30;
    lWidth = 3;
    E1_col = [229 245 224]./255;
    E2_col = [161 217 155]./255;
    E3_col = [49 163 84]./255;
end
for iT = 1:(nTasks+1)
    switch iT
        case {1,2}
            task_nm_tmp = tasks{iT};
        case 3
            task_nm_tmp = 'EpEmPool';
    end
    
    %% average data across subjects
    % all trials
    [m_fMRI_bins.(task_nm_tmp).allTrials,...
        sem_fMRI_bins.(task_nm_tmp).allTrials] = mean_sem_sd(fMRI_bins.(task_nm_tmp).allTrials,2);
    [m_choice_hE_bins.(task_nm_tmp).allTrials,...
        sem_choice_hE_bins.(task_nm_tmp).allTrials] = mean_sem_sd(choice_hE_bins.(task_nm_tmp).allTrials,2);
    [m_choice_hE_fit_bins.(task_nm_tmp).allTrials,...
        sem_choice_hE_fit_bins.(task_nm_tmp).allTrials] = mean_sem_sd(choice_hE_fit_bins.(task_nm_tmp).allTrials,2);
    [~,pval.(task_nm_tmp).choice_f_fMRI.allTrials] = ttest(b_choice_f_fMRI.(task_nm_tmp).allTrials(2,:));
    pval.(task_nm_tmp).choice_f_fMRI.perE = NaN(n_E_levels,1);
    for iE = 1:n_E_levels
        [~,pval.(task_nm_tmp).choice_f_fMRI.perE(iE)] = ttest(b_choice_f_fMRI.(task_nm_tmp).perElevel(2,:,iE));
    end
    % split per E level
    [m_fMRI_bins.(task_nm_tmp).perElevel,...
        sem_fMRI_bins.(task_nm_tmp).perElevel] = mean_sem_sd(fMRI_bins.(task_nm_tmp).perElevel,3);
    [m_choice_hE_bins.(task_nm_tmp).perElevel,...
        sem_choice_hE_bins.(task_nm_tmp).perElevel] = mean_sem_sd(choice_hE_bins.(task_nm_tmp).perElevel,3);
    [m_choice_hE_fit_bins.(task_nm_tmp).perElevel,...
        sem_choice_hE_fit_bins.(task_nm_tmp).perElevel] = mean_sem_sd(choice_hE_fit_bins.(task_nm_tmp).perElevel,3);
    
    %% figure
    if figDisp == 1
        %% choice = f(ROI)
        fig;
        choice_f_fMRI_hdl = errorbar(m_fMRI_bins.(task_nm_tmp).allTrials,...
            m_choice_hE_bins.(task_nm_tmp).allTrials.*100,...
            sem_choice_hE_bins.(task_nm_tmp).allTrials.*100);
        choice_f_fMRI_hdl.LineWidth = lWidth;
        choice_f_fMRI_hdl.Marker = 'o';
        choice_f_fMRI_hdl.Color = 'k';
        choice_f_fMRI_hdl.LineStyle = 'none';
        choice_fit_f_fMRI_hdl = plot(m_fMRI_bins.(task_nm_tmp).allTrials,...
            m_choice_hE_fit_bins.(task_nm_tmp).allTrials.*100);
        choice_fit_f_fMRI_hdl.LineWidth = lWidth;
        choice_fit_f_fMRI_hdl.Color = 'k';
        choice_fit_f_fMRI_hdl.LineStyle = '--';
        line([0 0],[0 100],'Color','k','LineStyle','-','LineWidth',lWidth);
        ylim([30 80]);
        xlabel([fMRI_ROI_short_nm,' BOLD during ',timePeriod_nm,' - ',task_nm_tmp]);
        ylabel(['Choice (%) - ',task_nm_tmp]);
        legend_size(pSize);
        
        %% choice = f(ROI*E_level)
        fig;
        [E_hdl, E_fit_hdl] = deal(cell(n_E_levels, 1));
        for iE = 1:n_E_levels
            E_hdl{iE} = errorbar(m_fMRI_bins.(task_nm_tmp).perElevel(:,iE),...
                m_choice_hE_bins.(task_nm_tmp).perElevel(:,iE).*100,...
                sem_choice_hE_bins.(task_nm_tmp).perElevel(:,iE).*100);
            E_fit_hdl{iE} = plot(m_fMRI_bins.(task_nm_tmp).perElevel(:,iE),...
                m_choice_hE_fit_bins.(task_nm_tmp).perElevel(:,iE).*100);
            E_hdl{iE}.Marker = 'o';
            E_hdl{iE}.LineWidth = lWidth;
            switch iE
                case 1
                    E_hdl{iE}.Color = E1_col;
                    E_fit_hdl{iE}.Color = E1_col;
                case 2
                    E_hdl{iE}.Color = E2_col;
                    E_fit_hdl{iE}.Color = E2_col;
                case 3
                    E_hdl{iE}.Color = E3_col;
                    E_fit_hdl{iE}.Color = E3_col;
            end
            E_hdl{iE}.LineStyle = 'none';
            E_fit_hdl{iE}.LineWidth = lWidth;
            E_fit_hdl{iE}.LineStyle = '--';
        end % effort level
        line([0 0],[0 100],'Color','k','LineStyle','-','LineWidth',lWidth);
        ylim([0 100]);
        legend([E_hdl{1},E_hdl{2},E_hdl{3}],...
            {'E1','E2','E3'});
        legend('Location','SouthEast');
        legend('boxoff');
        xlabel([fMRI_ROI_short_nm,' BOLD during ',timePeriod_nm,' - ',task_nm_tmp]);
        ylabel(['Choice (%) - ',task_nm_tmp]);
        legend_size(pSize);
    end % figure display
end % task loop
end % function