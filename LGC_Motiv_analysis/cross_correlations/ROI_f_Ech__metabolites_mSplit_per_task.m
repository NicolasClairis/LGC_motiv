%% extract the level of activity of a given ROI in function of each effort level chosen
% and split subjects with a median split based on the metabolites. the same
% script will also look at the choice depending on dmPFC activity and
% effort level proposed for the high effort option
%
% Same as ROI_f_Ech__metabolites_mSplit_per_task.m but spliting data
% between mental and physical effort task

%% general parameters
% display figures?
dispFig = true;

%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directory
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% subject definition
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load ROI
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);
% define which ROI, and which time period is of interest to you
% define ROI
ROI_names = fieldnames(ROI_trial_b_trial);
ROI_subList = ROI_trial_b_trial.subject_id;
ROI_names(strcmp(ROI_names,'subject_id')) = [];
if length(ROI_names) > 1
    error(['There should be only 1 ROI selected, not ',num2str(length(ROI_names))])
else
    ROI_nm = ROI_names;
end
ROI_short_nm = inputdlg('ROI short name?');
ROI_short_nm = ROI_short_nm{1};
% define task
task_names = {'Ep','Em'};
nTasks = length(task_names);
% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};

%% load behavior 
nRunsPerTask = 2;
nTrialsPerRun = 54;
nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
% main variables of interest
Echosen_possible = 0:3;
n_Ech = length(Echosen_possible);
hE_levels = 1:3;
n_hE_levels = length(hE_levels);
% data for all trials
[fMRI_ROI.Ep.allTrials, choice_hE.Ep.allTrials,...
    hE_level.Ep.allTrials, Echosen.Ep.allTrials,...
    fMRI_ROI.Em.allTrials, choice_hE.Em.allTrials,...
    hE_level.Em.allTrials, Echosen.Em.allTrials] = deal(NaN(nTrialsPerTask, NS));
% data split by effort chosen
[fMRI_ROI.Ep.Ech, fMRI_ROI.Em.Ech] = deal(NaN(n_Ech,NS));
% data split by level of effort proposed for the high effort option
[fMRI_ROI.Ep.hE_level, choice_hE.Ep.hE_level,...
    fMRI_ROI.Em.hE_level, choice_hE.Em.hE_level] = deal(NaN(n_hE_levels, NS));
[fMRI_ROI.Ep.choice_low.hE_level,...
    fMRI_ROI.Ep.choice_high.hE_level,...
    fMRI_ROI.Em.choice_low.hE_level,...
    fMRI_ROI.Em.choice_high.hE_level] = deal(NaN(n_hE_levels, NS));

% slopes
[b_fMRI_f_Ech.Ep.allSubs,...
    b_fMRI_f_E.Ep.lE_chosen.allSubs,...
    b_fMRI_f_E.Ep.hE_chosen.allSubs,...
    b_fMRI_f_Ech.Em.allSubs,...
    b_fMRI_f_E.Em.lE_chosen.Em.allSubs,...
    b_fMRI_f_E.Em.hE_chosen.Em.allSubs] = deal(NaN(2,NS));

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
        switch kRun
            case {1,2}
                runTrials_idx = 1:nTrialsPerRun;
            case {3,4}
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun;
        end
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
        
        %% load the data
        behaviorStruct_tmp = load([subBehaviorFolder,...
            'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
            '_task.mat']);
        choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
        switch task_nm_tmp
            case 'Em'
                choiceAndPerf_tmp = behaviorStruct_tmp.mentalE_perf;
            case 'Ep'
                choiceAndPerf_tmp = behaviorStruct_tmp.physicalPerf;
        end
        
        %% default side
        defaultSide_tmp = choiceOptions_tmp.default_LR;
        %% extract R or P
        RP_var_tmp = strcmp(choiceOptions_tmp.R_or_P,'R');
        %% choice
        choice_LR_tmp = choiceAndPerf_tmp.choice;
        % remove confidence info from choice:
        choice_LR_tmp(choice_LR_tmp == 2) = 1;
        choice_LR_tmp(choice_LR_tmp == -2) = -1;
        % extract high effort choice
        choice_highE_tmp = NaN(1,length(choice_LR_tmp));
        choice_highE_tmp(choice_LR_tmp == -defaultSide_tmp) = 1;
        choice_highE_tmp(choice_LR_tmp == defaultSide_tmp) = 0;
        
        %% effort level
        E_highE_tmp = (choiceOptions_tmp.E.left).*(defaultSide_tmp == 1) +...
            (choiceOptions_tmp.E.right).*(defaultSide_tmp == -1);
        
        %% effort chosen
        E_chosen_tmp = (choiceOptions_tmp.E.left).*(choice_LR_tmp == -1) +...
            (choiceOptions_tmp.E.right).*(choice_LR_tmp == 1);
        
        %% extract fMRI ROI mediator
        fMRI_ROI.(task_nm_tmp).allTrials(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
        %% extract all relevant variables
        choice_hE.(task_nm_tmp).allTrials(runTrials_idx, iS)  = choice_highE_tmp;
        hE_level.(task_nm_tmp).allTrials(runTrials_idx, iS)   = E_highE_tmp;
        Echosen.(task_nm_tmp).allTrials(runTrials_idx, iS)    = E_chosen_tmp;
    end % run loop
    
    for iT = 1:nTasks
        task_nm = task_names{iT};
        %% extract data per effort level and effort chosen
        % high effort level
        for iEff_level = 1:n_hE_levels
            jEff_idx = hE_level.(task_nm).allTrials(:, iS) == iEff_level;
            % fMRI ROI depending on effort level
            fMRI_ROI.(task_nm).hE_level(iEff_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jEff_idx, iS),1,'omitnan');
            choice_hE.(task_nm).hE_level(iEff_level,iS) = mean(choice_hE.(task_nm).allTrials(jEff_idx, iS),1,'omitnan');

            % fMRI ROI depending on effort level AND choice (like in Kurniawan
            % et al, 2021)
            jEff_highEchosen_idx = (hE_level.(task_nm).allTrials(:, iS) == iEff_level).*(choice_hE.(task_nm).allTrials(:, iS) == 1) == 1;
            jEff_lowEchosen_idx = (hE_level.(task_nm).allTrials(:, iS) == iEff_level).*(choice_hE.(task_nm).allTrials(:, iS) == 0) == 1;
            fMRI_ROI.(task_nm).choice_high.hE_level(iEff_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jEff_highEchosen_idx, iS),1,'omitnan');
            fMRI_ROI.(task_nm).choice_low.hE_level(iEff_level,iS) = mean(fMRI_ROI.(task_nm).allTrials(jEff_lowEchosen_idx, iS),1,'omitnan');
        end % high effort level
        % effort chosen
        for iEff_ch = 1:n_Ech
            jEff_ch_idx = Echosen.(task_nm).allTrials(:, iS) == (iEff_ch - 1);
            fMRI_ROI.(task_nm).Ech(iEff_ch,iS) = mean(fMRI_ROI.(task_nm).allTrials(jEff_ch_idx, iS),1,'omitnan');
        end % effort chosen

        %% test whether slopes are different (similar to what the mediation tests)
        b_fMRI_f_Ech.(task_nm).allTrials.allSubs(:,iS) = glmfit(Echosen_possible, fMRI_ROI.(task_nm).Ech(:,iS), 'normal');
        b_fMRI_f_E.(task_nm).lE_chosen.allSubs(:,iS) = glmfit(hE_levels, fMRI_ROI.(task_nm).choice_low.hE_level(:,iS), 'normal');
        b_fMRI_f_E.(task_nm).hE_chosen.allSubs(:,iS) = glmfit(hE_levels, fMRI_ROI.(task_nm).choice_high.hE_level(:,iS), 'normal');
    end % task loop
end % subject loop

%% median split based on metabolites
% extract of subjects based on the median split
[low_met_subs, high_met_subs, metabolite_nm, ROI_nm] = medSplit_metabolites(study_nm, subject_id);

%% perform extraction for each task independently
for iTask = 1:nTasks
    task_nm = task_names{iTask};
    %% average data independent of metabolites
    [fMRI_ROI.(task_nm).allSubs.choice_low.hE_level.mean,...
        fMRI_ROI.(task_nm).allSubs.choice_low.hE_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hE_level,2);
    [fMRI_ROI.(task_nm).allSubs.choice_high.hE_level.mean,...
        fMRI_ROI.(task_nm).allSubs.choice_high.hE_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hE_level,2);

    % perform the median split
    % median split on choices and fMRI for high effort level
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).hE_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).hE_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).hE_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).hE_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).hE_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).hE_level(:,high_met_subs),2);
    [choice_hE.(task_nm).(['l',metabolite_nm]).hE_level.mean,...
        choice_hE.(task_nm).(['l',metabolite_nm]).hE_level.sem] = mean_sem_sd(choice_hE.(task_nm).hE_level(:,low_met_subs),2);
    [choice_hE.(task_nm).(['h',metabolite_nm]).hE_level.mean,...
        choice_hE.(task_nm).(['h',metabolite_nm]).hE_level.sem] = mean_sem_sd(choice_hE.(task_nm).hE_level(:,high_met_subs),2);
    % median split on choices and fMRI for effort chosen
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).Ech.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).Ech.sem] = mean_sem_sd(fMRI_ROI.(task_nm).Ech(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).Ech.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).Ech.sem] = mean_sem_sd(fMRI_ROI.(task_nm).Ech(:,high_met_subs),2);
    % median split on fMRI for effort level and choice made
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hE_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hE_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hE_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hE_level.mean,...
        fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hE_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hE_level(:,low_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hE_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hE_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_low.hE_level(:,high_met_subs),2);
    [fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hE_level.mean,...
        fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hE_level.sem] = mean_sem_sd(fMRI_ROI.(task_nm).choice_high.hE_level(:,high_met_subs),2);

    % betas
    % global betas with effort chosen
    b_fMRI_f_Ech.(task_nm).(['low_',metabolite_nm]) = b_fMRI_f_Ech.(task_nm).allSubs(:,low_met_subs);
    b_fMRI_f_Ech.(task_nm).(['high_',metabolite_nm]) = b_fMRI_f_Ech.(task_nm).allSubs(:,high_met_subs);
    % betas split by effort level and choice
    b_fMRI_f_E.(task_nm).lE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_E.(task_nm).lE_chosen.allSubs(:,low_met_subs);
    b_fMRI_f_E.(task_nm).lE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_E.(task_nm).lE_chosen.allSubs(:,high_met_subs);
    b_fMRI_f_E.(task_nm).hE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_E.(task_nm).hE_chosen.allSubs(:,low_met_subs);
    b_fMRI_f_E.(task_nm).hE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_E.(task_nm).hE_chosen.allSubs(:,high_met_subs);

    %% ttest
    [~,pval.(task_nm).choice_hE_vs_lE.allSubs] = ttest(b_fMRI_f_E.(task_nm).lE_chosen.allSubs(2,:),...
        b_fMRI_f_E.(task_nm).hE_chosen.allSubs(2,:));
    [~,pval.(task_nm).choice_hE_vs_choice_hE_slope.(['low_vs_high_',metabolite_nm])] = ttest2(b_fMRI_f_E.(task_nm).hE_chosen.(['low_',metabolite_nm])(2,:),...
        b_fMRI_f_E.(task_nm).hE_chosen.(['high_',metabolite_nm])(2,:));

    %% figure
    if dispFig == true
        % general parameters
        lWidth = 3;
        black = [0 0 0];
        purple = [148 0 211]./255;
        orange = [241 163 64]./255;
        blue = [51 153 255]./255;
        pSize = 50;

        %% fMRI = f(E chosen) and metabolites
        fig;
        hold on;
        fMRI_Ech_low = errorbar(Echosen_possible-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).Ech.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).Ech.sem);
        fMRI_Ech_high = errorbar(Echosen_possible+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).Ech.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).Ech.sem);
        % curve parameters
        fMRI_Ech_low.LineStyle = '--';
        fMRI_Ech_low.LineWidth = lWidth;
        fMRI_Ech_low.MarkerEdgeColor = blue;
        fMRI_Ech_high.LineStyle = '-';
        fMRI_Ech_high.LineWidth = lWidth;
        fMRI_Ech_high.MarkerEdgeColor = purple;
        legend([fMRI_Ech_high, fMRI_Ech_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(Echosen_possible);
        xlabel('Echosen');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);

        %% fMRI and choice = f(high E proposed) and metabolites
        fig;
        % fMRI
        subplot(1,2,1);
        hold on;
        fMRI_hE_low = errorbar(hE_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).hE_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).hE_level.sem);
        fMRI_hE_high = errorbar(hE_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).hE_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).hE_level.sem);
        % curve parameters
        fMRI_hE_low.LineStyle = '--';
        fMRI_hE_low.LineWidth = lWidth;
        fMRI_hE_low.MarkerEdgeColor = blue;
        fMRI_hE_high.LineStyle = '-';
        fMRI_hE_high.LineWidth = lWidth;
        fMRI_hE_high.MarkerEdgeColor = purple;
        legend([fMRI_hE_high, fMRI_hE_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xlabel('high Effort level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);

        % choice
        subplot(1,2,2);
        hold on;
        choice_hE_low = errorbar(hE_levels-0.01,...
            choice_hE.(task_nm).(['l',metabolite_nm]).hE_level.mean,...
            choice_hE.(task_nm).(['l',metabolite_nm]).hE_level.sem);
        choice_hE_high = errorbar(hE_levels+0.01,...
            choice_hE.(task_nm).(['h',metabolite_nm]).hE_level.mean,...
            choice_hE.(task_nm).(['h',metabolite_nm]).hE_level.sem);
        % curve parameters
        choice_hE_low.LineStyle = '--';
        choice_hE_low.LineWidth = lWidth;
        choice_hE_low.MarkerEdgeColor = blue;
        choice_hE_high.LineStyle = '-';
        choice_hE_high.LineWidth = lWidth;
        choice_hE_high.MarkerEdgeColor = purple;
        legend([choice_hE_high, choice_hE_low],...
            {['high ',metabolite_nm],['low ', metabolite_nm]});
        legend('Location','NorthWest');
        legend('boxoff');
        xlabel('high Effort level');
        ylabel(['Choice = high effort ',task_nm,' (%)']);
        legend_size(pSize);

        %% fMRI = f(E level) according to choice made independent of metabolites
        fig;
        hold on;
        fMRI_lowEch = errorbar(hE_levels-0.01,...
            fMRI_ROI.(task_nm).allSubs.choice_low.hE_level.mean,...
            fMRI_ROI.(task_nm).allSubs.choice_low.hE_level.sem);
        fMRI_highEch = errorbar(hE_levels+0.01,...
            fMRI_ROI.(task_nm).allSubs.choice_high.hE_level.mean,...
            fMRI_ROI.(task_nm).allSubs.choice_high.hE_level.sem);
        % curve parameters
        fMRI_lowEch.LineStyle = '--';
        fMRI_lowEch.LineWidth = lWidth;
        fMRI_lowEch.Color = blue;
        fMRI_highEch.LineStyle = '-';
        fMRI_highEch.LineWidth = lWidth;
        fMRI_highEch.Color = purple;
        legend([fMRI_highEch, fMRI_lowEch],...
            {'high E chosen','low E chosen'});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(hE_levels);
        xlabel('high Effort level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);

        %% fMRI = f(E level) according to choice made
        fig;
        hold on;
        fMRI_lowEch_lowMet = errorbar(hE_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hE_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_low.hE_level.sem);
        fMRI_lowEch_highMet = errorbar(hE_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hE_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_low.hE_level.sem);
        fMRI_highEch_lowMet = errorbar(hE_levels-0.01,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hE_level.mean,...
            fMRI_ROI.(task_nm).(['l',metabolite_nm]).choice_high.hE_level.sem);
        fMRI_highEch_highMet = errorbar(hE_levels+0.01,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hE_level.mean,...
            fMRI_ROI.(task_nm).(['h',metabolite_nm]).choice_high.hE_level.sem);
        % curve parameters
        fMRI_lowEch_lowMet.LineStyle = '--';
        fMRI_lowEch_lowMet.LineWidth = lWidth;
        fMRI_lowEch_lowMet.Color = blue;
        fMRI_lowEch_highMet.LineStyle = '-';
        fMRI_lowEch_highMet.LineWidth = lWidth;
        fMRI_lowEch_highMet.Color = blue;
        fMRI_highEch_lowMet.LineStyle = '--';
        fMRI_highEch_lowMet.LineWidth = lWidth;
        fMRI_highEch_lowMet.Color = purple;
        fMRI_highEch_highMet.LineStyle = '-';
        fMRI_highEch_highMet.LineWidth = lWidth;
        fMRI_highEch_highMet.Color = purple;
        legend([fMRI_highEch_highMet, fMRI_highEch_lowMet,...
            fMRI_lowEch_highMet, fMRI_lowEch_lowMet],...
            {['high ',metabolite_nm, '- high E chosen'],['low ', metabolite_nm, '- high E chosen'],...
            ['high ',metabolite_nm, '- low E chosen'],['low ', metabolite_nm, '- low E chosen']});
        legend('Location','NorthWest');
        legend('boxoff');
        xticks(hE_levels);
        xlabel('high Effort level');
        ylabel([ROI_short_nm,' BOLD ',task_nm]);
        legend_size(pSize);
    end % figure display
end % task loop