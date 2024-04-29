%% extract the level of activity of a given ROI in function of each effort level chosen
% and split subjects with a median split based on the metabolites. the same
% script will also look at the choice depending on dmPFC activity and
% effort level proposed for the high effort option

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
task_names = {'Ep','Em','EpEmPool'};
which_task = listdlg('PromptString','Which task for ROI?','ListString',task_names);
ROI_task_to_look = task_names{which_task};
% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};

%% load behavior
which_bhv_task = listdlg('PromptString','Which task for behavior?','ListString',task_names);
behavioral_task_to_look = task_names{which_bhv_task};
nTrialsPerRun = 54;
switch behavioral_task_to_look
    case {'Ep','Em'}
        nRunsPerTask = 2;
        nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
        % data for all trials
        [fMRI_ROI.allTrials, choice_hE.allTrials,...
            hE_level.allTrials, Echosen.allTrials] = deal(NaN(nTrialsPerTask, NS));
    case 'EpEmPool'
        nRuns = 4;
        nTrials = nTrialsPerRun*nRuns;
        % data for all trials
        [fMRI_ROI.allTrials, choice_hE.allTrials,...
            hE_level.allTrials, Echosen.allTrials] = deal(NaN(nTrials, NS));
end
% main variables of interest
Echosen_possible = 0:3;
n_Ech = length(Echosen_possible);
hE_levels = 1:3;
n_hE_levels = length(hE_levels);
% data split by effort chosen
fMRI_ROI.Ech = NaN(n_Ech,NS);
% data split by level of effort proposed for the high effort option
[fMRI_ROI.hE_level, choice_hE.hE_level] = deal(NaN(n_hE_levels, NS));
[fMRI_ROI.choice_low.hE_level,...
    fMRI_ROI.choice_high.hE_level] = deal(NaN(n_hE_levels, NS));

% slopes
[b_fMRI_f_Ech.allSubs,...
    b_fMRI_f_E.lE_chosen.allSubs,...
    b_fMRI_f_E.hE_chosen.allSubs] = deal(NaN(2,NS));

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
        % define trial index for relevant variable to extract
        switch behavioral_task_to_look
            case {'Ep','Em'}
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(taskRun_idx-1);
            case 'EpEmPool'
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(kRun-1);
        end
        
        %% load the data
        % filter task based on what was selected in the inputs
        if strcmp(behavioral_task_to_look,'EpEmPool') ||...
                strcmp(behavioral_task_to_look, task_nm_tmp)
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
            
            %% high-effort money amount
            money_hE_tmp = ((choiceOptions_tmp.monetary_amount.left).*(defaultSide_tmp == 1) +...
                (choiceOptions_tmp.monetary_amount.right).*(defaultSide_tmp == -1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
            money_lE_tmp = ((choiceOptions_tmp.monetary_amount.left).*(defaultSide_tmp == -1) +...
                (choiceOptions_tmp.monetary_amount.right).*(defaultSide_tmp == 1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
            
            %% delta between high and low effort options
            deltaMoney_tmp = money_hE_tmp - money_lE_tmp;
            
            %% RT
            onsets_tmp = behaviorStruct_tmp.onsets;
            switch task_nm_tmp
                case 'Ep'
                    onsets_tmp = behaviorStruct_tmp.physicalPerf.onsets;
                    choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                case 'Em'
                    onsets_tmp = behaviorStruct_tmp.mentalE_perf.onsets;
                    choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
            end
            RT_tmp = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
            
            %% confidence rating
            uncertaintyRtg_tmp = NaN(1,length(choice_LR_tmp));
            uncertaintyRtg_tmp(abs(choice_LR_tmp) == 2) = 0; % high confidence = low uncertainty
            uncertaintyRtg_tmp(abs(choice_LR_tmp) == 1) = 1; % low confidence = high uncertainty
            
            %% extract all relevant variables
            choice_hE.allTrials(runTrials_idx, iS)  = choice_highE_tmp;
            hE_level.allTrials(runTrials_idx, iS)   = E_highE_tmp;
            Echosen.allTrials(runTrials_idx, iS)    = E_chosen_tmp;
        end % behavioral task filter
        %% extract fMRI ROI mediator
        if strcmp(ROI_task_to_look,'EpEmPool') ||...
                (strcmp(ROI_task_to_look, task_nm_tmp))
            fMRI_ROI.allTrials(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
        end
        
    end % run loop
    
    %% extract data per effort level and effort chosen
    % high effort level
    for iEff_level = 1:n_hE_levels
        jEff_idx = hE_level.allTrials(:, iS) == iEff_level;
        % fMRI ROI depending on effort level
        fMRI_ROI.hE_level(iEff_level,iS) = mean(fMRI_ROI.allTrials(jEff_idx, iS),1,'omitnan');
        choice_hE.hE_level(iEff_level,iS) = mean(choice_hE.allTrials(jEff_idx, iS),1,'omitnan');

        % fMRI ROI depending on effort level AND choice (like in Kurniawan
        % et al, 2021)
        jEff_highEchosen_idx = (hE_level.allTrials(:, iS) == iEff_level).*(choice_hE.allTrials(:, iS) == 1) == 1;
        jEff_lowEchosen_idx = (hE_level.allTrials(:, iS) == iEff_level).*(choice_hE.allTrials(:, iS) == 0) == 1;
        fMRI_ROI.choice_high.hE_level(iEff_level,iS) = mean(fMRI_ROI.allTrials(jEff_highEchosen_idx, iS),1,'omitnan');
        fMRI_ROI.choice_low.hE_level(iEff_level,iS) = mean(fMRI_ROI.allTrials(jEff_lowEchosen_idx, iS),1,'omitnan');
    end % high effort level
    % effort chosen
    for iEff_ch = 1:n_Ech
        jEff_ch_idx = Echosen.allTrials(:, iS) == (iEff_ch - 1);
        fMRI_ROI.Ech(iEff_ch,iS) = mean(fMRI_ROI.allTrials(jEff_ch_idx, iS),1,'omitnan');
    end % effort chosen
    
    %% test whether slopes are different (similar to what the mediation tests)
    b_fMRI_f_Ech.allTrials.allSubs(:,iS) = glmfit(Echosen_possible, fMRI_ROI.Ech(:,iS), 'normal');
    b_fMRI_f_E.lE_chosen.allSubs(:,iS) = glmfit(hE_levels, fMRI_ROI.choice_low.hE_level(:,iS), 'normal');
    b_fMRI_f_E.hE_chosen.allSubs(:,iS) = glmfit(hE_levels, fMRI_ROI.choice_high.hE_level(:,iS), 'normal');
end % subject loop

%% average data independent of metabolites
[fMRI_ROI.allSubs.choice_low.hE_level.mean,...
    fMRI_ROI.allSubs.choice_low.hE_level.sem] = mean_sem_sd(fMRI_ROI.choice_low.hE_level,2);
[fMRI_ROI.allSubs.choice_high.hE_level.mean,...
    fMRI_ROI.allSubs.choice_high.hE_level.sem] = mean_sem_sd(fMRI_ROI.choice_high.hE_level,2);

%% median split based on metabolites
% extract of subjects based on the median split
[low_met_subs, high_met_subs, metabolite_nm, MRS_ROI_nm] = medSplit_metabolites(study_nm, subject_id);

% perform the median split
% median split on choices and fMRI for high effort level
[fMRI_ROI.(['l',metabolite_nm]).hE_level.mean,...
    fMRI_ROI.(['l',metabolite_nm]).hE_level.sem] = mean_sem_sd(fMRI_ROI.hE_level(:,low_met_subs),2);
[fMRI_ROI.(['h',metabolite_nm]).hE_level.mean,...
    fMRI_ROI.(['h',metabolite_nm]).hE_level.sem] = mean_sem_sd(fMRI_ROI.hE_level(:,high_met_subs),2);
[choice_hE.(['l',metabolite_nm]).hE_level.mean,...
    choice_hE.(['l',metabolite_nm]).hE_level.sem] = mean_sem_sd(choice_hE.hE_level(:,low_met_subs),2);
[choice_hE.(['h',metabolite_nm]).hE_level.mean,...
    choice_hE.(['h',metabolite_nm]).hE_level.sem] = mean_sem_sd(choice_hE.hE_level(:,high_met_subs),2);
% median split on choices and fMRI for effort chosen
[fMRI_ROI.(['l',metabolite_nm]).Ech.mean,...
    fMRI_ROI.(['l',metabolite_nm]).Ech.sem] = mean_sem_sd(fMRI_ROI.Ech(:,low_met_subs),2);
[fMRI_ROI.(['h',metabolite_nm]).Ech.mean,...
    fMRI_ROI.(['h',metabolite_nm]).Ech.sem] = mean_sem_sd(fMRI_ROI.Ech(:,high_met_subs),2);
% median split on fMRI for effort level and choice made
[fMRI_ROI.(['l',metabolite_nm]).choice_low.hE_level.mean,...
    fMRI_ROI.(['l',metabolite_nm]).choice_low.hE_level.sem] = mean_sem_sd(fMRI_ROI.choice_low.hE_level(:,low_met_subs),2);
[fMRI_ROI.(['l',metabolite_nm]).choice_high.hE_level.mean,...
    fMRI_ROI.(['l',metabolite_nm]).choice_high.hE_level.sem] = mean_sem_sd(fMRI_ROI.choice_high.hE_level(:,low_met_subs),2);
[fMRI_ROI.(['h',metabolite_nm]).choice_low.hE_level.mean,...
    fMRI_ROI.(['h',metabolite_nm]).choice_low.hE_level.sem] = mean_sem_sd(fMRI_ROI.choice_low.hE_level(:,high_met_subs),2);
[fMRI_ROI.(['h',metabolite_nm]).choice_high.hE_level.mean,...
    fMRI_ROI.(['h',metabolite_nm]).choice_high.hE_level.sem] = mean_sem_sd(fMRI_ROI.choice_high.hE_level(:,high_met_subs),2);

% betas
% global betas with effort chosen
b_fMRI_f_Ech.(['low_',metabolite_nm]) = b_fMRI_f_Ech.allSubs(:,low_met_subs);
b_fMRI_f_Ech.(['high_',metabolite_nm]) = b_fMRI_f_Ech.allSubs(:,high_met_subs);
% betas split by effort level and choice
b_fMRI_f_E.lE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_E.lE_chosen.allSubs(:,low_met_subs);
b_fMRI_f_E.lE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_E.lE_chosen.allSubs(:,high_met_subs);
b_fMRI_f_E.hE_chosen.(['low_',metabolite_nm]) = b_fMRI_f_E.hE_chosen.allSubs(:,low_met_subs);
b_fMRI_f_E.hE_chosen.(['high_',metabolite_nm]) = b_fMRI_f_E.hE_chosen.allSubs(:,high_met_subs);

%% ttest
[~,pval.choice_hE_vs_lE.allSubs] = ttest(b_fMRI_f_E.lE_chosen.allSubs(2,:), b_fMRI_f_E.hE_chosen.allSubs(2,:));
[~,pval.choice_hE_vs_choice_hE_slope.(['low_vs_high_',metabolite_nm])] = ttest2(b_fMRI_f_E.hE_chosen.(['low_',metabolite_nm])(2,:),...
    b_fMRI_f_E.hE_chosen.(['high_',metabolite_nm])(2,:));

%% figure
if dispFig == true
    % general parameters
    lWidth = 3;
    black = [0 0 0];
    purple = [148 0 211]./255;
    orange = [241 163 64]./255;
    blue = [51 153 255]./255;
    pSize = 50;
    switch behavioral_task_to_look
        case 'Ep'
            full_bhv_taskName = 'physical task';
        case 'Em'
            full_bhv_taskName = 'mental task';
        case 'EpEmPool'
            full_bhv_taskName = 'both tasks';
    end
    switch ROI_task_to_look
        case 'Ep'
            full_ROI_taskName = 'physical task';
        case 'Em'
            full_ROI_taskName = 'mental task';
        case 'EpEmPool'
            full_ROI_taskName = 'both tasks';
    end
    
    %% fMRI = f(E chosen) and metabolites
    fig;
    hold on;
    fMRI_Ech_low = errorbar(Echosen_possible-0.01,...
        fMRI_ROI.(['l',metabolite_nm]).Ech.mean,...
        fMRI_ROI.(['l',metabolite_nm]).Ech.sem);
    fMRI_Ech_high = errorbar(Echosen_possible+0.01,...
        fMRI_ROI.(['h',metabolite_nm]).Ech.mean,...
        fMRI_ROI.(['h',metabolite_nm]).Ech.sem);
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
    xlabel({'Echosen';full_bhv_taskName});
    ylabel({[ROI_short_nm,' BOLD during ',timePeriod_nm];...
        full_ROI_taskName});
    legend_size(pSize);
    
    %% fMRI and choice = f(high E proposed) and metabolites
    fig;
    % fMRI
    subplot(1,2,1);
    hold on;
    fMRI_hE_low = errorbar(hE_levels-0.01,...
        fMRI_ROI.(['l',metabolite_nm]).hE_level.mean,...
        fMRI_ROI.(['l',metabolite_nm]).hE_level.sem);
    fMRI_hE_high = errorbar(hE_levels+0.01,...
        fMRI_ROI.(['h',metabolite_nm]).hE_level.mean,...
        fMRI_ROI.(['h',metabolite_nm]).hE_level.sem);
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
    xlabel({'high Effort level';...
        full_bhv_taskName});
    ylabel({[ROI_short_nm,' BOLD during ',timePeriod_nm],...
        full_ROI_taskName});
    legend_size(pSize);
    
    % choice
    subplot(1,2,2);
    hold on;
    choice_hE_low = errorbar(hE_levels-0.01,...
        choice_hE.(['l',metabolite_nm]).hE_level.mean,...
        choice_hE.(['l',metabolite_nm]).hE_level.sem);
    choice_hE_high = errorbar(hE_levels+0.01,...
        choice_hE.(['h',metabolite_nm]).hE_level.mean,...
        choice_hE.(['h',metabolite_nm]).hE_level.sem);
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
    xlabel({'high Effort level';...
        full_bhv_taskName});
    ylabel({'Choice = high effort (%)';...
        full_bhv_taskName});
    legend_size(pSize);

    %% fMRI = f(E level) according to choice made independent of metabolites
    fig;
    hold on;
    fMRI_lowEch = errorbar(hE_levels-0.01,...
        fMRI_ROI.allSubs.choice_low.hE_level.mean,...
        fMRI_ROI.allSubs.choice_low.hE_level.sem);
    fMRI_highEch = errorbar(hE_levels+0.01,...
        fMRI_ROI.allSubs.choice_high.hE_level.mean,...
        fMRI_ROI.allSubs.choice_high.hE_level.sem);
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
    xlabel({'high Effort level';...
        full_bhv_taskName});
    ylabel({[ROI_short_nm,' BOLD during ',timePeriod_nm];...
        full_ROI_taskName});
    legend_size(pSize);

    %% fMRI = f(E level) according to choice made
    fig;
    hold on;
    fMRI_lowEch_lowMet = errorbar(hE_levels-0.01,...
        fMRI_ROI.(['l',metabolite_nm]).choice_low.hE_level.mean,...
        fMRI_ROI.(['l',metabolite_nm]).choice_low.hE_level.sem);
    fMRI_lowEch_highMet = errorbar(hE_levels+0.01,...
        fMRI_ROI.(['h',metabolite_nm]).choice_low.hE_level.mean,...
        fMRI_ROI.(['h',metabolite_nm]).choice_low.hE_level.sem);
    fMRI_highEch_lowMet = errorbar(hE_levels-0.01,...
        fMRI_ROI.(['l',metabolite_nm]).choice_high.hE_level.mean,...
        fMRI_ROI.(['l',metabolite_nm]).choice_high.hE_level.sem);
    fMRI_highEch_highMet = errorbar(hE_levels+0.01,...
        fMRI_ROI.(['h',metabolite_nm]).choice_high.hE_level.mean,...
        fMRI_ROI.(['h',metabolite_nm]).choice_high.hE_level.sem);
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
    xlabel({'high Effort';full_bhv_taskName});
    ylabel({[ROI_short_nm,' BOLD during ',timePeriod_nm];...
        full_ROI_taskName});
    legend_size(pSize);
end % figure display