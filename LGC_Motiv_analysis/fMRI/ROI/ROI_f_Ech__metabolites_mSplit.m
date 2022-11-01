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
which_task = listdlg('PromptString','Which task?','ListString',task_names);
task_to_look = task_names{which_task};
% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};

%% load behavior 
nRuns = 4;
nTrialsPerRun = 54;
nTrials = nTrialsPerRun*nRuns;
% main variables of interest
Echosen_possible = 0:3;
n_Ech = length(Echosen_possible);
hE_levels = 1:3;
n_hE_levels = length(hE_levels);
% data for all trials
[fMRI_ROI.allTrials, choice_hE.allTrials,...
    hE_level.allTrials, Echosen.allTrials] = deal(NaN(nTrials, NS));
% data split by effort chosen
fMRI_ROI.Ech = NaN(n_Ech,NS);
% data split by level of effort proposed for the high effort option
[fMRI_ROI.hE_level, choice_hE.hE_level] = deal(NaN(n_hE_levels, NS));

% slopes
b_fMRI_f_Ech.allSubs = NaN(2,NS);

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
        
        %% extract fMRI ROI mediator
        if strcmp(task_to_look,'EpEmPool') ||...
                (strcmp(task_to_look, task_nm))
            fMRI_ROI.allTrials(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
        end
        %% extract all relevant variables
        choice_hE.allTrials(runTrials_idx, iS)  = choice_highE_tmp;
        hE_level.allTrials(runTrials_idx, iS)   = E_highE_tmp;
        Echosen.allTrials(runTrials_idx, iS)    = E_chosen_tmp;
    end % run loop
    
    %% extract data per effort level and effort chosen
    % high effort level
    for iEff_level = 1:n_hE_levels
        jEff_idx = hE_level.allTrials(:, iS) == iEff_level;
        fMRI_ROI.hE_level(iEff_level,iS) = mean(fMRI_ROI.allTrials(jEff_idx, iS),1,'omitnan');
        choice_hE.hE_level(iEff_level,iS) = mean(choice_hE.allTrials(jEff_idx, iS),1,'omitnan');
    end % high effort level
    % effort chosen
    for iEff_ch = 1:n_Ech
        jEff_ch_idx = Echosen.allTrials(:, iS) == (iEff_ch - 1);
        fMRI_ROI.Ech(iEff_ch,iS) = mean(fMRI_ROI.allTrials(jEff_ch_idx, iS),1,'omitnan');
    end % effort chosen
    
    %% test whether slopes are different (similar to what the mediation tests)
    b_fMRI_f_Ech.allSubs(:,iS) = glmfit(Echosen_possible, fMRI_ROI.Ech(:,iS), 'normal');
end % subject loop

%% median split based on metabolites
[low_met_subs, high_met_subs, metabolite_nm, ROI_nm] = medSplit_metabolites(study_nm, subject_id);
[fMRI_ROI.(['l',metabolite_nm]).hE_level.mean,...
    fMRI_ROI.(['l',metabolite_nm]).hE_level.sem] = mean_sem_sd(fMRI_ROI.hE_level(:,low_met_subs),2);
[fMRI_ROI.(['h',metabolite_nm]).hE_level.mean,...
    fMRI_ROI.(['h',metabolite_nm]).hE_level.sem] = mean_sem_sd(fMRI_ROI.hE_level(:,high_met_subs),2);
[choice_hE.(['l',metabolite_nm]).hE_level.mean,...
    choice_hE.(['l',metabolite_nm]).hE_level.sem] = mean_sem_sd(choice_hE.hE_level(:,low_met_subs),2);
[choice_hE.(['h',metabolite_nm]).hE_level.mean,...
    choice_hE.(['h',metabolite_nm]).hE_level.sem] = mean_sem_sd(choice_hE.hE_level(:,high_met_subs),2);
[fMRI_ROI.(['l',metabolite_nm]).Ech.mean,...
    fMRI_ROI.(['l',metabolite_nm]).Ech.sem] = mean_sem_sd(fMRI_ROI.Ech(:,low_met_subs),2);
[fMRI_ROI.(['h',metabolite_nm]).Ech.mean,...
    fMRI_ROI.(['h',metabolite_nm]).Ech.sem] = mean_sem_sd(fMRI_ROI.Ech(:,high_met_subs),2);
% betas
b_fMRI_f_Ech.(['low_',metabolite_nm]) = b_fMRI_f_Ech.allSubs(:,low_met_subs);
b_fMRI_f_Ech.(['high_',metabolite_nm]) = b_fMRI_f_Ech.allSubs(:,high_met_subs);


%% figure
if dispFig == true
    % general parameters
    lWidth = 3;
    black = [0 0 0];
    purple = [153 142 195]./255;
    orange = [241 163 64]./255;
    blue = [0 255 0]./255;
    pSize = 50;
    
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
    fMRI_Ech_high.LineStyle = '--';
    fMRI_Ech_high.LineWidth = lWidth;
    fMRI_Ech_high.MarkerEdgeColor = purple;
    legend([fMRI_Ech_high, fMRI_Ech_low],...
        {['high ',metabolite_nm],['low ', metabolite_nm]});
    legend('Location','NorthWest');
    legend('boxoff');
    xlabel('Echosen');
    ylabel([ROI_short_nm,' BOLD']);
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
    fMRI_hE_high.LineStyle = '--';
    fMRI_hE_high.LineWidth = lWidth;
    fMRI_hE_high.MarkerEdgeColor = purple;
    legend([fMRI_hE_high, fMRI_hE_low],...
        {['high ',metabolite_nm],['low ', metabolite_nm]});
    legend('Location','NorthWest');
    legend('boxoff');
    xlabel('high Effort level');
    ylabel([ROI_short_nm,' BOLD']);
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
    choice_hE_high.LineStyle = '--';
    choice_hE_high.LineWidth = lWidth;
    choice_hE_high.MarkerEdgeColor = purple;
    legend([choice_hE_high, choice_hE_low],...
        {['high ',metabolite_nm],['low ', metabolite_nm]});
    legend('Location','NorthWest');
    legend('boxoff');
    xlabel('high Effort level');
    ylabel('Choice = high effort (%)');
    legend_size(pSize);
end % figure display