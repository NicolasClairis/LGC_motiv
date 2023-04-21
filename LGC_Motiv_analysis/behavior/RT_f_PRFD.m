function[stats_meanRT, stats_medianRT] = RT_f_PRFD(figDisp, n_bins)
% [stats_meanRT, stats_medianRT] = RT_f_PRFD(figDisp, n_bins)
% RT_f_PRFD tests whether there is a link between the score in the social
% dominance questionnaire (PRF-D) and reaction times across subjects
%
% INPUTS
% figDisp: display output figure (1) or not (0)
%
% n_bins: number of bins?
%
% OUTPUTS
% stats_meanRT: statistics with p.value for mean RT
%
% stats_medianRT: statistics with p.value for median RT

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% by default display figure
if ~exist('figDisp','var') || isempty(figDisp) || ~islogical(figDisp)
    figDisp = true;
end

%% define default number of bins
if ~exist('n_bins','var') || isempty(n_bins) || ~isnumeric(n_bins)
    n_bins = 6;
end

%% subject selection
[study_nm, ~, ~, subject_id, NS] = sub_id;

%% define main variables
nTrialsPerRun = 54;
tasks = {'Ep','Em'};
nTasks = length(tasks);
[PRFDscore,...
    meanRT.Ep, meanRT.Em,...
    medianRT.Ep, medianRT.Em] = deal(NaN(1,NS));
nInc = 6;
nInc_chosen = 8;
nEff = 3;
nEff_chosen = 4;
for iTask = 1:nTasks
    task_id = tasks{iTask};
    [RT_f_inc.(task_id)] = deal(NaN(nInc, NS));
    [RT_f_inc_ch.(task_id)] = deal(NaN(nInc_chosen, NS));
    [RT_f_E.(task_id)] = deal(NaN(nEff, NS));
    [RT_f_E_ch.(task_id)] = NaN(nEff_chosen, NS);
    [RT_f_R_E.(task_id),...
        RT_f_R_E_hE_chosen.(task_id),...
        RT_f_R_E_lE_chosen.(task_id)] = deal(NaN(nInc, nEff, NS));
end

%% load questionnaire results
[excelReadQuestionnairesFile] = load_questionnaires_data();
% extract CID
quest_CID_Table = excelReadQuestionnairesFile.CID;
CID_table = cell(1,length(quest_CID_Table));
for iS = 1:length(quest_CID_Table)
    if quest_CID_Table(iS) < 10
        CID_table{iS} = ['00',num2str(quest_CID_Table(iS))];
    elseif quest_CID_Table(iS) >= 10 && quest_CID_Table(iS) < 100
        CID_table{iS} = ['0',num2str(quest_CID_Table(iS))];
    elseif quest_CID_Table(iS) >= 100
        CID_table{iS} = num2str(quest_CID_Table(iS));
    end
end
PRFD_scoreTable = excelReadQuestionnairesFile.PRF_DScore;
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];
    
    %% extract PRFD score
    sub_quest_idx = strcmp(CID_table, sub_nm);
    PRFDscore(iS) = PRFD_scoreTable(sub_quest_idx);
    
    %% extract average choice RT across all tasks
    % extract information of runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    nTotalRuns = length(runsStruct.tasks);
    for iTask = 1:nTasks
        task_id = tasks{iTask};
        switch task_id
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end
        runs.(task_id) = strcmp(runsStruct.tasks,task_id);
        nRuns.(task_id) = sum(runs.(task_id));
        nSubjectTotalTrials.(task_id) = nTrialsPerRun*nRuns.(task_id);
        % prepare data to be extracted across runs
        [RT_perTrial_tmp.(task_id),...
            inc_perTrial_tmp.(task_id),...
            eff_perTrial_tmp.(task_id),...
            choice_hE_perTrial_tmp.(task_id)] = deal(NaN(nSubjectTotalTrials.(task_id), 1));
        
        jRun = 0;
        for iRun = 1:nTotalRuns
            runToInclude = 0;
            if runs.(task_id)(iRun) == 1
                jRun = jRun + 1;
                runToInclude = 1;
            end
            run_nm = num2str(runs.runsToKeep(iRun));
            
            if runToInclude == 1
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);
                
                %% load data
                behaviorStruct_tmp = load([subBehaviorFolder,...
                    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
                    '_task.mat']);
                % load relevant data
                choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
                switch task_id
                    case 'Ep'
                        onsets_tmp = behaviorStruct_tmp.physicalPerf.onsets;
                        choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                    case 'Em'
                        onsets_tmp = behaviorStruct_tmp.mentalE_perf.onsets;
                        choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
                end
                % remove confidence information from choice
                choice_LR_tmp(choice_LR_tmp == 2) = 1;
                choice_LR_tmp(choice_LR_tmp == -2) = -1;
                % extract RT
                RT_tmp = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
                % replace by NaN when no choice was performed
                RT_tmp(choice_LR_tmp == 0) = NaN;
                RT_perTrial_tmp.(task_id)(runTrials_idx) = RT_tmp;
                % extract choice made
                defaultSide_tmp = choiceOptions_tmp.default_LR;
                choice_hE_tmp = choice_LR_tmp == -defaultSide_tmp; % high effort chosen = non-default chosen
                choice_hE_tmp = double(choice_hE_tmp);
                choice_hE_tmp(choice_LR_tmp == 0) = NaN; % place NaN when no choice was made
                choice_hE_perTrial_tmp.(task_id)(runTrials_idx) = choice_hE_tmp;
                % extract level of monetary incentive proposed for high
                % effort option for each trial
                RP_var = strcmp(choiceOptions_tmp.R_or_P, 'R');
                RP_var = double(RP_var);
                RP_var(strcmp(choiceOptions_tmp.R_or_P, 'P')) = -1;
                money_level_left_tmp    = choiceOptions_tmp.monetary_amount.left.*RP_var;
                money_level_right_tmp   = choiceOptions_tmp.monetary_amount.right.*RP_var;
                inc_perTrial_tmp.(task_id)(runTrials_idx) = money_level_left_tmp.*(defaultSide_tmp == 1) +...
                    money_level_right_tmp.*(defaultSide_tmp == -1);
                % extract level of high effort option for each trial
                hE_eff_proposed_tmp = choiceOptions_tmp.E.left.*(defaultSide_tmp == 1) +...
                    choiceOptions_tmp.E.right.*(defaultSide_tmp == -1);
                eff_perTrial_tmp.(task_id)(runTrials_idx) = hE_eff_proposed_tmp;
                
            end % run to include?
        end % run loop
        
        % average the RT within each subject
        meanRT.(task_id)(iS) = mean(RT_perTrial_tmp.(task_id),1,'omitnan');
        medianRT.(task_id)(iS) = median(RT_perTrial_tmp.(task_id),1,'omitnan');
        
        % extract RT according to incentives and efforts
        for iEff = 1:nEff
            RT_f_E.(task_id)(iEff,iS) = mean(RT_perTrial_tmp.(task_id)(eff_perTrial_tmp.(task_id) == iEff),'omitnan');
        end
        for iEffch = 1:nEff_chosen
            if iEffch == 1
                RT_f_E_ch.(task_id)(iEffch,iS) = mean(RT_perTrial_tmp.(task_id)(choice_hE_perTrial_tmp.(task_id) == 0),'omitnan');
            else
                trial_idx = ((eff_perTrial_tmp.(task_id) == iEff).*(choice_hE_perTrial_tmp.(task_id) == 1)) == 1;
                RT_f_E_ch.(task_id)(iEffch,iS) = mean(RT_perTrial_tmp.(task_id)(trial_idx),'omitnan');
            end
        end
    end % physical/mental loop
end % subject loop

%% average across subjects
% mean_sem_sd()

%% test link between PRFD and RT
for iTask = 1:nTasks
    task_id = tasks{iTask};
    goodSubs = ~isnan(PRFDscore);
    [beta_meanRT.(task_id), ~,stats_meanRT.(task_id)] = glmfit(PRFDscore(goodSubs), meanRT.(task_id)(goodSubs), 'normal');
    meanRT_f_PRFD.(task_id) = glmval(beta_meanRT.(task_id), PRFDscore(goodSubs), 'identity');
    [beta_medianRT.(task_id), ~,stats_medianRT.(task_id)] = glmfit(PRFDscore(goodSubs), medianRT.(task_id)(goodSubs), 'normal');
    medianRT_f_PRFD.(task_id) = glmval(beta_medianRT.(task_id), PRFDscore(goodSubs), 'identity');
end
%% figure display
if figDisp == true
    pSize = 50;
    lWidth = 3;
    
    % look at mean and median RT = f(PRFD)
    fig;
    jPlot = 0;
    for iTask = 1:nTasks
        task_id = tasks{iTask};
        
        jPlot = jPlot + 1;
        subplot(1,4,jPlot);
        hold on;
        scatter(PRFDscore(goodSubs), meanRT.(task_id)(goodSubs));
        plot(PRFDscore(goodSubs), meanRT_f_PRFD.(task_id),...
            'LineWidth',lWidth);
        xlabel('PRF-D score');
        ylabel([task_id,' - mean RT (s)']);
        legend_size(pSize);
        
        jPlot = jPlot + 1;
        subplot(1,4,jPlot);
        hold on;
        scatter(PRFDscore(goodSubs), medianRT.(task_id)(goodSubs));
        plot(PRFDscore(goodSubs), medianRT_f_PRFD.(task_id),...
            'LineWidth',lWidth);
        xlabel('PRF-D score');
        ylabel([task_id,' - median RT (s)']);
        legend_size(pSize);
    end % task loop
end % figure display

end % function