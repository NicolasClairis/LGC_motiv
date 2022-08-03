function[] = RT_f_PRFD(figDisp, n_bins)
% RT_f_PRFD(figDisp, n_bins)
% RT_f_PRFD tests whether there is a link between the score in the social
% dominance questionnaire (PRF-D) and reaction times across subjects

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% by default display figure
if ~exist('figDisp','var') || isempty(figDisp) || ~islogical(figDisp)
    figDisp = true;
end

%% define default number of bins
if ~exist(n_bins,'var') || isempty(n_bins) || ~isnumeric(n_bins)
    n_bins = 6;
end

%% subject selection
study_nm = 'study1';
condition = subject_condition();
gender = 'all';
[subject_id, NS] = LGCM_subject_selection('study1',condition,'all');

%% define main variables
nTrialsPerRun = 54;
[PRFDscore, mRT.Ep, mRT.Em] = deal(NaN(1,NS));
nInc = 6;
nInc_chosen = 8;
nEff = 3;
nEff_chosen = 4;
[RT_f_inc.Ep, RT_f_inc.Em] = deal(NaN(nInc, NS));
[RT_f_inc_ch.Ep, RT_f_inc_ch.Em] = deal(NaN(nInc_chosen, NS));
[RT_f_E.Ep, RT_f_E.Em] = deal(NaN(nEff, NS));
[RT_f_E_ch.Ep, RT_f_E_ch.Em] = NaN(nEff_chosen, NS);
[RT_f_R_E.Ep, RT_f_R_E.Em,...
    RT_f_R_E_hE_chosen.Ep, RT_f_R_E_hE_chosen.Em,...
    RT_f_R_E_lE_chosen.Ep, RT_f_R_E_lE_chosen.Em] = deal(NaN(nInc, nEff, NS));

%% load questionnaire results
[excelReadQuestionnairesFile] = load_questionnaires_data();
% extract CID
quest_CID_Table = excelReadQuestionnairesFile.CID;
PRFD_scoreTable = excelReadQuestionnairesFile.PRF_DScore;
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];
    
    %% extract PRFD score
    sub_quest_idx = strcmp(quest_CID_Table, sub_nm);
    PRFDscore(iS) = PRFD_scoreTable(sub_quest_idx);
    
    %% extract average choice RT across all tasks
    % extract information of runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    nTotalRuns = length(runsStruct.tasks);
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
                task_fullName = 'physical';
            case 2
                task_id = 'Em';
                task_fullName = 'mental';
        end
        runs.(task_id) = strcmp(runsStruct.tasks,task_id);
        nRuns.(task_id) = sum(runs.(task_id));
        nSubjectTotalTrials.(task_id) = nTrialsPerRun*nRuns.(task_id);
        % prepare data to be extracted across runs
        [RT_perTrial_tmp.(task_id),...
            inc_perTrial_tmp.(task_id),...
            eff_perTrial_tmp.(task_id),...
            choice_hE_perTrial_tmp.(task_id)] = NaN(nSubjectTotalTrials.(task_id), 1);
        
        jRun = 0;
        for iRun = 1:nTotalRuns
            runToInclude = 0;
            if runs.(task_id)(iRun) == 1
                jRun = jRun + 1;
                runToInclude = 1;
            end
            
            if runToInclude == 1
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);
                
                %% load data
                behaviorStruct_tmp = load([subBehaviorFolder,...
                    'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                    '_task.mat']);
                % load relevant data
                onsets_tmp = behaviorStruct_tmp.onsets;
                choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
                switch task_id
                    case 'Ep'
                        onsets_tmp = behaviorStruct_tmp.physicalPerf.onsets;
                        choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                    case 'Em'
                        onsets_tmp = behaviorStruct_tmp.mentalE_perf.onsets;
                        choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
                end
                RT_tmp = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
                RT_perTrial_tmp.(task_id)(runTrials_idx) = RT_tmp;
                % extract choice made
                defaultSide_tmp = choiceOptions_tmp.default_LR;
                choice_hE_tmp = choice_LR_tmp == -defaultSide_tmp; % high effort chosen = non-default chosen
                choice_hE_tmp(choice_LR_tmp == 0) = NaN;% place NaN when no choice was made
                choice_hE_perTrial_tmp.(task_id)(runTrials_idx) = choice_hE_tmp;
                % extract level of monetary incentive proposed for high
                % effort option for each trial
                error('extraction of incentive value missing');
                % extract level of high effort option for each trial
                hE_eff_proposed_tmp = choiceOptions_tmp.E.left.*(defaultSide_tmp == 1) +...
                    choiceOptions_tmp.E.right.*(defaultSide_tmp == -1);
                eff_perTrial_tmp.(task_id)(runTrials_idx) = hE_eff_proposed_tmp;
                
            end % run to include?
        end % run loop
        
        % average the subject within each subject
        error('average RT per subject per task, per money and per effort');
    end % physical/mental loop
end % subject loop

%% average across subjects
mean_sem_sd

%% figure display
if figDisp == true

end % figure display

end % function