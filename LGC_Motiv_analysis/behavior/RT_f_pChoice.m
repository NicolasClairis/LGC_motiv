function[] = RT_f_pChoice()
% [] = RT_f_pChoice(study_nm, condition, subject_id)
% RT_f_pChoice will display the average RT (across subjects) in seconds in
% function of the probability of choosing the high effort option based on
% the bayesian model.
%
% OUTPUTS
% empty for now but ideally should also include a GLM with a test


%% subject selection
[study_nm, condition, subject_id, NS] = subject_selection;

%% working directory
computerRoot = LGCM_root_paths;
studyResultsFolder = [computerRoot, study_nm, filesep];
bayesianResultsFolder = [fullfile('C:','Users','clairis','Desktop','Github',...
    'LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'), filesep];
%% general variables
tasks = {'Ep','Em'};
nTasks = length(tasks);
nRunsPerTask = 2;
nBins = 6;
for iT = 1:nTasks
    task_nm1 = tasks{iT};
    [RT.(task_nm1).bins_perSubPerRun, pHighE.(task_nm1).bins_perSubPerRun] = deal(NaN(nBins, NS, nRunsPerTask));
    [RT.(task_nm1).bins_perSub, pHighE.(task_nm1).bins_perSub] = deal(NaN(nBins, NS));
end % task loop
mdl_nm = 'mdl_3';

%% loop across subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyResultsFolder, 'CID',sub_nm,filesep,'behavior',filesep];
    [runsStruct, n_runs] = runs_definition(study_nm, sub_nm, condition);
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    for iR = 1:n_runs
        kRun = okRuns(iR);
        run_nm = num2str(kRun);
        run_fullNm = ['run',num2str(run_nm)];
        task_nm2 = taskNames{iR};
        [task_fullName] = task_fullName_extraction(task_nm2);
        switch kRun
            case {1,2}
                jRun = 1;
            case {3,4}
                jRun = 2;
        end
        
        %% load RT
        RT_allTrials_tmp = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        
        %% load p(choice = high E) based on model
        [~, ~, ~, pChoice_allTrials_tmp] = extract_bayesian_mdl(bayesianResultsFolder, subBehaviorFolder,...
            sub_nm, run_nm, task_fullName, mdl_nm);
        
        %% extract bins for the current run
        [RT.(task_nm2).bins_perSubPerRun(:,iS,jRun),...
            pHighE.(task_nm2).bins_perSubPerRun(:,iS,jRun)] = do_bin2(RT_allTrials_tmp,...
            pChoice_allTrials_tmp, nBins, 0);
    end % run loop
    
    %% average runs within each subject
    for iT = 1:nTasks
        task_nm3 = tasks{iT};
        RT.(task_nm3).bins_perSub(:,iS) = mean(RT.(task_nm3).bins_perSubPerRun(:,iS,:),3,'omitnan');
        pHighE.(task_nm3).bins_perSub(:,iS) = mean(pHighE.(task_nm3).bins_perSubPerRun(:,iS,:),3,'omitnan');
    end % task loop
end % subject loop

%% average across subjects + figure
fig;
[pSize, lWidth, col, mSize] = general_fig_prm;
for iT = 1:nTasks
    task_nm4 = tasks{iT};
    % average + SEM
    [RT.(task_nm4).bins_m_aSubs,...
        RT.(task_nm4).bins_sem_aSubs] = mean_sem_sd(RT.(task_nm4).bins_perSub, 2);
    [pHighE.(task_nm4).bins_m_aSubs,...
        pHighE.(task_nm4).bins_sem_aSubs] = mean_sem_sd(pHighE.(task_nm4).bins_perSub, 2);
    
    % figure
    subplot(1,nTasks,iT); hold on;
    er_hdl = errorbar(pHighE.(task_nm4).bins_m_aSubs.*100,...
        RT.(task_nm4).bins_m_aSubs,...
        RT.(task_nm4).bins_sem_aSubs);
    errorbar_hdl_upgrade(er_hdl);
    xlim([0 100]);
    ylim([1.8 2.5]);
    xlabel('p(high E) (%)');
    ylabel('RT (s)');
    legend_size(pSize);
end % task loop
end % function