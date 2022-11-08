%% script to check whether people with a higher trait anxiety assessed 
% through the STAI-T questionnaire also overshoot more in the physical
% effort task.

%% subject definition
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% load STAI-T questionnaire data
[excelReadQuestionnairesFile, sub_STAI_CID_list] = load_questionnaires_data;

%% extract grip overshoot and STAI-T
[STAI_T,...
    grip_avg_AUC_overshoot.allTrials] = deal(NaN(1,NS));
n_Ech = 4;
grip_avg_AUC_overshoot.per_Ech = NaN(n_Ech, NS);
nGripRuns = 2;
nTrialsPerRun = 54;
nGripTrials = nTrialsPerRun*nGripRuns;
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    %% extract STAI-T data
    STAI_ixd = strcmp(sub_STAI_CID_list, sub_nm);
    STAI_T(iS) = excelReadQuestionnairesFile.STAITraitScore(STAI_ixd);

    %% extract grip overshoot
    AUC_overshoot_perSub.allTrials = NaN(nGripTrials, 1);
    AUC_overshoot_perSub.per_Ech = NaN(n_Ech, nGripRuns);
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
                jGripRun = 1;
            case {3,4}
                jGripRun = 2;
        end
        run_nm_bis = ['run',num2str(jGripRun)];
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jGripRun-1);
        switch task_nm_tmp
            case 'Ep'
                [~, ~, ~, AUC_overshoot_tmp] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
                AUC_overshoot_perSub.allTrials(runTrials_idx) = AUC_overshoot_tmp.allTrials;
                AUC_overshoot_perSub.per_Ech(:,jGripRun) = AUC_overshoot_tmp.per_Ech;
        end
    end % run loop

    %% average data per subject
    grip_avg_AUC_overshoot.allTrials(iS) = mean(AUC_overshoot_perSub.allTrials,1,'omitnan');
    grip_avg_AUC_overshoot.per_Ech(:,iS) = mean(AUC_overshoot_perSub.per_Ech,2,'omitnan');
end % subject loop

%% test correlation
[b_STAI_avg_AUC_overshoot, ~,stats_STAI_avg_overshoot] = glmfit(STAI_T, grip_avg_AUC_overshoot.allTrials, 'normal');
STAI_T_fit = sort(STAI_T);
avg_AUC_overshoot_fit = glmval(b_STAI_avg_AUC_overshoot, STAI_T_fit, 'identity');

%% figure display
pSize = 30;
lWidth = 3;
black = [0 0 0];

fig;
hold on;
scat_hdl = scatter(STAI_T, grip_avg_AUC_overshoot.allTrials);
fit_hdl = plot(STAI_T_fit, avg_AUC_overshoot_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '--';
xlabel('STAI-T');
ylabel('average AUC overshoot');
legend_size(pSize);

% fig;
% scatter(STAI_T, grip_avg_AUC_overshoot.allTrials);
% xlabel('STAI-T');
% ylabel('average AUC overshoot');
% legend_size(pSize);