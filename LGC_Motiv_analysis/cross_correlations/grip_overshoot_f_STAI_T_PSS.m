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
[excelReadQuestionnairesFile, sub_questionnaire_CID_list] = load_questionnaires_data;

%% extract grip overshoot and STAI-T
[STAI_T, PSS14,...
    grip_avg_AUC_overshoot.allTrials,...
    grip_avg_peakF.allTrials] = deal(NaN(1,NS));
n_Ech = 4;
[grip_avg_AUC_overshoot.per_Ech,...
    grip_avg_peakF.per_Ech] = deal(NaN(n_Ech, NS));
nGripRuns = 2;
nTrialsPerRun = 54;
nGripTrials = nTrialsPerRun*nGripRuns;
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    %% extract STAI-T data
    questionnaire_idx = strcmp(sub_questionnaire_CID_list, sub_nm);
    STAI_T(iS) = excelReadQuestionnairesFile.STAITraitScore(questionnaire_idx);
    PSS14(iS) = excelReadQuestionnairesFile.PSS_14Score(questionnaire_idx);
    %% extract grip overshoot
    [AUC_overshoot_perSub.allTrials,...
        peakF_perSub.allTrials] = deal(NaN(nGripTrials, 1));
    [AUC_overshoot_perSub.per_Ech,...
        peakF_perSub.per_Ech]  = deal(NaN(n_Ech, nGripRuns));
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
                [~, ~, forcePeak_tmp, AUC_overshoot_tmp] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
                AUC_overshoot_perSub.allTrials(runTrials_idx) = AUC_overshoot_tmp.allTrials;
                AUC_overshoot_perSub.per_Ech(:,jGripRun) = AUC_overshoot_tmp.per_Ech;
                peakF_perSub.allTrials(runTrials_idx) = forcePeak_tmp.allTrials;
                peakF_perSub.per_Ech(:,jGripRun) = forcePeak_tmp.per_Ech;
        end
    end % run loop

    %% average data per subject
    grip_avg_AUC_overshoot.allTrials(iS) = mean(AUC_overshoot_perSub.allTrials,1,'omitnan');
    grip_avg_AUC_overshoot.per_Ech(:,iS) = mean(AUC_overshoot_perSub.per_Ech,2,'omitnan');
    grip_avg_peakF.allTrials(iS) = mean(peakF_perSub.allTrials,1,'omitnan');
    grip_avg_peakF.per_Ech(:,iS) = mean(peakF_perSub.per_Ech,2,'omitnan');
end % subject loop

%% test correlation
% STAI-T
[b_STAI_avg_AUC_overshoot, ~,stats_STAI_avg_overshoot] = glmfit(STAI_T, grip_avg_AUC_overshoot.allTrials, 'normal');
STAI_T_fit = sort(STAI_T);
avg_AUC_overshoot_STAI_T_fit = glmval(b_STAI_avg_AUC_overshoot, STAI_T_fit, 'identity');
pval.STAI.AUC_overshoot = stats_STAI_avg_overshoot.p;

[b_STAI_avg_peakF, ~,stats_STAI_avg_peakF] = glmfit(STAI_T, grip_avg_peakF.allTrials, 'normal');
avg_peakF_STAI_T_fit = glmval(b_STAI_avg_peakF, STAI_T_fit, 'identity');
pval.STAI.peakF = stats_STAI_avg_peakF.p;

% PSS-14
[b_PSS14_avg_AUC_overshoot, ~,stats_PSS14_avg_overshoot] = glmfit(PSS14, grip_avg_AUC_overshoot.allTrials, 'normal');
PSS14_fit = sort(PSS14);
avg_AUC_overshoot_PSS14_fit = glmval(b_PSS14_avg_AUC_overshoot, PSS14_fit, 'identity');
pval.PSS14.AUC_overshoot = stats_PSS14_avg_overshoot.p;

[b_PSS14_avg_peakF, ~,stats_PSS14_avg_peakF] = glmfit(PSS14, grip_avg_peakF.allTrials, 'normal');
avg_peakF_PSS14_fit = glmval(b_PSS14_avg_peakF, PSS14_fit, 'identity');
pval.PSS14.peakF = stats_PSS14_avg_peakF.p;

%% figure display
pSize = 30;
lWidth = 3;
black = [0 0 0];

%% STAI-T
% AUC overshoot
fig;
hold on;
scat_hdl = scatter(STAI_T, grip_avg_AUC_overshoot.allTrials);
fit_hdl = plot(STAI_T_fit, avg_AUC_overshoot_STAI_T_fit);
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

% peak force
fig;
hold on;
scat_hdl = scatter(STAI_T, grip_avg_peakF.allTrials);
fit_hdl = plot(STAI_T_fit, avg_peakF_STAI_T_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '--';
xlabel('STAI-T');
ylabel('peak force (%)');
legend_size(pSize);

%% PSS-14
% AUC overshoot
fig;
hold on;
scat_hdl = scatter(PSS14, grip_avg_AUC_overshoot.allTrials);
fit_hdl = plot(PSS14_fit, avg_AUC_overshoot_PSS14_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '--';
xlabel('PSS14');
ylabel('average AUC overshoot');
legend_size(pSize);

% peak force
fig;
hold on;
scat_hdl = scatter(PSS14, grip_avg_peakF.allTrials);
fit_hdl = plot(PSS14_fit, avg_peakF_PSS14_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '--';
xlabel('PSS14');
ylabel('peak force (%)');
legend_size(pSize);