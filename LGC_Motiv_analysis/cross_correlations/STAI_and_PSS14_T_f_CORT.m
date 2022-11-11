%% compare MADRS-S and interleukin levels

%% load subject
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load questionnaires
[excelReadQuestionnairesFile, quest_sub_CID_list] = load_questionnaires_data;

%% load cortisol
[CORT_data] = load_CORT;

%% extract the data for the selected subjects
[CORT_AUCg, STAI_T_score, PSS14_score] = deal(NaN(1, NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    CORT_sub_idx = strcmp(['CID',sub_nm], CORT_data.CID);
    quest_sub_idx = strcmp(sub_nm, quest_sub_CID_list);

    CORT_AUCg(iS) = CORT_data.AUCg(CORT_sub_idx);
    STAI_T_score(iS) = excelReadQuestionnairesFile.STAITraitScore(quest_sub_idx);
    PSS14_score(iS) = excelReadQuestionnairesFile.PSS_14Score(quest_sub_idx);
end % subject loop

%% do they correlate with each other
% STAI-T vs CORT
CORT_STAI_goodSubs = ~isnan(CORT_AUCg).*~isnan(STAI_T_score) == 1;
[b_CORT_STAI_T,~,stats_tmp] = glmfit(CORT_AUCg(CORT_STAI_goodSubs),...
    STAI_T_score(CORT_STAI_goodSubs), 'normal');
pval.CORT_STAI_T = stats_tmp.p;
CORT_for_fit = sort(CORT_AUCg(CORT_STAI_goodSubs));
CORT_STAI_T_fit = glmval(b_CORT_STAI_T, CORT_for_fit, 'identity');

% PSS14 vs CORT
CORT_PSS14_goodSubs = ~isnan(CORT_AUCg).*~isnan(PSS14_score) == 1;
[b_CORT_PSS14,~,stats_tmp] = glmfit(CORT_AUCg(CORT_PSS14_goodSubs),...
    PSS14_score(CORT_PSS14_goodSubs), 'normal');
pval.CORT_PSS14 = stats_tmp.p;
CORT_for_fit_bis = sort(CORT_AUCg(CORT_PSS14_goodSubs));
CORT_PSS14_fit = glmval(b_CORT_PSS14, CORT_for_fit_bis, 'identity');

%% figures
pSize = 50;
lWidth = 3;
black = [0 0 0];
grey = [143 143 143]./255;
%% test linear regression
fig;

% STAI-T = f(CORT)
subplot(1,2,1);
hold on;
scat_hdl = scatter(CORT_AUCg(CORT_STAI_goodSubs),...
    STAI_T_score(CORT_STAI_goodSubs));
fit_hdl = plot(CORT_for_fit, CORT_STAI_T_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
xlabel('Cortisol AUCg');
ylabel('STAI-T');
legend_size(pSize);

% PSS14 = f(CORT)
subplot(1,2,2);
hold on;
scat_hdl = scatter(CORT_AUCg(CORT_PSS14_goodSubs),...
    PSS14_score(CORT_PSS14_goodSubs));
fit_hdl = plot(CORT_for_fit_bis, CORT_PSS14_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
xlabel('Cortisol AUCg');
ylabel('PSS-14');
legend_size(pSize);
