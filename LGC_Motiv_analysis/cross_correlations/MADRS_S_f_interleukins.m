%% compare MADRS-S and interleukin levels

%% load subject
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load MADRS-S
[excelReadQuestionnairesFile, MADRS_sub_CID_list] = load_questionnaires_data;

%% load interleukins
[IL_data] = load_IL(study_nm);

%% extract the data for the selected subjects
[IL1b, IL6, IL18, ILsum, MADRS_S_score] = deal(NaN(1, NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    IL_sub_idx = strcmp(['CID',sub_nm], IL_data.CID);
    MADRS_sub_idx = strcmp(sub_nm, MADRS_sub_CID_list);

    IL1b(iS) = IL_data.IL1b(IL_sub_idx);
    IL6(iS) = IL_data.IL6(IL_sub_idx);
    IL18(iS) = IL_data.IL18(IL_sub_idx);
    ILsum(iS) = IL1b(iS) + IL6(iS) + IL18(iS);
    MADRS_S_score(iS) = excelReadQuestionnairesFile.MADRS_SCorrected(MADRS_sub_idx);
end % subject loop

%% correlate with each other
IL1b_goodSubs = ~isnan(IL1b).*~isnan(MADRS_S_score) == 1;
[b_IL1b,~,stats_tmp] = glmfit(IL1b(IL1b_goodSubs),...
    MADRS_S_score(IL1b_goodSubs), 'normal');
pval.IL1b = stats_tmp.p;
IL1b_fit = sort(IL1b(IL1b_goodSubs));
IL1b_prm_fit = glmval(b_IL1b, IL1b_fit, 'identity');

IL6_goodSubs = ~isnan(IL6).*~isnan(MADRS_S_score) == 1;
[b_IL6,~,stats_tmp] = glmfit(IL1b(IL6_goodSubs),...
    MADRS_S_score(IL6_goodSubs), 'normal');
pval.IL6 = stats_tmp.p;
IL6_fit = sort(IL6(IL6_goodSubs));
IL6_prm_fit = glmval(b_IL6, IL6_fit, 'identity');

IL18_goodSubs = ~isnan(IL18).*~isnan(MADRS_S_score) == 1;
[b_IL18,~,stats_tmp] = glmfit(IL18(IL18_goodSubs),...
    MADRS_S_score(IL18_goodSubs), 'normal');
pval.IL18 = stats_tmp.p;
IL18_fit = sort(IL18(IL18_goodSubs));
IL18_prm_fit = glmval(b_IL18, IL18_fit, 'identity');

ILsum_goodSubs = ~isnan(ILsum).*~isnan(MADRS_S_score) == 1;
[b_ILsum,~,stats_tmp] = glmfit(ILsum(ILsum_goodSubs),...
    MADRS_S_score(ILsum_goodSubs), 'normal');
pval.IL_sum = stats_tmp.p;
ILsum_fit = sort(ILsum(ILsum_goodSubs));
ILsum_prm_fit = glmval(b_ILsum, ILsum_fit, 'identity');

%% test splitting MADRS-S into 3 groups (<4 healthy; [4-6] low/mild depressed; >= 6 depressed)
MADRS_S_idx.low = MADRS_S_score < 4;
MADRS_S_idx.mid = (MADRS_S_score >= 4) & (MADRS_S_score < 6);
MADRS_S_idx.high = MADRS_S_score >= 6;

% average the data for each group
MADRS_S_groups = {'low','mid','high'};
n_MADRS_S = length(MADRS_S_groups);
for iMADRS = 1:n_MADRS_S
    MADRS_S_grp_nm = MADRS_S_groups{iMADRS};
    MADRS_S_grp_nm_bis = [MADRS_S_grp_nm,'_MADRS_S'];
    sub_idx = MADRS_S_idx.(MADRS_S_grp_nm);
    [m_IL1b.(MADRS_S_grp_nm_bis),...
        sem_IL1b.(MADRS_S_grp_nm_bis)] = mean_sem_sd(IL1b(sub_idx), 2);
    [m_IL6.(MADRS_S_grp_nm_bis),...
        sem_IL6.(MADRS_S_grp_nm_bis)] = mean_sem_sd(IL6(sub_idx), 2);
    [m_IL18.(MADRS_S_grp_nm_bis),...
        sem_IL18.(MADRS_S_grp_nm_bis)] = mean_sem_sd(IL18(sub_idx), 2);
    [m_ILsum.(MADRS_S_grp_nm_bis),...
        sem_ILsum.(MADRS_S_grp_nm_bis)] = mean_sem_sd(ILsum(sub_idx), 2);
end

%% figures
pSize = 50;
lWidth = 3;
black = [0 0 0];
grey = [143 143 143]./255;
%% test linear regression
fig;

% IL 1b
subplot(2,2,1);
hold on;
scat_hdl = scatter(IL1b(IL1b_goodSubs),...
    MADRS_S_score(IL1b_goodSubs));
fit_hdl = plot(IL1b_fit, IL1b_prm_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
ylabel('MADRS-S');
xlabel('IL-1b (pg/mL)');
legend_size(pSize);

% IL 6
subplot(2,2,2);
hold on;
scat_hdl = scatter(IL6(IL6_goodSubs),...
    MADRS_S_score(IL6_goodSubs));
fit_hdl = plot(IL6_fit, IL6_prm_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
ylabel('MADRS-S');
xlabel('IL-6 (pg/mL)');
legend_size(pSize);

% IL 18
subplot(2,2,3);
hold on;
scat_hdl = scatter(IL18(IL18_goodSubs),...
    MADRS_S_score(IL18_goodSubs));
fit_hdl = plot(IL18_fit, IL18_prm_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
ylabel('MADRS-S');
xlabel('IL-18 (pg/mL)');
legend_size(pSize);

% sum IL
subplot(2,2,4);
hold on;
scat_hdl = scatter(ILsum(ILsum_goodSubs),...
    MADRS_S_score(ILsum_goodSubs));
fit_hdl = plot(ILsum_fit, ILsum_prm_fit);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = black;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
ylabel('MADRS-S');
xlabel('IL-1b + IL6 + IL18 (pg/mL)');
legend_size(pSize);

%% test split per MADRS-S level
fig;

% IL 1b
subplot(2,2,1);
hold on;
err_hdl = errorbar(1:3,...
    [m_IL1b.low_MADRS_S,...
    m_IL1b.mid_MADRS_S,...
    m_IL1b.high_MADRS_S],...
    [sem_IL1b.low_MADRS_S,...
    sem_IL1b.mid_MADRS_S,...
    sem_IL1b.high_MADRS_S]);
err_hdl.LineWidth = lWidth;
xticks(1:3);
xticklabels({'low','mid','high'});
xlabel('MADRS-S');
ylabel('IL-1b (pg/mL)');
legend_size(pSize);

% IL 6
subplot(2,2,2);
hold on;
err_hdl = errorbar(1:3,...
    [m_IL6.low_MADRS_S,...
    m_IL6.mid_MADRS_S,...
    m_IL6.high_MADRS_S],...
    [sem_IL6.low_MADRS_S,...
    sem_IL6.mid_MADRS_S,...
    sem_IL6.high_MADRS_S]);
err_hdl.LineWidth = lWidth;
xticks(1:3);
xticklabels({'low','mid','high'});
xlabel('MADRS-S');
ylabel('IL-6 (pg/mL)');
legend_size(pSize);

% IL 18
subplot(2,2,3);
hold on;
err_hdl = errorbar(1:3,...
    [m_IL18.low_MADRS_S,...
    m_IL18.mid_MADRS_S,...
    m_IL18.high_MADRS_S],...
    [sem_IL18.low_MADRS_S,...
    sem_IL18.mid_MADRS_S,...
    sem_IL18.high_MADRS_S]);
err_hdl.LineWidth = lWidth;
xticks(1:3);
xticklabels({'low','mid','high'});
xlabel('MADRS-S');
ylabel('IL-18 (pg/mL)');
legend_size(pSize);

% sum IL
subplot(2,2,4);
hold on;
err_hdl = errorbar(1:3,...
    [m_ILsum.low_MADRS_S,...
    m_ILsum.mid_MADRS_S,...
    m_ILsum.high_MADRS_S],...
    [sem_ILsum.low_MADRS_S,...
    sem_ILsum.mid_MADRS_S,...
    sem_ILsum.high_MADRS_S]);
err_hdl.LineWidth = lWidth;
xticks(1:3);
xticklabels({'low','mid','high'});
xlabel('MADRS-S');
ylabel('IL-1b + IL6 + IL18 (pg/mL)');
legend_size(pSize);