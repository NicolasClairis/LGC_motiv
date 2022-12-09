% test whether aversion to efforts extracted with the modeling is related
% to the JPI-R score (i.e. energy)

%% working directories
computerRoot = LGCM_root_paths;
%% define subjects
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load questionnaire scores
[excelReadQuestionnairesFile, quest_S_sub_CID_list] = load_questionnaires_data();
MADRS_S_fullList = excelReadQuestionnairesFile.MADRS_SCorrected;
JPI_R_fullList = excelReadQuestionnairesFile.JPI_RScore;
[MADRS_S_score, JPI_R_score] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    MADRS_S_score(iS) = MADRS_S_fullList(strcmp(quest_S_sub_CID_list, sub_nm));
    JPI_R_score(iS) = JPI_R_fullList(strcmp(quest_S_sub_CID_list, sub_nm));
end % subject list

%% load parameters
[mdlType, mdlN] = behavioral_model_selection;
behavioralPrm = prm_extraction(study_nm, subject_id,mdlType, mdlN);
kEp = behavioralPrm.kEp;
kEm = behavioralPrm.kEm;

%% perform correlation
goodSubs_MADRS_S = ~isnan(kEp.*kEm).*(~isnan(MADRS_S_score)) == 1;
goodSubs_JPI_R = ~isnan(kEp.*kEm).*(~isnan(JPI_R_score)) == 1;
[betas.kEp_MADRS,~,stats.kEp_MADRS] = glmfit(kEp(goodSubs_MADRS_S), MADRS_S_score(goodSubs_MADRS_S),'normal');
[betas.kEm_MADRS,~,stats.kEm_MADRS] = glmfit(kEm(goodSubs_MADRS_S), MADRS_S_score(goodSubs_MADRS_S),'normal');
[betas.kEp_JPI,~,stats.kEp_JPI] = glmfit(kEp(goodSubs_JPI_R), JPI_R_score(goodSubs_JPI_R),'normal');
[betas.kEm_JPI,~,stats.kEm_JPI] = glmfit(kEm(goodSubs_JPI_R), JPI_R_score(goodSubs_JPI_R),'normal');
% extract p.values
pval.kEp_MADRS = stats.kEp_MADRS.p;
pval.kEm_MADRS = stats.kEm_MADRS.p;
pval.kEp_JPIR = stats.kEp_JPIR.p;
pval.kEm_JPIR = stats.kEm_JPIR.p;
% extract a fit
kEp_MADRS_sort = sort(kEp(goodSubs_MADRS_S));
[MADRS_kEp_fit] = glmval(betas.kEp_MADRS, kEp_MADRS_sort, 'identity');
kEm_MADRS_sort = sort(kEm(goodSubs_MADRS_S));
[MADRS_kEm_fit] = glmval(betas.kEm_MADRS, kEm_MADRS_sort, 'identity');
kEp_JPI_sort = sort(kEp(goodSubs_JPI_R));
[JPI_kEp_fit] = glmval(betas.kEp_JPI, kEp_JPI_sort, 'identity');
kEm_JPI_sort = sort(kEm(goodSubs_JPI_R));
[JPI_kEm_fit] = glmval(betas.kEm_JPI, kEm_JPI_sort, 'identity');
%% display figure
lWidth = 3;
pSize = 30;
grey = [143 143 143]./255;

fig;
% Ep/MADRS-S
subplot(2,2,1);
hold on;
scat_hdl = scatter(kEp(goodSubs_MADRS_S),...
    MADRS_S_score(goodSubs_MADRS_S));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kEp_MADRS_sort,...
    MADRS_kEp_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel('kEp');
ylabel('MADRS-S');
legend_size(pSize);

% Em/MADRS-S
subplot(2,2,2);
hold on;
scat_hdl = scatter(kEm(goodSubs_MADRS_S),...
    MADRS_S_score(goodSubs_MADRS_S));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kEm_MADRS_sort,...
    MADRS_kEm_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel('kEm');
ylabel('MADRS-S');
legend_size(pSize);

% Ep/JPI-R
subplot(2,2,3);
hold on;
scat_hdl = scatter(kEp(goodSubs_JPI_R),...
    JPI_R_score(goodSubs_JPI_R));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kEp_JPI_sort,...
    JPI_kEp_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel('kEp');
ylabel('JPI-R');
legend_size(pSize);

% Em/JPI-R
subplot(2,2,4);
hold on;
scat_hdl = scatter(kEm(goodSubs_JPI_R),...
    JPI_R_score(goodSubs_JPI_R));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kEm_JPI_sort,...
    JPI_kEm_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel('kEm');
ylabel('JPI-R');
legend_size(pSize);