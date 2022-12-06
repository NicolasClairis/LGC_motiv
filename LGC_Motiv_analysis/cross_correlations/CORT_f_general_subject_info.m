
%% define subjects of interest
condition = subject_condition;
study_nm = 'study1';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load cortisol
[CORT_data] = load_CORT(study_nm);

%% load age and gender information
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
gal_sub_List = excelReadGeneralFile.CID;
age_list = excelReadGeneralFile.Age_yearsOld_;
gender_list = strcmp(excelReadGeneralFile.Sexe_femaleF_maleM_,'F');
BMI_list = excelReadGeneralFile.BMI;

%% extract correspondency
[CORT.AUCg, CORT.avgCORT,...
    age, gender, BMI] = deal(NaN(1,NS));
for iS = 1:NS
    % subject info
    sub_nm = ['CID',subject_id{iS}];
    sub_age_idx = strcmp(sub_nm, gal_sub_List);
    sub_CORT_idx = strcmp(sub_nm, CORT_data.CID);
    % extract general informations
    age(iS) = age_list(sub_age_idx);
    gender(iS) = gender_list(sub_age_idx);
    BMI(iS) = BMI_list(sub_age_idx);
    % extract CORT data
    CORT.AUCg(iS) = CORT_data.AUCg(sub_CORT_idx);
    CORT.avgCORT(iS) = mean(CORT_data.CORT(:,sub_CORT_idx),1,'omitnan');
end % subject loop

%% test correlations
% age = f(CORT AUC)
goodSubs.age_f_AUCg = (~isnan(CORT.AUCg)).*(~isnan(age)) == 1;
[betas.age_f_AUCg,~,stats.age_f_AUCg] = glmfit(CORT.AUCg(goodSubs.age_f_AUCg),...
    age(goodSubs.age_f_AUCg),'normal');
pval.age_f_AUCg = stats.age_f_AUC.p;
CORT_fit.age_f_AUCg = sort(CORT.AUCg(goodSubs.age_f_AUCg));
age_fit.age_f_AUCg = glmval(betas.age_f_AUCg,...
    CORT_fit.age_f_AUCg, 'identity');

% BMI = f(CORT AUC)
goodSubs.BMI_f_AUCg = (~isnan(CORT.AUCg)).*(~isnan(BMI)) == 1;
[betas.BMI_f_AUCg,~,stats.BMI_f_AUCg] = glmfit(CORT.AUCg(goodSubs.BMI_f_AUCg),...
    BMI(goodSubs.BMI_f_AUCg),'normal');
pval.BMI_f_AUCg = stats.BMI_f_AUCg.p;
CORT_fit.BMI_f_AUCg = sort(CORT.AUCg(goodSubs.BMI_f_AUCg));
BMI_fit.BMI_f_AUCg = glmval(betas.BMI_f_AUCg,...
    CORT_fit.BMI_f_AUCg, 'identity');

% age = f(avg(CORT))
goodSubs.age_f_avgCORT = (~isnan(CORT.avgCORT)).*(~isnan(age)) == 1;
[betas.age_f_avgCORT,~,stats.age_f_avgCORT] = glmfit(CORT.avgCORT(goodSubs.age_f_avgCORT),...
    age(goodSubs.age_f_avgCORT),'normal');
pval.age_f_avgCORT = stats.age_f_avgCORT.p;
CORT_fit.age_f_avgCORT = sort(CORT.avgCORT(goodSubs.age_f_avgCORT));
age_fit.age_f_avgCORT = glmval(betas.age_f_avgCORT,...
    CORT_fit.age_f_avgCORT, 'identity');

% BMI = f(avg(CORT))
goodSubs.BMI_f_avgCORT = (~isnan(CORT.avgCORT)).*(~isnan(BMI)) == 1;
[betas.BMI_f_avgCORT,~,stats.BMI_f_avgCORT] = glmfit(CORT.avgCORT(goodSubs.BMI_f_avgCORT),...
    BMI(goodSubs.BMI_f_avgCORT),'normal');
pval.BMI_f_avgCORT = stats.BMI_f_avgCORT.p;
CORT_fit.BMI_f_avgCORT = sort(CORT.avgCORT(goodSubs.BMI_f_avgCORT));
BMI_fit.BMI_f_avgCORT = glmval(betas.BMI_f_avgCORT,...
    CORT_fit.BMI_f_avgCORT, 'identity');

% extract results with gender
[~,pval.gender_f_AUCg] = ttest2(CORT.AUCg(gender == 1),...
    CORT.AUCg(gender == 0));
[~,pval.gender_f_avgCORT] = ttest2(CORT.avgCORT(gender == 1),...
    CORT.avgCORT(gender == 0));

%% display figures
lWidth = 3;
pSize = 30;
grey = [143 143 143]./255;

fig;
%% age
% AUC
subplot(2,3,1);
hold on;
scat_hdl = scatter(CORT.AUCg, age);
fit_hdl = plot(CORT_fit.age_f_AUCg, age_fit.age_f_AUCg);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
xlabel('Cortisol AUC');
ylabel('Age (y.o.)');
legend_size(pSize);

% average CORT
subplot(2,3,4);
hold on;
scat_hdl = scatter(CORT.avgCORT, age);
fit_hdl = plot(CORT_fit.age_f_avgCORT, age_fit.age_f_avgCORT);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
xlabel('Cortisol (μg/dL)');
ylabel('Age (y.o.)');
legend_size(pSize);

%% BMI
% AUC
subplot(2,3,2);
hold on;
scat_hdl = scatter(CORT.AUCg, BMI);
fit_hdl = plot(CORT_fit.BMI_f_AUCg, BMI_fit.BMI_f_AUCg);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
xlabel('Cortisol AUC');
ylabel('BMI');
legend_size(pSize);

% average CORT
subplot(2,3,5);
hold on;
scat_hdl = scatter(CORT.avgCORT, BMI);
fit_hdl = plot(CORT_fit.BMI_f_avgCORT, BMI_fit.BMI_f_avgCORT);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
xlabel('Cortisol (μg/dL)');
ylabel('BMI');
legend_size(pSize);

%% gender
% AUC
subplot(2,3,3);
hold on;
Violin({CORT.AUCg(gender == 1)},1);
Violin({CORT.AUCg(gender == 0)},2);
ylabel('Cortisol AUC');
xticks(1:2);
xticklabels({'females','males'});
legend_size(pSize);

% average CORT
subplot(2,3,6);
hold on;
Violin({CORT.avgCORT(gender == 1)},1);
Violin({CORT.avgCORT(gender == 0)},2);
ylabel('Cortisol (μg/dL)');
xticks(1:2);
xticklabels({'females','males'});
legend_size(pSize);