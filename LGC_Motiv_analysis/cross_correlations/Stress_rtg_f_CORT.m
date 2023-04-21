function[betas, pval] = Stress_rtg_f_CORT()
% [betas, pval] = Stress_rtg_f_CORT()
% Stress_rtg_f_CORT will extract the stress levels and the cortisol
% of each timepoint.
%
% INPUTS
%
% OUTPUTS
% betas: individual beta for the slope between stress rating and cortisol
%
% pval: p.value for an unpaired t.test against zero of the betas

%% load subject
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load stress questionnaires
[excelReadGeneralFile] = load_gal_data_bis(study_nm);

%% load cortisol
[CORT_data] = load_CORT(study_nm);

%% extract the data for the selected subjects
nTimePoints = 4;
[CORT, stressRtg, stressRtg_f_CORT_fit] = deal(NaN(nTimePoints, NS));
betas.stressRtg_f_CORT = NaN(2,NS);
for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    CORT_sub_idx = strcmp(sub_nm, CORT_data.CID);
    quest_sub_idx = strcmp(sub_nm, excelReadGeneralFile.CID);

    CORT(:,iS) = CORT_data.CORT(:,CORT_sub_idx);
    stressRtg(1,iS) = excelReadGeneralFile.StressPr__MRS(quest_sub_idx);
    stressRtg(2,iS) = excelReadGeneralFile.StressPost_MRS(quest_sub_idx);
    stressRtg(3,iS) = excelReadGeneralFile.StressPr__IRMf(quest_sub_idx);
    stressRtg(4,iS) = excelReadGeneralFile.StressPost_IRMf(quest_sub_idx);

    %% correlate data together if dataset is complete
    if sum(~isnan(CORT(:,iS))) == 4
        betas.stressRtg_f_CORT(:,iS) = glmfit(CORT(:,iS),...
            stressRtg(:,iS),...
            'normal');
        stressRtg_f_CORT_fit(:,iS) = glmval(betas.stressRtg_f_CORT(:,iS),...
            CORT(:,iS), 'identity');
    end
end % subject loop

%% try to perform a correlation with fixed effects 
%(ignoring inter-individual differences)
[megaCORT, megaStressRtg] = deal(NaN(4*NS,1));
for iS = 1:NS
    sub_idx = (1:4) + 4*(iS - 1);
    megaCORT(sub_idx) = CORT(:,iS);
    megaStressRtg(sub_idx) = stressRtg(:,iS);
end % subject loop
okTrials = (~isnan(megaCORT)).*(~isnan(megaStressRtg)) == 1;
[betas.stressRtg_f_CORT_allTogether,~,stats_tmp] = glmfit(megaCORT(okTrials),...
            megaStressRtg(okTrials),...
            'normal');
pval.stressRtg_f_CORT_allTogether = stats_tmp.p;
megaCort_ascending = sort(megaCORT(okTrials));
fitStressRtg_f_megaCort = glmval(betas.stressRtg_f_CORT_allTogether,...
    megaCort_ascending, 'identity');
        
%% correlation across individuals averaging data per subject
mCort_perSub = mean(CORT,1,'omitnan');
mStressRtg_perSub = mean(stressRtg,1,'omitnan');
okSubs = (~isnan(mCort_perSub)).*(~isnan(mStressRtg_perSub)) == 1;
[betas.avg_stressRtg_f_avg_CORT,~,stats_tmp] = glmfit(mCort_perSub(okSubs),...
            mStressRtg_perSub(okSubs),...
            'normal');
pval.avg_stressRtg_f_avg_CORT = stats_tmp.p;
mCort_perSub_ascending = sort(mCort_perSub(okSubs));
fitStressRtg_f_mCort_perSub = glmval(betas.avg_stressRtg_f_avg_CORT,...
    mCort_perSub_ascending, 'identity');

%% average
mCORT = mean(CORT,2,'omitnan');
[mStressRtg,...
    semStressRtg] = mean_sem_sd(stressRtg,2);
mStressRtg_fit = mean(stressRtg_f_CORT_fit,2,'omitnan');

%% are the betas significantly different from zero?
[~,pval.stressRtg_f_CORT] = ttest(betas.stressRtg_f_CORT(2,:));

%% figures
pSize = 50;
lWidth = 3;
black = [0 0 0];
grey = [143 143 143]./255;

%% test linear regression per subject averaged (random effects)
fig;

% STAI-T = f(CORT)
hold on;
plot_hdl = errorbar(mCORT, mStressRtg, semStressRtg);
fit_hdl = plot(mCORT, mStressRtg_fit);

plot_hdl.LineStyle = '--';
plot_hdl.LineWidth = lWidth;
plot_hdl.Color = black;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
xlabel('Cortisol');
ylabel('Stress rating');
legend_size(pSize);

%% test linear regression across all subjects (fixed effects)
fig;

% STAI-T = f(CORT)
hold on;
scatter_hdl = scatter(megaCORT, megaStressRtg);
fit_hdl = plot(megaCort_ascending, fitStressRtg_f_megaCort);

scatter_hdl.LineWidth = lWidth;
scatter_hdl.MarkerEdgeColor = black;
fit_hdl.LineStyle = '--';
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
xlabel('Cortisol (µg/dL)');
ylabel('Stress rating');
legend_size(pSize);

%% test linear regression across all subjects (fixed effects) bis
fig;

% STAI-T = f(CORT)
hold on;
scatter_hdl = scatter(mCort_perSub, mStressRtg_perSub);
fit_hdl = plot(mCort_perSub_ascending, fitStressRtg_f_mCort_perSub);

scatter_hdl.LineWidth = lWidth;
scatter_hdl.MarkerEdgeColor = black;
fit_hdl.LineStyle = '--';
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
xlabel('Cortisol (µg/dL)');
ylabel('Stress rating');
legend_size(pSize);

end % function