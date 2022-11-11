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
[excelReadGeneralFile] = load_gal_data_bis;

%% load cortisol
[CORT_data] = load_CORT;

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

%% average
mCORT = mean(CORT,2,'omitnan');
[mStressRtg,...
    semStressRtg] = mean_sem_sd(stressRtg,2);
mStressRtg_fit = mean(stressRtg_f_CORT_fit,2,'omitnan');

%% are the betas significant?
[~,pval.stressRtg_f_CORT] = ttest(betas.stressRtg_f_CORT(2,:));

%% figures
pSize = 50;
lWidth = 3;
black = [0 0 0];
grey = [143 143 143]./255;

%% test linear regression
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

end % function