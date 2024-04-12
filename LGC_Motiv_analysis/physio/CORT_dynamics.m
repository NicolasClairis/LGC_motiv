%% load subject
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load cortisol
[CORT_data] = load_CORT(study_nm);

%% extract the data for the selected subjects
nTimePoints = 4;
[CORT, time,...
    time_sorted, CORT_fit] = deal(NaN(nTimePoints, NS));
betas = NaN(2, NS);

for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    CORT_sub_idx = strcmp(sub_nm, CORT_data.CID);

    CORT(:,iS) = CORT_data.CORT(:,CORT_sub_idx);
    time(:,iS) = CORT_data.timings(:,CORT_sub_idx);
    
    % GLM
    [~, betas(:,iS), ~, ~,...
        time_sorted(:,iS), CORT_fit(:,iS)] = glm_package(time(:,iS), CORT(:,iS), 'normal');
end % subject loop

%% test if betas significantly different from zero
[~,pval_betas] = ttest(betas(2,:));

%% average
[mCort, semCORT] = mean_sem_sd(CORT,2);
[mTime, semTime] = mean_sem_sd(time,2);
[mTime_sorted, semTime_sorted] = mean_sem_sd(time_sorted,2);
[mCort_fit, semCORT_fit] = mean_sem_sd(CORT_fit,2);

%% plot
pSize = 50;
[~, lWidth, col] = general_fig_prm;

fig;
er_hdl = errorbar(mTime, mCort, semCORT, semCORT, semTime, semTime);
er_hdl.LineWidth = lWidth;
% er_hdl.LineStyle = ':';
er_hdl.LineStyle = 'none';
er_hdl.Color = col.grey;
fit_hdl = plot(mTime_sorted, mCort_fit, 'LineStyle','-');
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = col.black;
xlabel('Time of the day (hours)');
ylabel('Cortisol (µg/dL)');