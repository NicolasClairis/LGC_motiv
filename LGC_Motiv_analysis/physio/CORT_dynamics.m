%% load subject
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load cortisol
[CORT_data] = load_CORT(study_nm);

%% extract the data for the selected subjects
nTimePoints = 4;
[CORT, time] = deal(NaN(nTimePoints, NS));

for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    CORT_sub_idx = strcmp(sub_nm, CORT_data.CID);

    CORT(:,iS) = CORT_data.CORT(:,CORT_sub_idx);
    time(:,iS) = CORT_data.timings(:,CORT_sub_idx);
end % subject loop

%% average
[mCort, semCORT] = mean_sem_sd(CORT,2);
[mTime, semTime] = mean_sem_sd(time,2);

%% plot
lWidth = 3;
pSize = 50;

fig;
hdl = errorbar(mTime, mCort, semCORT, semCORT, semTime, semTime);
hdl.LineWidth = lWidth;
hdl.LineStyle = ':';
hdl.Color = 'k';
xlabel('Time of the day (hours)');
ylabel('CORTISOL (µg/dL)');
legend_size(pSize);