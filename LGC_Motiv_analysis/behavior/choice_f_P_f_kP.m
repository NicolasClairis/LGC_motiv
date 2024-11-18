function[choice_hEhP, kBias] = choice_f_P_f_kP(study_nm, subject_id, condition)
%[choice_hE, kBias] = choice_f_P_f_kP(study_nm, subject_id, condition)
% choice_f_P_f_kP will look at the average proportion of choices per
% punishment level and also by splitting on the kP parameter from the model
% looking at both global choices and choices in each task separately.
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% OUTPUTS
% choice_hE: structure with proportion of high-effort choices per subject
%
% kBias: structure with bias parameter averaged according to the different
% grouping methods (to see if there is any change with the kP parameter

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id] = LGCM_subject_selection(study_nm,condition);
end
%% main parameters
tasks = {'Ep','Em','all'};
nTasks = length(tasks);
n_hP_levels = 3;
x_hP = 1:n_hP_levels;
nGroups = 3; % divide in 3 groups (low/medium/high)
nGroups2 = 4; % divide in 4 groups (very low/low/high/very high)

% median split colours
low_col1 =[127 191 123]./255;
high_col1 =[175 141 195]./255;
% 3 groups (low/medium/high) colours
low_col2 =[239 237 245]./255;
med_col2 =[188 189 220]./255;
high_col2 =[106 81 103]./255;
% 4 groups (low/medium/high) colours
very_low_col3   = [242 240 247]./255;
low_col3        = [203 201 226]./255;
high_col3       = [158 154 200]./255;
very_high_col3  = [106 81 163]./255;
pSize = 30;
lWidth = 3;

%% load parameter
prm = prm_extraction(study_nm, subject_id, 'bayesian');

%% extract parameter
kP_allSubs = prm.kP;
kBias_allSubs = prm.kBias;

%% initialize figures
% figure/task
fig_choice_f_P_f_kP_violin_perTask = fig;
fig_choice_f_P_f_kP_mean_perTask = fig;
fig_choice_f_P_f_kP_3_groups_mean_perTask = fig;
fig_choice_f_P_f_kP_4_groups_mean_perTask = fig;
% figure pooling tasks
fig_choice_f_P_f_kP_violin = fig;
fig_choice_f_P_f_kP_mean = fig;
fig_choice_f_P_f_kP_3_groups_mean = fig;
fig_choice_f_P_f_kP_4_groups_mean = fig;

%% extract data
for iT = 1:nTasks
    task_nm = tasks{iT};
    
    %% extract proportion of choices
    [choice_hEhP_tmp] = extract_proportion_choice_hP_perSub(study_nm, subject_id, condition);
    choice_hEhP.(task_nm) = choice_hEhP_tmp.(task_nm);
    
    %% 1) median split on the data
    [low_kP_subs1, high_kP_subs1] = medSplit_prm(study_nm, subject_id, 'kP');
    choice_hEhP.(task_nm).mSplit.low_prm.subject_id = subject_id(low_kP_subs1);
    choice_hEhP.(task_nm).mSplit.high_prm.subject_id = subject_id(high_kP_subs1);

    %% 2) perform 3 groups low/medium/high
    [~, ~, group_idx1] = do_bin2(kP_allSubs, kP_allSubs, nGroups, 0);
    low_kP_group1 = group_idx1 == 1;
    med_kP_group1 = group_idx1 == 2;
    high_kP_group1 = group_idx1 == 3;
    choice_hEhP.(task_nm).groups1.low_prm.subject_id = subject_id(low_kP_group1);
    choice_hEhP.(task_nm).groups1.med_prm.subject_id = subject_id(med_kP_group1);
    choice_hEhP.(task_nm).groups1.high_prm.subject_id = subject_id(high_kP_group1);

    %% 3) perform 4 groups very low/low/high/very high
    [~, ~, group_idx2] = do_bin2(kP_allSubs, kP_allSubs, nGroups2, 0);
    very_low_kP_group2 = group_idx2 == 1;
    low_kP_group2 = group_idx2 == 2;
    high_kP_group2 = group_idx2 == 3;
    very_high_kP_group2 = group_idx2 == 4;
    choice_hEhP.(task_nm).groups2.very_low_prm.subject_id = subject_id(very_low_kP_group2);
    choice_hEhP.(task_nm).groups2.low_prm.subject_id = subject_id(low_kP_group2);
    choice_hEhP.(task_nm).groups2.high_prm.subject_id = subject_id(high_kP_group2);
    choice_hEhP.(task_nm).groups2.very_high_prm.subject_id = subject_id(very_high_kP_group2);

    %% split people according to parameter
    % median split
    choice_hEhP.(task_nm).mSplit.low_kP = choice_hEhP.(task_nm).allSubs(:,low_kP_subs1);
    choice_hEhP.(task_nm).mSplit.high_kP = choice_hEhP.(task_nm).allSubs(:,high_kP_subs1);
    kBias.mSplit.low_kP.allSubs = kBias_allSubs(low_kP_subs1);
    kBias.mSplit.high_kP.allSubs = kBias_allSubs(high_kP_subs1);

    % 3 groups
    choice_hEhP.(task_nm).groups1.low_kP = choice_hEhP.(task_nm).allSubs(:,low_kP_group1);
    choice_hEhP.(task_nm).groups1.med_kP = choice_hEhP.(task_nm).allSubs(:,med_kP_group1);
    choice_hEhP.(task_nm).groups1.high_kP = choice_hEhP.(task_nm).allSubs(:,high_kP_group1);
    kBias.groups1.low_kP.allSubs = kBias_allSubs(low_kP_group1);
    kBias.groups1.med_kP.allSubs = kBias_allSubs(med_kP_group1);
    kBias.groups1.high_kP.allSubs = kBias_allSubs(high_kP_group1);
    % 4 groups
    choice_hEhP.(task_nm).groups2.very_low_kP = choice_hEhP.(task_nm).allSubs(:,very_low_kP_group2);
    choice_hEhP.(task_nm).groups2.low_kP = choice_hEhP.(task_nm).allSubs(:,low_kP_group2);
    choice_hEhP.(task_nm).groups2.high_kP = choice_hEhP.(task_nm).allSubs(:,high_kP_group2);
    choice_hEhP.(task_nm).groups2.very_high_kP = choice_hEhP.(task_nm).allSubs(:,very_high_kP_group2);
    kBias.groups2.very_low_kP.allSubs = kBias_allSubs(very_low_kP_group2);
    kBias.groups2.low_kP.allSubs = kBias_allSubs(low_kP_group2);
    kBias.groups2.high_kP.allSubs = kBias_allSubs(high_kP_group2);
    kBias.groups2.very_high_kP.allSubs = kBias_allSubs(very_high_kP_group2);
    %% average
    % median split
    [m_choice_hEhP.(task_nm).mSplit.low_kP,...
        sem_choice_hEhP.(task_nm).mSplit.low_kP] = mean_sem_sd(choice_hEhP.(task_nm).mSplit.low_kP,2);
    [m_choice_hEhP.(task_nm).mSplit.high_kP,...
        sem_choice_hEhP.(task_nm).mSplit.high_kP] = mean_sem_sd(choice_hEhP.(task_nm).mSplit.high_kP,2);
    [kBias.mSplit.low_kP.mean,...
        kBias.mSplit.low_kP.sem] = mean_sem_sd(kBias.mSplit.low_kP.allSubs,2);
    [kBias.mSplit.high_kP.mean,...
        kBias.mSplit.high_kP.sem] = mean_sem_sd(kBias.mSplit.high_kP.allSubs,2);
    % 3 groups
    [m_choice_hEhP.(task_nm).groups1.low_kP,...
        sem_choice_hEhP.(task_nm).groups1.low_kP] = mean_sem_sd(choice_hEhP.(task_nm).groups1.low_kP,2);
    [m_choice_hEhP.(task_nm).groups1.med_kP,...
        sem_choice_hEhP.(task_nm).groups1.med_kP] = mean_sem_sd(choice_hEhP.(task_nm).groups1.med_kP,2);
    [m_choice_hEhP.(task_nm).groups1.high_kP,...
        sem_choice_hEhP.(task_nm).groups1.high_kP] = mean_sem_sd(choice_hEhP.(task_nm).groups1.high_kP,2);
    [kBias.groups1.low_kP.mean,...
        kBias.groups1.low_kP.sem] = mean_sem_sd(kBias.groups1.low_kP.allSubs,2);
    [kBias.groups1.med_kP.mean,...
        kBias.groups1.med_kP.sem] = mean_sem_sd(kBias.groups1.med_kP.allSubs,2);
    [kBias.groups1.high_kP.mean,...
        kBias.groups1.high_kP.sem] = mean_sem_sd(kBias.groups1.high_kP.allSubs,2);
    % 4 groups
    [m_choice_hEhP.(task_nm).groups2.very_low_kP,...
        sem_choice_hEhP.(task_nm).groups2.very_low_kP] = mean_sem_sd(choice_hEhP.(task_nm).groups2.very_low_kP,2);
    [m_choice_hEhP.(task_nm).groups2.low_kP,...
        sem_choice_hEhP.(task_nm).groups2.low_kP] = mean_sem_sd(choice_hEhP.(task_nm).groups2.low_kP,2);
    [m_choice_hEhP.(task_nm).groups2.high_kP,...
        sem_choice_hEhP.(task_nm).groups2.high_kP] = mean_sem_sd(choice_hEhP.(task_nm).groups2.high_kP,2);
    [m_choice_hEhP.(task_nm).groups2.very_high_kP,...
        sem_choice_hEhP.(task_nm).groups2.very_high_kP] = mean_sem_sd(choice_hEhP.(task_nm).groups2.very_high_kP,2);
    [kBias.groups2.very_low_kP.mean,...
        kBias.groups2.very_low_kP.sem] = mean_sem_sd(kBias.groups2.very_low_kP.allSubs,2);
    [kBias.groups2.low_kP.mean,...
        kBias.groups2.low_kP.sem] = mean_sem_sd(kBias.groups2.low_kP.allSubs,2);
    [kBias.groups2.high_kP.mean,...
        kBias.groups2.high_kP.sem] = mean_sem_sd(kBias.groups2.high_kP.allSubs,2);
    [kBias.groups2.very_high_kP.mean,...
        kBias.groups2.very_high_kP.sem] = mean_sem_sd(kBias.groups2.very_high_kP.allSubs,2);
    %% figure with violin plot for low vs high median split
    switch task_nm
        case {'Ep','Em'}
            figure(fig_choice_f_P_f_kP_violin_perTask);
            subplot(1,2,iT);
        case 'all'
            figure(fig_choice_f_P_f_kP_violin);
    end
    hold on;
    if sum(low_kP_subs1) == sum(high_kP_subs1)
        choice_violin_mtrx = [choice_hEhP.(task_nm).mSplit.low_kP(1,:)',...
            choice_hEhP.(task_nm).mSplit.high_kP(1,:)',...
            choice_hEhP.(task_nm).mSplit.low_kP(2,:)',...
            choice_hEhP.(task_nm).mSplit.high_kP(2,:)',...
            choice_hEhP.(task_nm).mSplit.low_kP(3,:)',...
            choice_hEhP.(task_nm).mSplit.high_kP(3,:)'];
    elseif sum(low_kP_subs1) == sum(high_kP_subs1) + 1 % in case median split on number that is not odd
        choice_violin_mtrx = [choice_hEhP.(task_nm).mSplit.low_kP(1,:)',...
            [choice_hEhP.(task_nm).mSplit.high_kP(1,:),NaN]',...
            choice_hEhP.(task_nm).mSplit.low_kP(2,:)',...
            [choice_hEhP.(task_nm).mSplit.high_kP(2,:),NaN]',...
            choice_hEhP.(task_nm).mSplit.low_kP(3,:)',...
            [choice_hEhP.(task_nm).mSplit.high_kP(3,:),NaN]'];
    end
    % replace 0 by eps to avoid bugs
    choice_violin_mtrx(choice_violin_mtrx == 0) = eps;
    violinplot(choice_violin_mtrx,...
        {'P1_l','P1_h','P2_l','P2_h','P3_l','P3_h'},...
        'ViolinColor',[low_col1;high_col1;low_col1;high_col1;low_col1;high_col1]);
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);

    %% figure with average for low vs high median split
    switch task_nm
        case {'Ep','Em'}
            figure(fig_choice_f_P_f_kP_mean_perTask);
            subplot(1,2,iT);
        case 'all'
            figure(fig_choice_f_P_f_kP_mean);
    end
    hold on;
    very_low_hdl = errorbar(x_hP-0.02,...
        m_choice_hEhP.(task_nm).mSplit.low_kP,...
        sem_choice_hEhP.(task_nm).mSplit.low_kP);
    very_low_hdl.LineWidth = lWidth;
    very_low_hdl.Color = low_col1;
    very_low_hdl.LineStyle = '-';
    high_hdl = errorbar(x_hP+0.02,...
        m_choice_hEhP.(task_nm).mSplit.high_kP,...
        sem_choice_hEhP.(task_nm).mSplit.high_kP);
    high_hdl.LineWidth = lWidth;
    high_hdl.Color = high_col1;
    high_hdl.LineStyle = '-';
    legend([high_hdl, very_low_hdl],{'high kP','low kP'});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    xticks(1:n_hP_levels);
    xlabel('Punishment level');
    legend_size(pSize);

    %% figure with mean for 3 groups (low/medium/high parameter)
    switch task_nm
        case {'Ep','Em'}
            figure(fig_choice_f_P_f_kP_3_groups_mean_perTask);
            subplot(1,2,iT);
        case 'all'
            figure(fig_choice_f_P_f_kP_3_groups_mean);
    end
    hold on;
    very_low_hdl = errorbar(x_hP-0.02,...
        m_choice_hEhP.(task_nm).groups1.low_kP,...
        sem_choice_hEhP.(task_nm).groups1.low_kP);
    very_low_hdl.LineWidth = lWidth;
    very_low_hdl.Color = low_col2;
    very_low_hdl.LineStyle = '-';
    med_hdl = errorbar(x_hP-0.02,...
        m_choice_hEhP.(task_nm).groups1.med_kP,...
        sem_choice_hEhP.(task_nm).groups1.med_kP);
    med_hdl.LineWidth = lWidth;
    med_hdl.Color = med_col2;
    med_hdl.LineStyle = '-';
    high_hdl = errorbar(x_hP+0.02,...
        m_choice_hEhP.(task_nm).groups1.high_kP,...
        sem_choice_hEhP.(task_nm).groups1.high_kP);
    high_hdl.LineWidth = lWidth;
    high_hdl.Color = high_col2;
    high_hdl.LineStyle = '-';
    legend([high_hdl, med_hdl, very_low_hdl],...
        {'high kP','medium kP','low kP'});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0.2 1]);
    ylabel([task_nm,' choices (%)']);
    xticks(1:n_hP_levels);
    xlabel('Punishment level');
    legend_size(pSize);
    %% figure with mean for 4 groups (very low/low/high/very high parameter)
    switch task_nm
        case {'Ep','Em'}
            figure(fig_choice_f_P_f_kP_4_groups_mean_perTask);
            subplot(1,2,iT);
        case 'all'
            figure(fig_choice_f_P_f_kP_4_groups_mean);
    end
    hold on;
    very_low_hdl = errorbar(x_hP-0.02,...
        m_choice_hEhP.(task_nm).groups2.very_low_kP,...
        sem_choice_hEhP.(task_nm).groups2.very_low_kP);
    very_low_hdl.LineWidth = lWidth;
    very_low_hdl.Color = very_low_col3;
    very_low_hdl.LineStyle = '-';
    low_hdl = errorbar(x_hP-0.02,...
        m_choice_hEhP.(task_nm).groups2.low_kP,...
        sem_choice_hEhP.(task_nm).groups2.low_kP);
    low_hdl.LineWidth = lWidth;
    low_hdl.Color = low_col3;
    low_hdl.LineStyle = '-';
    high_hdl = errorbar(x_hP+0.02,...
        m_choice_hEhP.(task_nm).groups2.high_kP,...
        sem_choice_hEhP.(task_nm).groups2.high_kP);
    high_hdl.LineWidth = lWidth;
    high_hdl.Color = high_col3;
    high_hdl.LineStyle = '-';
    very_high_hdl = errorbar(x_hP+0.02,...
        m_choice_hEhP.(task_nm).groups2.very_high_kP,...
        sem_choice_hEhP.(task_nm).groups2.very_high_kP);
    very_high_hdl.LineWidth = lWidth;
    very_high_hdl.Color = very_high_col3;
    very_high_hdl.LineStyle = '-';
    legend([very_high_hdl, high_hdl, low_hdl, very_low_hdl],...
        {'very high kP','high kP','low kP','very low kP'});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0.2 1]);
    ylabel([task_nm,' choices (%)']);
    xticks(1:n_hP_levels);
    xlabel('Punishment level');
    legend_size(pSize);

end % task loop
%% show kBias according to different groups
fig;
% median split
subplot(1,3,1);
bar_hdl = bar(1:2,[kBias.mSplit.low_kP.mean, kBias.mSplit.high_kP.mean]);
erbar_hdl = errorbar(1:2,...
    [kBias.mSplit.low_kP.mean, kBias.mSplit.high_kP.mean],...
    [kBias.mSplit.low_kP.sem, kBias.mSplit.high_kP.sem]);
errorbar_hdl_upgrade(erbar_hdl, 'k');
xticks(1:2);
xticklabels({'low','high'});
xlim([0.8 2.2]);
xlabel('kP');
ylabel('kBias');


% 3 groups
subplot(1,3,2);
bar_hdl = bar(1:3,[kBias.groups1.low_kP.mean, kBias.groups1.med_kP.mean, kBias.groups1.high_kP.mean]);
erbar_hdl = errorbar(1:3,...
    [kBias.groups1.low_kP.mean, kBias.groups1.med_kP.mean, kBias.groups1.high_kP.mean],...
    [kBias.groups1.low_kP.sem, kBias.groups1.med_kP.sem, kBias.groups1.high_kP.sem]);
errorbar_hdl_upgrade(erbar_hdl, 'k');
xticks(1:3);
xticklabels({'low','med','high'});
xlim([0.8 3.2]);
xlabel('kP');
ylabel('kBias');

% 4 groups
subplot(1,3,3);
bar_hdl = bar(1:4,[kBias.groups2.very_low_kP.mean, kBias.groups2.low_kP.mean, kBias.groups2.high_kP.mean, kBias.groups2.very_high_kP.mean]);
erbar_hdl = errorbar(1:4,...
    [kBias.groups2.very_low_kP.mean, kBias.groups2.low_kP.mean, kBias.groups2.high_kP.mean, kBias.groups2.very_high_kP.mean],...
    [kBias.groups2.very_low_kP.sem, kBias.groups2.low_kP.sem, kBias.groups2.high_kP.sem, kBias.groups2.very_high_kP.sem]);
errorbar_hdl_upgrade(erbar_hdl, 'k');
xticks(1:4);
xticklabels({'vlow','low','high','vhigh'});
xlim([0.8 4.2]);
xlabel('kP');
ylabel('kBias');

end % function