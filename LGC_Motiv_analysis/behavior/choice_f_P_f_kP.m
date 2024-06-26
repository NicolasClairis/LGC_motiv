function[choice_hEhP] = choice_f_P_f_kP(study_nm, subject_id, condition)
%[choice_hE] = choice_f_P_f_kP(study_nm, subject_id, condition)
% choice_f_P_f_kP will look at the average proportion of choices per
% punishment level and also by splitting on the kP parameter from the model
% looking at both global choices and choices in each task separately.
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% OUTPUTS
% choice_hE: structure with proportion of high-effort choices per subject

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
fig_choice_f_P_f_kP_violin = fig;
fig_choice_f_P_f_kP_mean = fig;
% fig_choice_f_P_f_kP_all = fig;
% fig_choice_f_P_f_kP_low_kP = fig;
fig_choice_f_P_f_kP_3_groups_mean = fig;
n_hP_levels = 3;
x_hP = 1:n_hP_levels;
nGroups = 3; % divide in 3 groups (low/medium/high)

% median split colours
low_col1 =[127 191 123]./255;
high_col1 =[175 141 195]./255;
% low/medium/high colours
low_col2 =[239 237 245]./255;
med_col2 =[188 189 220]./255;
high_col2 =[106 81 103]./255;
pSize = 30;
lWidth = 3;

%% load parameter
prm = prm_extraction(study_nm, subject_id, 'bayesian', '3');

%% extract parameter
prm_allSubs = prm.kP;

%% extract data
for iT = 1:nTasks
    task_nm = tasks{iT};
    
    %% extract proportion of choices
    [choice_hEhP_tmp] = extract_proportion_choice_hP_perSub(study_nm, subject_id, condition, task_nm);
    choice_hEhP.(task_nm) = choice_hEhP_tmp.(task_nm);
    
    %% 1) median split on the data
    [low_prm_subs1, high_prm_subs1] = medSplit_prm(study_nm, subject_id, kP);
    choice_hEhP.(task_nm).mSplit.low_prm.subject_id = subject_id(low_prm_subs1);
    choice_hEhP.(task_nm).mSplit.high_prm.subject_id = subject_id(high_prm_subs1);
    %% 2) perform 3 groups low/medium/high
    [~, ~, group_idx] = do_bin2(prm_allSubs, prm_allSubs, nGroups, 0);
    choice_hEhP.(task_nm).groups.low_prm.subject_id = subject_id(group_idx == 1);
    choice_hEhP.(task_nm).groups.med_prm.subject_id = subject_id(group_idx == 2);
    choice_hEhP.(task_nm).groups.high_prm.subject_id = subject_id(group_idx == 3);
    
    %% split people according to parameter
    % median split
    choice_hEhP.(task_nm).mSplit.(['low_',kP]) = choice_hEhP.(task_nm).allSubs(:,low_prm_subs1);
    choice_hEhP.(task_nm).mSplit.(['high_',kP]) = choice_hEhP.(task_nm).allSubs(:,high_prm_subs1);
    % more groups
    choice_hEhP.(task_nm).groups.(['low_',kP]) = choice_hEhP.(task_nm).allSubs(:,group_idx == 1);
    choice_hEhP.(task_nm).groups.(['med_',kP]) = choice_hEhP.(task_nm).allSubs(:,group_idx == 2);
    choice_hEhP.(task_nm).groups.(['high_',kP]) = choice_hEhP.(task_nm).allSubs(:,group_idx == 3);
    %% average
    % median split
    [m_choice_hEhP.(task_nm).mSplit.(['low_',kP]),...
        sem_choice_hEhP.(task_nm).mSplit.(['low_',kP])] = mean_sem_sd(choice_hEhP.(task_nm).mSplit.(['low_',kP]),2);
    [m_choice_hEhP.(task_nm).mSplit.(['high_',kP]),...
        sem_choice_hEhP.(task_nm).mSplit.(['high_',kP])] = mean_sem_sd(choice_hEhP.(task_nm).mSplit.(['high_',kP]),2);
    % more groups
    [m_choice_hEhP.(task_nm).groups.(['low_',kP]),...
        sem_choice_hEhP.(task_nm).groups.(['low_',kP])] = mean_sem_sd(choice_hEhP.(task_nm).groups.(['low_',kP]),2);
    [m_choice_hEhP.(task_nm).groups.(['med_',kP]),...
        sem_choice_hEhP.(task_nm).groups.(['med_',kP])] = mean_sem_sd(choice_hEhP.(task_nm).groups.(['med_',kP]),2);
    [m_choice_hEhP.(task_nm).groups.(['high_',kP]),...
        sem_choice_hEhP.(task_nm).groups.(['high_',kP])] = mean_sem_sd(choice_hEhP.(task_nm).groups.(['high_',kP]),2);
    
    %% figure with violin plot for low vs high median split
    figure(fig_choice_f_P_f_kP_violin);
    subplot(1,2,iT);
    hold on;
    if sum(low_prm_subs1) == sum(high_prm_subs1)
        choice_violin_mtrx = [choice_hEhP.(task_nm).mSplit.(['low_',kP])(1,:)',...
            choice_hEhP.(task_nm).mSplit.(['high_',kP])(1,:)',...
            choice_hEhP.(task_nm).mSplit.(['low_',kP])(2,:)',...
            choice_hEhP.(task_nm).mSplit.(['high_',kP])(2,:)',...
            choice_hEhP.(task_nm).mSplit.(['low_',kP])(3,:)',...
            choice_hEhP.(task_nm).mSplit.(['high_',kP])(3,:)'];
    elseif sum(low_prm_subs1) == sum(high_prm_subs1) + 1 % in case median split on number that is not odd
        choice_violin_mtrx = [choice_hEhP.(task_nm).mSplit.(['low_',kP])(1,:)',...
            [choice_hEhP.(task_nm).mSplit.(['high_',kP])(1,:),NaN]',...
            choice_hEhP.(task_nm).mSplit.(['low_',kP])(2,:)',...
            [choice_hEhP.(task_nm).mSplit.(['high_',kP])(2,:),NaN]',...
            choice_hEhP.(task_nm).mSplit.(['low_',kP])(3,:)',...
            [choice_hEhP.(task_nm).mSplit.(['high_',kP])(3,:),NaN]'];
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
    figure(fig_choice_f_P_f_kP_mean);
    subplot(1,2,iT);
    hold on;
    low_hdl = errorbar(x_hP-0.02,...
        m_choice_hEhP.(task_nm).mSplit.(['low_',kP]),...
        sem_choice_hEhP.(task_nm).mSplit.(['low_',kP]));
    low_hdl.LineWidth = lWidth;
    low_hdl.Color = low_col1;
    low_hdl.LineStyle = '-';
    high_hdl = errorbar(x_hP+0.02,...
        m_choice_hEhP.(task_nm).mSplit.(['high_',kP]),...
        sem_choice_hEhP.(task_nm).mSplit.(['high_',kP]));
    high_hdl.LineWidth = lWidth;
    high_hdl.Color = high_col1;
    high_hdl.LineStyle = '-';
    legend([high_hdl, low_hdl],{['high ',kP],['low ',kP]});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    xticks(1:n_hP_levels);
    xlabel('Punishment level');
    legend_size(pSize);

%% figure with mean for 3 groups (low/medium/high parameter)
figure(fig_choice_f_P_f_kP_3_groups_mean);
subplot(1,2,iT);
hold on;
low_hdl = errorbar(x_hP-0.02,...
    m_choice_hEhP.(task_nm).groups.(['low_',kP]),...
    sem_choice_hEhP.(task_nm).groups.(['low_',kP]));
low_hdl.LineWidth = lWidth;
low_hdl.Color = low_col2;
low_hdl.LineStyle = '-';
med_hdl = errorbar(x_hP-0.02,...
    m_choice_hEhP.(task_nm).groups.(['med_',kP]),...
    sem_choice_hEhP.(task_nm).groups.(['med_',kP]));
med_hdl.LineWidth = lWidth;
med_hdl.Color = med_col2;
med_hdl.LineStyle = '-';
high_hdl = errorbar(x_hP+0.02,...
    m_choice_hEhP.(task_nm).groups.(['high_',kP]),...
    sem_choice_hEhP.(task_nm).groups.(['high_',kP]));
high_hdl.LineWidth = lWidth;
high_hdl.Color = high_col2;
high_hdl.LineStyle = '-';
legend([high_hdl, med_hdl, low_hdl],...
    {['high ',kP],['medium ',kP],['low ',kP]});
legend('boxoff');
legend('Location','NorthEast');
ylim([0.2 1]);
ylabel([task_nm,' choices (%)']);
xticks(1:n_hP_levels);
xlabel('Punishment level');
legend_size(pSize);
end % task loop

end % function