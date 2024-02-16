function[choice_hE] = choice_f_E_f_kE(study_nm, subject_id, condition)
%[choice_hE] = choice_f_E_f_kE(study_nm, subject_id, condition)
% choice_f_E_f_kE will look at the average proportion of choices per effort
% level and also by splitting on the kEp and kEm parameters from the model
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% OUTPUTS
% choice_hE: structure with proportion of choice per subject

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
tasks = {'Ep','Em'};
nTasks = length(tasks);
fig_choice_f_E_f_kE_violin = fig;
fig_choice_f_E_f_kE_mean = fig;
% fig_choice_f_E_f_kE_all = fig;
% fig_choice_f_E_f_kE_low_kE = fig;
fig_choice_f_E_f_kE_3_groups_mean = fig;
n_hE_levels = 3;
x_hE = 1:n_hE_levels;
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

%% extract data
for iT = 1:nTasks
    task_nm = tasks{iT};
    
    %% extract proportion of choices
    [choice_hE_tmp] = extract_proportion_choice_hE_perSub(study_nm, subject_id, condition, task_nm);
    choice_hE.(task_nm) = choice_hE_tmp.(task_nm);
    %% extract parameter
    prm_nm = ['k',task_nm];
    prm_allSubs = prm.(prm_nm);
    
    %% 1) median split on the data
    [low_prm_subs1, high_prm_subs1] = medSplit_prm(study_nm, subject_id, prm_nm);
    choice_hE.(task_nm).mSplit.low_prm.subject_id = subject_id(low_prm_subs1);
    choice_hE.(task_nm).mSplit.high_prm.subject_id = subject_id(high_prm_subs1);
    %% 2) perform 3 groups low/medium/high
    [~, ~, group_idx] = do_bin2(prm_allSubs, prm_allSubs, nGroups, 0);
    choice_hE.(task_nm).groups.low_prm.subject_id = subject_id(group_idx == 1);
    choice_hE.(task_nm).groups.med_prm.subject_id = subject_id(group_idx == 2);
    choice_hE.(task_nm).groups.high_prm.subject_id = subject_id(group_idx == 3);
    
    %% split people according to parameter
    % median split
    choice_hE.(task_nm).mSplit.(['low_',prm_nm]) = choice_hE.(task_nm).allSubs(:,low_prm_subs1);
    choice_hE.(task_nm).mSplit.(['high_',prm_nm]) = choice_hE.(task_nm).allSubs(:,high_prm_subs1);
    % more groups
    choice_hE.(task_nm).groups.(['low_',prm_nm]) = choice_hE.(task_nm).allSubs(:,group_idx == 1);
    choice_hE.(task_nm).groups.(['med_',prm_nm]) = choice_hE.(task_nm).allSubs(:,group_idx == 2);
    choice_hE.(task_nm).groups.(['high_',prm_nm]) = choice_hE.(task_nm).allSubs(:,group_idx == 3);
    %% average
    % median split
    [m_choice_hE.(task_nm).mSplit.(['low_',prm_nm]),...
        sem_choice_hE.(task_nm).mSplit.(['low_',prm_nm])] = mean_sem_sd(choice_hE.(task_nm).mSplit.(['low_',prm_nm]),2);
    [m_choice_hE.(task_nm).mSplit.(['high_',prm_nm]),...
        sem_choice_hE.(task_nm).mSplit.(['high_',prm_nm])] = mean_sem_sd(choice_hE.(task_nm).mSplit.(['high_',prm_nm]),2);
    % more groups
    [m_choice_hE.(task_nm).groups.(['low_',prm_nm]),...
        sem_choice_hE.(task_nm).groups.(['low_',prm_nm])] = mean_sem_sd(choice_hE.(task_nm).groups.(['low_',prm_nm]),2);
    [m_choice_hE.(task_nm).groups.(['med_',prm_nm]),...
        sem_choice_hE.(task_nm).groups.(['med_',prm_nm])] = mean_sem_sd(choice_hE.(task_nm).groups.(['med_',prm_nm]),2);
    [m_choice_hE.(task_nm).groups.(['high_',prm_nm]),...
        sem_choice_hE.(task_nm).groups.(['high_',prm_nm])] = mean_sem_sd(choice_hE.(task_nm).groups.(['high_',prm_nm]),2);
    
    %% figure with violin plot for low vs high median split
    figure(fig_choice_f_E_f_kE_violin);
    subplot(1,2,iT);
    hold on;
    if sum(low_prm_subs1) == sum(high_prm_subs1)
        choice_violin_mtrx = [choice_hE.(task_nm).mSplit.(['low_',prm_nm])(1,:)',...
            choice_hE.(task_nm).mSplit.(['high_',prm_nm])(1,:)',...
            choice_hE.(task_nm).mSplit.(['low_',prm_nm])(2,:)',...
            choice_hE.(task_nm).mSplit.(['high_',prm_nm])(2,:)',...
            choice_hE.(task_nm).mSplit.(['low_',prm_nm])(3,:)',...
            choice_hE.(task_nm).mSplit.(['high_',prm_nm])(3,:)'];
    elseif sum(low_prm_subs1) == sum(high_prm_subs1) + 1 % in case median split on number that is not odd
        choice_violin_mtrx = [choice_hE.(task_nm).mSplit.(['low_',prm_nm])(1,:)',...
            [choice_hE.(task_nm).mSplit.(['high_',prm_nm])(1,:),NaN]',...
            choice_hE.(task_nm).mSplit.(['low_',prm_nm])(2,:)',...
            [choice_hE.(task_nm).mSplit.(['high_',prm_nm])(2,:),NaN]',...
            choice_hE.(task_nm).mSplit.(['low_',prm_nm])(3,:)',...
            [choice_hE.(task_nm).mSplit.(['high_',prm_nm])(3,:),NaN]'];
    end
    % replace 0 by eps to avoid bugs
    choice_violin_mtrx(choice_violin_mtrx == 0) = eps;
    violinplot(choice_violin_mtrx,...
        {'E1_l','E1_h','E2_l','E2_h','E3_l','E3_h'},...
        'ViolinColor',[low_col1;high_col1;low_col1;high_col1;low_col1;high_col1]);
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);
    
    %% figure with average for low vs high median split
    figure(fig_choice_f_E_f_kE_mean);
    subplot(1,2,iT);
    hold on;
    low_hdl = errorbar(x_hE-0.02,...
        m_choice_hE.(task_nm).mSplit.(['low_',prm_nm]),...
        sem_choice_hE.(task_nm).mSplit.(['low_',prm_nm]));
    low_hdl.LineWidth = lWidth;
    low_hdl.Color = low_col1;
    low_hdl.LineStyle = '-';
    high_hdl = errorbar(x_hE+0.02,...
        m_choice_hE.(task_nm).mSplit.(['high_',prm_nm]),...
        sem_choice_hE.(task_nm).mSplit.(['high_',prm_nm]));
    high_hdl.LineWidth = lWidth;
    high_hdl.Color = high_col1;
    high_hdl.LineStyle = '-';
    legend([high_hdl, low_hdl],{['high ',prm_nm],['low ',prm_nm]});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    xticks(1:n_hE_levels);
    xlabel('Effort level');
    legend_size(pSize);
    
%     %% figure with everybody
%     figure(fig_choice_f_E_f_kE_all);
%     subplot(1,2,iT);
%     hold on;
%     for iS = 1:NS
%         plot_hdl = plot(x_hE, choice_hE.(task_nm).allSubs(:,iS));
%         plot_hdl.LineWidth = 1;
%         plot_hdl.LineStyle = '-';
%     end % loop through subjects
%     ylim([0 1]);
%     xticks(1:3);
%     xticklabels({'E1','E2','E3'});
%     ylabel([task_nm,' choices (%)']);
%     legend_size(pSize);
    
%     %% figure only with low kE
%     figure(fig_choice_f_E_f_kE_low_kE);
%     subplot(1,2,iT);
%     title(['low ',prm_nm]);
%     hold on;
%     for iS = 1:size(choice_hE.(task_nm).(['high_',prm_nm]),2)
%         plot_hdl = plot(x_hE, choice_hE.(task_nm).(['high_',prm_nm])(:,iS));
%         plot_hdl.LineWidth = 1;
%         plot_hdl.LineStyle = '-';
%     end % loop through subjects
%     ylim([0 1]);
%     xticks(1:3);
%     xticklabels({'E1','E2','E3'});
%     ylabel([task_nm,' choices (%)']);
%     legend_size(pSize);

%% figure with mean for 3 groups (low/medium/high parameter)
figure(fig_choice_f_E_f_kE_3_groups_mean);
subplot(1,2,iT);
hold on;
low_hdl = errorbar(x_hE-0.02,...
    m_choice_hE.(task_nm).groups.(['low_',prm_nm]),...
    sem_choice_hE.(task_nm).groups.(['low_',prm_nm]));
low_hdl.LineWidth = lWidth;
low_hdl.Color = low_col2;
low_hdl.LineStyle = '-';
med_hdl = errorbar(x_hE-0.02,...
    m_choice_hE.(task_nm).groups.(['med_',prm_nm]),...
    sem_choice_hE.(task_nm).groups.(['med_',prm_nm]));
med_hdl.LineWidth = lWidth;
med_hdl.Color = med_col2;
med_hdl.LineStyle = '-';
high_hdl = errorbar(x_hE+0.02,...
    m_choice_hE.(task_nm).groups.(['high_',prm_nm]),...
    sem_choice_hE.(task_nm).groups.(['high_',prm_nm]));
high_hdl.LineWidth = lWidth;
high_hdl.Color = high_col2;
high_hdl.LineStyle = '-';
legend([high_hdl, med_hdl, low_hdl],...
    {['high ',prm_nm],['medium ',prm_nm],['low ',prm_nm]});
legend('boxoff');
legend('Location','NorthEast');
ylim([0.2 1]);
ylabel([task_nm,' choices (%)']);
xticks(1:n_hE_levels);
xlabel('Effort level');
legend_size(pSize);
end % task loop

end % function