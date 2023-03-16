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
    [subject_id, NS] = LGCM_subject_selection(study_nm,condition);
else
    NS = length(subject_id);
end
%% main parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
fig_choice_f_E_f_kE_violin = fig;
fig_choice_f_E_f_kE_mean = fig;
fig_choice_f_E_f_kE_all = fig;
fig_choice_f_E_f_kE_low_kE = fig;
n_hE_levels = 3;
x_hE = 1:n_hE_levels;

low_col =[127 191 123]./255;
high_col =[175 141 195]./255;
pSize = 30;
lWidth = 3;

%% extract data
for iT = 1:nTasks
    task_nm = tasks{iT};
    
    %% extract proportion of choices
    [choice_hE_tmp] = extract_proportion_choice_hE_perSub(study_nm, subject_id, condition, task_nm);
    choice_hE.(task_nm) = choice_hE_tmp.(task_nm);
    %% extract parameter
    prm_nm = ['k',task_nm];
    [low_prm_subs, high_prm_subs] = medSplit_prm(study_nm, subject_id, prm_nm);
    choice_hE.(task_nm).low_prm.subject_id = subject_id(low_prm_subs);
    choice_hE.(task_nm).high_prm.subject_id = subject_id(high_prm_subs);
    
    %% split people according to parameter
    choice_hE.(task_nm).(['low_',prm_nm]) = choice_hE.(task_nm).allSubs(:,low_prm_subs);
    choice_hE.(task_nm).(['high_',prm_nm]) = choice_hE.(task_nm).allSubs(:,high_prm_subs);
    %% average
    [m_choice_hE.(task_nm).(['low_',prm_nm]),...
        sem_choice_hE.(task_nm).(['low_',prm_nm])] = mean_sem_sd(choice_hE.(task_nm).(['low_',prm_nm]),2);
    [m_choice_hE.(task_nm).(['high_',prm_nm]),...
        sem_choice_hE.(task_nm).(['high_',prm_nm])] = mean_sem_sd(choice_hE.(task_nm).(['high_',prm_nm]),2);
    
    %% figure with violin plot
    figure(fig_choice_f_E_f_kE_violin);
    subplot(1,2,iT);
    hold on;
    choice_violin_mtrx = [choice_hE.(task_nm).(['low_',prm_nm])(1,:)',...
        choice_hE.(task_nm).(['high_',prm_nm])(1,:)',...
        choice_hE.(task_nm).(['low_',prm_nm])(2,:)',...
        choice_hE.(task_nm).(['high_',prm_nm])(2,:)',...
        choice_hE.(task_nm).(['low_',prm_nm])(3,:)',...
        choice_hE.(task_nm).(['high_',prm_nm])(3,:)'];
    % replace 0 by eps to avoid bugs
    choice_violin_mtrx(choice_violin_mtrx == 0) = eps;
    violinplot(choice_violin_mtrx,...
        {'E1_l','E1_h','E2_l','E2_h','E3_l','E3_h'},...
        'ViolinColor',[low_col;high_col;low_col;high_col;low_col;high_col]);
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);
    
    %% figure with mean
    figure(fig_choice_f_E_f_kE_mean);
    subplot(1,2,iT);
    hold on;
    low_hdl = errorbar(x_hE-0.02,...
        m_choice_hE.(task_nm).(['low_',prm_nm]),...
        sem_choice_hE.(task_nm).(['low_',prm_nm]));
    low_hdl.LineWidth = lWidth;
    low_hdl.Color = low_col;
    low_hdl.LineStyle = '-';
    high_hdl = errorbar(x_hE+0.02,...
        m_choice_hE.(task_nm).(['high_',prm_nm]),...
        sem_choice_hE.(task_nm).(['high_',prm_nm]));
    high_hdl.LineWidth = lWidth;
    high_hdl.Color = high_col;
    high_hdl.LineStyle = '-';
    legend([high_hdl, low_hdl],{['high ',prm_nm],['low ',prm_nm]});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);
    
    %% figure with everybody
    figure(fig_choice_f_E_f_kE_all);
    subplot(1,2,iT);
    hold on;
    for iS = 1:NS
        plot_hdl = plot(x_hE, choice_hE.(task_nm).allSubs(:,iS));
        plot_hdl.LineWidth = 1;
        plot_hdl.LineStyle = '-';
    end % loop through subjects
    ylim([0 1]);
    xticks(1:3);
    xticklabels({'E1','E2','E3'});
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);
    
    %% figure only with low kE
    figure(fig_choice_f_E_f_kE_low_kE);
    subplot(1,2,iT);
    title(['low ',prm_nm]);
    hold on;
    for iS = 1:size(choice_hE.(task_nm).(['high_',prm_nm]),2)
        plot_hdl = plot(x_hE, choice_hE.(task_nm).(['high_',prm_nm])(:,iS));
        plot_hdl.LineWidth = 1;
        plot_hdl.LineStyle = '-';
    end % loop through subjects
    ylim([0 1]);
    xticks(1:3);
    xticklabels({'E1','E2','E3'});
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);
end % task loop

end % function