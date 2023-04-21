function[NS] = choice_f_E_fkE_and_kF(study_nm)
% [NS] = choice_f_E_fkE_and_kF(study_nm)
% choice_f_E_fkE_and_kF will split the choices according to effort level,
% kEp and kFp (and kEm, kFm respectively) in order to see whether there is 
% an interaction between kE and kF.
%
% INPUTS
% study_nm: study name
%
% OUTPUTS
% NS: number of subjects used

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% main parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
n_hE_levels = 3;
x_hE = 1:n_hE_levels;
low_kE_low_kF_col =[223 194 125]./255;
low_kE_high_kF_col =[166 97 26]./255;
high_kE_low_kF_col =[128 205 193]./255;
high_kE_high_kF_col =[1 133 113]./255;
pSize = 30;
lWidth = 3;


%% extract choices = f(E) and f(kE)
[choice_hE] = choice_f_E_f_kE(study_nm, subject_id, condition);
close all;
fig_choice_f_E_f_kE_kF_mean = fig;
fig_choice_f_E_f_kE_kF_mean_focus_low_kE = fig;
%% loop through tasks
for iT = 1:nTasks
    task_nm = tasks{iT};
    %% extract time effect
    switch task_nm
        case 'Ep'
            kE_prm_nm = 'kEp';
            time_prm_nm = 'kFp';
            time_prm_nm_bis = 'kFp';
        case 'Em'
            kE_prm_nm = 'kEm';
            time_prm_nm = 'kLm';
            time_prm_nm_bis = 'kFm';
    end
    
    %% extract index of subject depending on low vs high kE and kF
    low_kE_subs = choice_hE.(task_nm).low_prm.subject_id;
    high_kE_subs = choice_hE.(task_nm).high_prm.subject_id;
    [low_kE_low_kF_subs, low_kE_high_kF_subs] = medSplit_prm(study_nm, low_kE_subs, time_prm_nm);
    [high_kE_low_kF_subs, high_kE_high_kF_subs] = medSplit_prm(study_nm, high_kE_subs, time_prm_nm);
    
    %% extract data for subjects depending on kE/kF
    choice_hE.(task_nm).low_kE_low_kF = choice_hE.(task_nm).(['low_',kE_prm_nm])(:,low_kE_low_kF_subs);
    choice_hE.(task_nm).low_kE_high_kF = choice_hE.(task_nm).(['low_',kE_prm_nm])(:,low_kE_high_kF_subs);
    choice_hE.(task_nm).high_kE_low_kF = choice_hE.(task_nm).(['high_',kE_prm_nm])(:,high_kE_low_kF_subs);
    choice_hE.(task_nm).high_kE_high_kF = choice_hE.(task_nm).(['high_',kE_prm_nm])(:,high_kE_high_kF_subs);
    
    %% average data within each group
    [m_choice_hE.(task_nm).low_kE_low_kF,...
        sem_choice_hE.(task_nm).low_kE_low_kF] = mean_sem_sd(choice_hE.(task_nm).low_kE_low_kF,2);
    [m_choice_hE.(task_nm).low_kE_high_kF,...
        sem_choice_hE.(task_nm).low_kE_high_kF] = mean_sem_sd(choice_hE.(task_nm).low_kE_high_kF,2);
    [m_choice_hE.(task_nm).high_kE_low_kF,...
        sem_choice_hE.(task_nm).high_kE_low_kF] = mean_sem_sd(choice_hE.(task_nm).high_kE_low_kF,2);
    [m_choice_hE.(task_nm).high_kE_high_kF,...
        sem_choice_hE.(task_nm).high_kE_high_kF] = mean_sem_sd(choice_hE.(task_nm).high_kE_high_kF,2);
    
    %% figures
    %% figure with mean
    figure(fig_choice_f_E_f_kE_kF_mean);
    subplot(1,2,iT);
    hold on;
    low_kE_low_kF_hdl = errorbar(x_hE-0.02,...
        m_choice_hE.(task_nm).low_kE_low_kF,...
        sem_choice_hE.(task_nm).low_kE_low_kF);
    low_kE_low_kF_hdl.LineWidth = lWidth;
    low_kE_low_kF_hdl.Color = low_kE_low_kF_col;
    low_kE_low_kF_hdl.LineStyle = '-';
    
    
    low_kE_high_kF_hdl = errorbar(x_hE+0.02,...
        m_choice_hE.(task_nm).low_kE_high_kF,...
        sem_choice_hE.(task_nm).low_kE_high_kF);
    low_kE_high_kF_hdl.LineWidth = lWidth;
    low_kE_high_kF_hdl.Color = low_kE_high_kF_col;
    low_kE_high_kF_hdl.LineStyle = '-';
    
    high_kE_low_kF_hdl = errorbar(x_hE-0.02,...
        m_choice_hE.(task_nm).high_kE_low_kF,...
        sem_choice_hE.(task_nm).high_kE_low_kF);
    high_kE_low_kF_hdl.LineWidth = lWidth;
    high_kE_low_kF_hdl.Color = high_kE_low_kF_col;
    high_kE_low_kF_hdl.LineStyle = '-';
    
    high_kE_high_kF_hdl = errorbar(x_hE+0.02,...
        m_choice_hE.(task_nm).high_kE_high_kF,...
        sem_choice_hE.(task_nm).high_kE_high_kF);
    high_kE_high_kF_hdl.LineWidth = lWidth;
    high_kE_high_kF_hdl.Color = high_kE_high_kF_col;
    high_kE_high_kF_hdl.LineStyle = '-';
    
    legend([low_kE_high_kF_hdl, low_kE_low_kF_hdl,...
        high_kE_high_kF_hdl, high_kE_low_kF_hdl],...
        {['low ',kE_prm_nm,' high ',time_prm_nm_bis],...
        ['low ',kE_prm_nm,' low ',time_prm_nm_bis],...
        ['high ',kE_prm_nm,' high ',time_prm_nm_bis],...
        ['high ',kE_prm_nm,' low ',time_prm_nm_bis]});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    xticks(1:n_hE_levels);
    xlabel('high E level');
    legend_size(pSize);
    
    %% figure focusing on low kE
    figure(fig_choice_f_E_f_kE_kF_mean_focus_low_kE);
    subplot(1,2,iT);
    hold on;
    low_kE_low_kF_hdl = errorbar(x_hE-0.02,...
        m_choice_hE.(task_nm).low_kE_low_kF,...
        sem_choice_hE.(task_nm).low_kE_low_kF);
    low_kE_low_kF_hdl.LineWidth = lWidth;
    low_kE_low_kF_hdl.Color = low_kE_low_kF_col;
    low_kE_low_kF_hdl.LineStyle = '-';
    
    
    low_kE_high_kF_hdl = errorbar(x_hE+0.02,...
        m_choice_hE.(task_nm).low_kE_high_kF,...
        sem_choice_hE.(task_nm).low_kE_high_kF);
    low_kE_high_kF_hdl.LineWidth = lWidth;
    low_kE_high_kF_hdl.Color = low_kE_high_kF_col;
    low_kE_high_kF_hdl.LineStyle = '-';
    
    legend([low_kE_high_kF_hdl, low_kE_low_kF_hdl],...
        {['low ',kE_prm_nm,' high ',time_prm_nm_bis],...
        ['low ',kE_prm_nm,' low ',time_prm_nm_bis]});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0 1]);
    xticks(1:n_hE_levels);
    xlabel('high E level');
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);
end % task loop

end % function