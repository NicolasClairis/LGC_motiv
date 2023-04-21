function[NS] = choice_f_E_fkE_and_bias(study_nm)
% [NS] = choice_f_E_fkE_and_bias(study_nm)
% choice_f_E_fkE_and_bias will split the choices according to effort level,
% kEp and bias (and kEm with bias respectively) in order to see whether there is 
% an interaction between kE and bias.
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
low_kE_low_bias_col =[223 194 125]./255;
low_kE_high_bias_col =[166 97 26]./255;
high_kE_low_bias_col =[128 205 193]./255;
high_kE_high_bias_col =[1 133 113]./255;
pSize = 30;
lWidth = 3;


%% extract choices = f(E) and f(kE)
[choice_hE] = choice_f_E_f_kE(study_nm, subject_id, condition);
close all;
fig_choice_f_E_f_kE_bias_mean = fig;
fig_choice_f_E_f_kE_bias_mean_focus_low_kE = fig;
%% loop through tasks
for iT = 1:nTasks
    task_nm = tasks{iT};
    %% extract time effect
    switch task_nm
        case 'Ep'
            kE_prm_nm = 'kEp';
        case 'Em'
            kE_prm_nm = 'kEm';
    end
    bias_prm_nm = 'kBiasM';
    bias_prm_nm_bis = 'bias';
    
    %% extract index of subject depending on low vs high kE and bias
    low_kE_subs = choice_hE.(task_nm).low_prm.subject_id;
    high_kE_subs = choice_hE.(task_nm).high_prm.subject_id;
    [low_kE_low_bias_subs, low_kE_high_bias_subs] = medSplit_prm(study_nm, low_kE_subs, bias_prm_nm);
    [high_kE_low_bias_subs, high_kE_high_bias_subs] = medSplit_prm(study_nm, high_kE_subs, bias_prm_nm);
    
    %% extract data for subjects depending on kE/bias
    choice_hE.(task_nm).low_kE_low_bias = choice_hE.(task_nm).(['low_',kE_prm_nm])(:,low_kE_low_bias_subs);
    choice_hE.(task_nm).low_kE_high_bias = choice_hE.(task_nm).(['low_',kE_prm_nm])(:,low_kE_high_bias_subs);
    choice_hE.(task_nm).high_kE_low_bias = choice_hE.(task_nm).(['high_',kE_prm_nm])(:,high_kE_low_bias_subs);
    choice_hE.(task_nm).high_kE_high_bias = choice_hE.(task_nm).(['high_',kE_prm_nm])(:,high_kE_high_bias_subs);
    
    %% average data within each group
    [m_choice_hE.(task_nm).low_kE_low_bias,...
        sem_choice_hE.(task_nm).low_kE_low_bias] = mean_sem_sd(choice_hE.(task_nm).low_kE_low_bias,2);
    [m_choice_hE.(task_nm).low_kE_high_bias,...
        sem_choice_hE.(task_nm).low_kE_high_bias] = mean_sem_sd(choice_hE.(task_nm).low_kE_high_bias,2);
    [m_choice_hE.(task_nm).high_kE_low_bias,...
        sem_choice_hE.(task_nm).high_kE_low_bias] = mean_sem_sd(choice_hE.(task_nm).high_kE_low_bias,2);
    [m_choice_hE.(task_nm).high_kE_high_bias,...
        sem_choice_hE.(task_nm).high_kE_high_bias] = mean_sem_sd(choice_hE.(task_nm).high_kE_high_bias,2);
    
    %% figures
    %% figure with mean
    figure(fig_choice_f_E_f_kE_bias_mean);
    subplot(1,2,iT);
    hold on;
    low_kE_low_bias_hdl = errorbar(x_hE-0.02,...
        m_choice_hE.(task_nm).low_kE_low_bias,...
        sem_choice_hE.(task_nm).low_kE_low_bias);
    low_kE_low_bias_hdl.LineWidth = lWidth;
    low_kE_low_bias_hdl.Color = low_kE_low_bias_col;
    low_kE_low_bias_hdl.LineStyle = '-';
    
    
    low_kE_high_bias_hdl = errorbar(x_hE+0.02,...
        m_choice_hE.(task_nm).low_kE_high_bias,...
        sem_choice_hE.(task_nm).low_kE_high_bias);
    low_kE_high_bias_hdl.LineWidth = lWidth;
    low_kE_high_bias_hdl.Color = low_kE_high_bias_col;
    low_kE_high_bias_hdl.LineStyle = '-';
    
    high_kE_low_bias_hdl = errorbar(x_hE-0.02,...
        m_choice_hE.(task_nm).high_kE_low_bias,...
        sem_choice_hE.(task_nm).high_kE_low_bias);
    high_kE_low_bias_hdl.LineWidth = lWidth;
    high_kE_low_bias_hdl.Color = high_kE_low_bias_col;
    high_kE_low_bias_hdl.LineStyle = '-';
    
    high_kE_high_bias_hdl = errorbar(x_hE+0.02,...
        m_choice_hE.(task_nm).high_kE_high_bias,...
        sem_choice_hE.(task_nm).high_kE_high_bias);
    high_kE_high_bias_hdl.LineWidth = lWidth;
    high_kE_high_bias_hdl.Color = high_kE_high_bias_col;
    high_kE_high_bias_hdl.LineStyle = '-';
    
    legend([low_kE_high_bias_hdl, low_kE_low_bias_hdl,...
        high_kE_high_bias_hdl, high_kE_low_bias_hdl],...
        {['low ',kE_prm_nm,' high ',bias_prm_nm_bis],...
        ['low ',kE_prm_nm,' low ',bias_prm_nm_bis],...
        ['high ',kE_prm_nm,' high ',bias_prm_nm_bis],...
        ['high ',kE_prm_nm,' low ',bias_prm_nm_bis]});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    xticks(1:n_hE_levels);
    xlabel('high E level');
    legend_size(pSize);
    
    %% figure focusing on low kE
    figure(fig_choice_f_E_f_kE_bias_mean_focus_low_kE);
    subplot(1,2,iT);
    hold on;
    low_kE_low_bias_hdl = errorbar(x_hE-0.02,...
        m_choice_hE.(task_nm).low_kE_low_bias,...
        sem_choice_hE.(task_nm).low_kE_low_bias);
    low_kE_low_bias_hdl.LineWidth = lWidth;
    low_kE_low_bias_hdl.Color = low_kE_low_bias_col;
    low_kE_low_bias_hdl.LineStyle = '-';
    
    
    low_kE_high_bias_hdl = errorbar(x_hE+0.02,...
        m_choice_hE.(task_nm).low_kE_high_bias,...
        sem_choice_hE.(task_nm).low_kE_high_bias);
    low_kE_high_bias_hdl.LineWidth = lWidth;
    low_kE_high_bias_hdl.Color = low_kE_high_bias_col;
    low_kE_high_bias_hdl.LineStyle = '-';
    
    legend([low_kE_high_bias_hdl, low_kE_low_bias_hdl],...
        {['low ',kE_prm_nm,' high ',bias_prm_nm_bis],...
        ['low ',kE_prm_nm,' low ',bias_prm_nm_bis]});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0 1]);
    xticks(1:n_hE_levels);
    xlabel('high E level');
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);
end % task loop

end % function