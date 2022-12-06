% script in order to extract different groups based on metabolite levels
% and see whether clear-cut clusters appear.

%% working directories
computerRoot = LGCM_root_paths;
%% define subjects
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load metabolites
[metabolite_allSubs, ROI_nm,...
    metabolite_nm] = metabolite_extraction(study_nm, subject_id);

%% load parameters
[mdlType, mdlN] = behavioral_model_selection;
behavioralPrm = prm_extraction(study_nm, subject_id,mdlType, mdlN);
behavioral_prm_names = fieldnames(behavioralPrm);
behavioral_prm_names(strcmp(behavioral_prm_names,'CID')) = [];
n_prm = length(behavioral_prm_names);

%% split data into three groups
nb_bins = 3;
[~, Bin_val, bin_idx] = do_bin2(1:length(metabolite_allSubs),...
    metabolite_allSubs, nb_bins, 0);
metab_low = bin_idx == 1;
metab_mid = bin_idx == 2;
metab_high = bin_idx == 3;

%% extract parameters for each group
for iPrm = 1:n_prm
    prm_nm = behavioral_prm_names{iPrm};
    % zscore each parameter across individuals
    norm_prm_tmp = nanzscore(behavioralPrm.(prm_nm));
%     norm_prm_tmp = behavioralPrm.(prm_nm);
    % split according to metabolite levels
    [m_prm.(prm_nm).(metabolite_nm).low,...
        sem_prm.(prm_nm).(metabolite_nm).low] = mean_sem_sd(norm_prm_tmp(metab_low),2);
    [m_prm.(prm_nm).(metabolite_nm).mid,...
        sem_prm.(prm_nm).(metabolite_nm).mid] = mean_sem_sd(norm_prm_tmp(metab_mid),2);
    [m_prm.(prm_nm).(metabolite_nm).high,...
        sem_prm.(prm_nm).(metabolite_nm).high] = mean_sem_sd(norm_prm_tmp(metab_high),2);
end % parameter loop

%% display figure
lWidth = 3;
pSize = 30;
%% MADRS-S
fig;
hold on;
x_main = 1:3; % 3 groups
xPrm = (-0.4):(1/(n_prm+1)):0.4;
bar_hdl_leg = NaN(1,n_prm);
for iPrm = 1:n_prm
    prm_nm = behavioral_prm_names{iPrm};
    bar_hdl = bar([x_main+xPrm(iPrm)],...
        [m_prm.(prm_nm).(metabolite_nm).low,...
        m_prm.(prm_nm).(metabolite_nm).mid,...
        m_prm.(prm_nm).(metabolite_nm).high]);
    bar_hdl.BarWidth = (1/(n_prm+1));
    bar_hdl_leg(iPrm) = bar_hdl;
    er_hdl = errorbar([x_main+xPrm(iPrm)],...
        [m_prm.(prm_nm).(metabolite_nm).low,...
        m_prm.(prm_nm).(metabolite_nm).mid,...
        m_prm.(prm_nm).(metabolite_nm).high],...
        [sem_prm.(prm_nm).(metabolite_nm).low,...
        sem_prm.(prm_nm).(metabolite_nm).mid,...
        sem_prm.(prm_nm).(metabolite_nm).high]);
    er_hdl.LineStyle = 'none';
    er_hdl.LineWidth = lWidth;
    er_hdl.Color = 'k';
end
legend(bar_hdl_leg, behavioral_prm_names);
legend('Location','northwest');
legend('boxoff');
xticks(x_main);
xticklabels({[metabolite_nm,' low'],...
    [metabolite_nm,' mid'],...
    [metabolite_nm,' high']});
xlim([0.4 3.6]);
ylabel('z(parameter)');
legend_size(pSize);