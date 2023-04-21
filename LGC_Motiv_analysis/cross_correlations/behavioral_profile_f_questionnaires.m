% split subjects according to MADRS-S and JPI-R levels and look at the
% profiles accordingly

%% working directories
computerRoot = LGCM_root_paths;
%% define subjects
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load questionnaire scores
[excelReadQuestionnairesFile, quest_S_sub_CID_list] = load_questionnaires_data();
MADRS_S_fullList = excelReadQuestionnairesFile.MADRS_SCorrected;
JPI_R_fullList = excelReadQuestionnairesFile.JPI_RScore;
[MADRS_S_score, JPI_R_score] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    MADRS_S_score(iS) = MADRS_S_fullList(strcmp(quest_S_sub_CID_list, sub_nm));
    JPI_R_score(iS) = JPI_R_fullList(strcmp(quest_S_sub_CID_list, sub_nm));
end % subject list

%% define groups based on scores
MADRS_S.low = MADRS_S_score < 4;
MADRS_S.mid = (MADRS_S_score >= 4).*(MADRS_S_score < 6) == 1;
MADRS_S.high = MADRS_S_score >= 6;

JPI_R.low = JPI_R_score < 10;
JPI_R.high = JPI_R_score >= 10;

%% load parameters
[mdlType, mdlN] = behavioral_model_selection;
behavioralPrm = prm_extraction(study_nm, subject_id,mdlType, mdlN);
behavioral_prm_names = fieldnames(behavioralPrm);
behavioral_prm_names(strcmp(behavioral_prm_names,'CID')) = [];
n_prm = length(behavioral_prm_names);

%% extract parameters for each group
for iPrm = 1:n_prm
    prm_nm = behavioral_prm_names{iPrm};
    % MADRS-S
    [m_prm.(prm_nm).MADRS_S.low,...
        sem_prm.(prm_nm).MADRS_S.low] = mean_sem_sd(behavioralPrm.(prm_nm)(MADRS_S.low),2);
    [m_prm.(prm_nm).MADRS_S.mid,...
        sem_prm.(prm_nm).MADRS_S.mid] = mean_sem_sd(behavioralPrm.(prm_nm)(MADRS_S.mid),2);
    [m_prm.(prm_nm).MADRS_S.high,...
        sem_prm.(prm_nm).MADRS_S.high] = mean_sem_sd(behavioralPrm.(prm_nm)(MADRS_S.high),2);
    
    % JPI-R
    [m_prm.(prm_nm).JPI_R.low,...
        sem_prm.(prm_nm).JPI_R.low] = mean_sem_sd(behavioralPrm.(prm_nm)(JPI_R.low),2);
    [m_prm.(prm_nm).JPI_R.high,...
        sem_prm.(prm_nm).JPI_R.high] = mean_sem_sd(behavioralPrm.(prm_nm)(JPI_R.high),2);
end % parameter loop

%% display figure
lWidth = 3;
pSize = 30;
%% MADRS-S
fig;
hold on;
x_MADRS_S = 1:3; % 3 groups
xPrm = (-0.4):(1/(n_prm+1)):0.4;
bar_hdl_leg = NaN(1,n_prm);
for iPrm = 1:n_prm
    prm_nm = behavioral_prm_names{iPrm};
    bar_hdl = bar([x_MADRS_S+xPrm(iPrm)],...
        [m_prm.(prm_nm).MADRS_S.low,...
        m_prm.(prm_nm).MADRS_S.mid,...
        m_prm.(prm_nm).MADRS_S.high]);
    bar_hdl.BarWidth = (1/(n_prm+1));
    bar_hdl_leg(iPrm) = bar_hdl;
    er_hdl = errorbar([x_MADRS_S+xPrm(iPrm)],...
        [m_prm.(prm_nm).MADRS_S.low,...
        m_prm.(prm_nm).MADRS_S.mid,...
        m_prm.(prm_nm).MADRS_S.high],...
        [sem_prm.(prm_nm).MADRS_S.low,...
        sem_prm.(prm_nm).MADRS_S.mid,...
        sem_prm.(prm_nm).MADRS_S.high]);
    er_hdl.LineStyle = 'none';
    er_hdl.LineWidth = lWidth;
    er_hdl.Color = 'k';
end
legend(bar_hdl_leg, behavioral_prm_names);
legend('boxoff');
xticks(x_MADRS_S);
xticklabels({'MADRS-S low','MADRS-S mid','MADRS-S high'});
xlim([0.4 3.6]);
ylabel('Parameter');
legend_size(pSize);

%% JPI-R
fig;
hold on;
x_JPI_R = 1:2; % 2 groups
xPrm = (-0.4):(1/(n_prm+1)):0.4;
bar_hdl_leg = NaN(1,n_prm);
for iPrm = 1:n_prm
    prm_nm = behavioral_prm_names{iPrm};
    bar_hdl = bar([x_JPI_R+xPrm(iPrm)],...
        [m_prm.(prm_nm).JPI_R.low,...
        m_prm.(prm_nm).JPI_R.high]);
    bar_hdl.BarWidth = (1/(n_prm+1));
    bar_hdl_leg(iPrm) = bar_hdl;
    er_hdl = errorbar([x_JPI_R+xPrm(iPrm)],...
        [m_prm.(prm_nm).JPI_R.low,...
        m_prm.(prm_nm).JPI_R.high],...
        [sem_prm.(prm_nm).JPI_R.low,...
        sem_prm.(prm_nm).JPI_R.high]);
    er_hdl.LineStyle = 'none';
    er_hdl.LineWidth = lWidth;
    er_hdl.Color = 'k';
end
legend(bar_hdl_leg, behavioral_prm_names);
legend('boxoff');
xticks(x_JPI_R);
xticklabels({'JPI-R low','JPI-R high'});
xlim([0.4 2.6]);
ylabel('Parameter');
legend_size(pSize);