% check correlation between stress through questionnaires( STAI-T, PSS-14),
% VAS scale the day of the experiment and hormones (Cortisol,
% Testosterone), gender and model parameters

%% load subject
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);
% load gender
males_id = LGCM_subject_selection(study_nm,condition,'males');
females_id = LGCM_subject_selection(study_nm,condition,'females');

%% load stress-related data
[excelReadQuestionnairesFile, quest_sub_CID_list] = load_questionnaires_data;

%% load hormonal data
[TESTO_data] = load_TESTO(study_nm);
[CORT_data] = load_CORT(study_nm);

[excelReadGeneralFile] = load_gal_data_bis(study_nm);

%% load model parameters
[prm] = prm_extraction(study_nm, subject_id);
prm_sub_CID_list = prm.CID;

%% extract the data for the selected subjects
[STAI_T_score, PSS14_score,...
    CORT_AUCg, TESTO_AUCg,...
    kEp, kEm, kFp, kFm] = deal(NaN(NS,1));
for iS = 1:NS
    sub_nm = subject_id{iS};
    CORT_sub_idx = strcmp(['CID',sub_nm], CORT_data.CID);
    TESTO_sub_idx = strcmp(['CID',sub_nm], TESTO_data.CID);
    quest_sub_idx = strcmp(sub_nm, quest_sub_CID_list);
    prm_sub_idx = strcmp(sub_nm, prm_sub_CID_list);

    CORT_AUCg(iS) = CORT_data.AUCg(CORT_sub_idx);
    TESTO_AUCg(iS) = TESTO_data.AUCg(TESTO_sub_idx);
    STAI_T_score(iS) = excelReadQuestionnairesFile.STAITraitScore(quest_sub_idx);
    PSS14_score(iS) = excelReadQuestionnairesFile.PSS_14Score(quest_sub_idx);
    kEp(iS) = prm.kEp(prm_sub_idx);
    kEm(iS) = prm.kEm(prm_sub_idx);
    kFp(iS) = prm.kFp(prm_sub_idx);
    kFm(iS) = prm.kLm(prm_sub_idx);
end % subject loop

%% check stress vs model parameters
% kEp
[rho.kEp_f_STAI, pval.kEp_f_STAI] = corr(kEp, STAI_T_score);
[rho.kEp_f_PSS14, pval.kEp_f_PSS14] = corr(kEp, PSS14_score);
% [rho.kEp_f_stress_rtg, pval.kEp_f_stress_rtg] = corr(kEp, );
% kEm
[rho.kEm_f_STAI, pval.kEm_f_STAI] = corr(kEm, STAI_T_score);
[rho.kEm_f_PSS14, pval.kEm_f_PSS14] = corr(kEm, PSS14_score);
% [rho.kEm_f_stress_rtg, pval.kEm_f_stress_rtg] = corr(kEm, );
% kFp
[rho.kFp_f_STAI, pval.kFp_f_STAI] = corr(kFp, STAI_T_score);
[rho.kFp_f_PSS14, pval.kFp_f_PSS14] = corr(kFp, PSS14_score);
% [rho.kFp_f_stress_rtg, pval.kFp_f_stress_rtg] = corr(kFp, );
% kFm
[rho.kFm_f_STAI, pval.kFm_f_STAI] = corr(kFm, STAI_T_score);
[rho.kFm_f_PSS14, pval.kFm_f_PSS14] = corr(kFm, PSS14_score);
% [rho.kFm_f_stress_rtg, pval.kFm_f_stress_rtg] = corr(kFm, );