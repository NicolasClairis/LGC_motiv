% function[]= questionnaires_f_sex()

%% subject selection
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

%% extract questionnaires
% prepare data
[STAI_T.male, SIAS.male, PSS14.male,...
    CTQ_emotionalA.male,CTQ_physicalA.male,...
    CTQ_sexA.male,CTQ_minDenial.male,...
    CTQ_emotionalN.male,CTQ_physicalN.male,...
    MADRS_S.male, Lars_e.male,...
    MPSTEFS_physical.male, MPSTEFS_mental.male,...
    JPIR.male, SHAP.male, SPSRQ_R.male, SPSRQ_P.male,...
    PRF_D.male, CI.male,...
    n_children.male] = deal(NaN(1,male_NS));
[STAI_T.female, SIAS.female, PSS14.female,...
    CTQ_emotionalA.female,CTQ_physicalA.female,...
    CTQ_sexA.female,CTQ_minDenial.female,...
    CTQ_emotionalN.female,CTQ_physicalN.female,...
    MADRS_S.female, Lars_e.female,...
    MPSTEFS_physical.female, MPSTEFS_mental.female,...
    JPIR.female, SHAP.female, SPSRQ_R.female, SPSRQ_P.female,...
    PRF_D.female, CI.female,...
    n_children.female] = deal(NaN(1,female_NS));
% load all questionnaires
[excelReadQuestionnairesFile, sub_CID_list] = load_questionnaires_data;

%% extract relevant data
for iSex = 1:2
    switch iSex
        case 1
            sex_nm = 'male';
            subject_id = male_CIDS;
            NS = male_NS;
        case 2
            sex_nm = 'female';
            subject_id = female_CIDS;
            NS = female_NS;
    end
    
    %% loop over subjects
    for iS = 1:NS
        sub_nm = subject_id{iS};
        sub_idx = strcmp(sub_CID_list,sub_nm);
        
        % stress/anxiety questionnaires
        STAI_T.(sex_nm)(iS) = excelReadQuestionnairesFile.STAITraitScore(sub_idx);
        SIAS.(sex_nm)(iS) = excelReadQuestionnairesFile.SIASScore(sub_idx);
        PSS14.(sex_nm)(iS) = excelReadQuestionnairesFile.PSS_14Score(sub_idx);
        
        % motivation questionnaires
        
        % dominance
        
        % CTQ
        
        % Hexaco
    end % subject loop
end % loop over sex

%% perform test for each
% perform comparison on raw data first
[~,pval.uncorrected.STAI_T] = ttest2(STAI_T.male, STAI_T.female);
[~,pval.uncorrected.SIAS] = ttest2(SIAS.male, SIAS.female);
[~,pval.uncorrected.PSS14] = ttest2(PSS14.male, PSS14.female);

% use pval_adjust if you want to correct for multiple comparisons

% end % function