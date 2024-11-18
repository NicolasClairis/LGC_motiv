% function[]= questionnaires_f_sex()
% questionnaires_f_sex will compare the different questionnaires used in
% the study between male and female, after grouping them by category.


%% subject selection
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

%% extract questionnaires
[excelReadQuestionnairesFile, sub_CID_list] = load_questionnaires_data;
[excelReadQuestionnairesFile2] = load_gal_data_bis(study_nm);

%% filter weird IPAQ values
IPAQ_all = excelReadQuestionnairesFile.IPAQ;
[~,~,cleaned_IPAQ] = rmv_outliers_3sd(IPAQ_all);

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
    
    %% prepare variables of interest
    [n_covid.(sex_nm), ISCE.(sex_nm),...
        money.(sex_nm), age.(sex_nm), weight.(sex_nm), height.(sex_nm), BMI.(sex_nm),...
        avg_sleep.(sex_nm), prevDay_sleep.(sex_nm), avg_min_prevDay_sleep.(sex_nm),...
        STAI_T.(sex_nm), SIAS.(sex_nm), PSS14.(sex_nm),...
        CTQ_emotionalA.(sex_nm),CTQ_physicalA.(sex_nm),...
        CTQ_sexA.(sex_nm),CTQ_minDenial.(sex_nm),...
        CTQ_emotionalN.(sex_nm),CTQ_physicalN.(sex_nm),...
        MADRS_S.(sex_nm), BIS_NPI.(sex_nm), BIS_MI.(sex_nm), BIS_AI.(sex_nm),...
        Lars_e_AI_EverydayProd.(sex_nm), Lars_e_AI_Init.(sex_nm),...
        Lars_e_IC_Interest.(sex_nm), Lars_e_IC_Novelty.(sex_nm), Lars_e_IC_Motiv.(sex_nm), Lars_e_IC_Social.(sex_nm),...
        Lars_e_ActionInit.(sex_nm), Lars_e_IntellectCuriosity.(sex_nm), Lars_e_EmotResp.(sex_nm), Lars_e_SelfAwareness.(sex_nm),...
        MPSTEFS_physical.(sex_nm), MPSTEFS_mental.(sex_nm),...
        JPIR.(sex_nm), SHAP.(sex_nm), SPSRQ_R.(sex_nm), SPSRQ_P.(sex_nm),...
        IPAQ.(sex_nm), IPAQinactivity.(sex_nm),...
        PRF_D.(sex_nm), CI_enjCompet.(sex_nm), CI_contentiousness.(sex_nm),...
        socialLadder.(sex_nm),...
        hexaco1_HonestyHumility.(sex_nm), hexaco2_emotion.(sex_nm), hexaco3_extraversion.(sex_nm),...
        hexaco4_agreeableness.(sex_nm), hexaco5_consciousness.(sex_nm), hexaco6_openness.(sex_nm)] = deal(NaN(1,NS));
    
    % extract sleep
    [avg_sleep.(sex_nm),...
        prevDay_sleep.(sex_nm),...
        avg_min_prevDay_sleep.(sex_nm)] = extract_sleep(study_nm, subject_id, NS);
    
    %% loop over subjects
    for iS = 1:NS
        sub_nm = subject_id{iS};
        sub_idx = strcmp(sub_CID_list,sub_nm);
        sub_idx2 = strcmp(excelReadQuestionnairesFile2.CID, ['CID',sub_nm]);
        
        % general
        n_covid.(sex_nm)(iS) = excelReadQuestionnairesFile.NombreD_infectionsAuCOVID(sub_idx);
        age.(sex_nm)(iS) = excelReadQuestionnairesFile2.Age_yearsOld_(sub_idx2);
        weight.(sex_nm)(iS) = excelReadQuestionnairesFile2.Weight(sub_idx2);
        height.(sex_nm)(iS) = excelReadQuestionnairesFile2.Height(sub_idx2).*100; % convert in cm (instead of m)
        BMI.(sex_nm)(iS) = excelReadQuestionnairesFile2.BMI(sub_idx2);
        
        % socio-economic status (education, money, etc.)
        ISCE.(sex_nm)(iS) = excelReadQuestionnairesFile.niveauD__ducation_bas_SurISCE_(sub_idx);
        money.(sex_nm)(iS) = str2double(excelReadQuestionnairesFile.moyenneRevenuBrutAnnuel(sub_idx))./1000; % convert in thousand CHF (not just CHF)
        socialLadder.(sex_nm)(iS) = excelReadQuestionnairesFile2.SocialComparisonLadder(sub_idx2);
        
        % stress/anxiety questionnaires
        STAI_T.(sex_nm)(iS) = excelReadQuestionnairesFile.STAITraitScore(sub_idx);
        SIAS.(sex_nm)(iS) = excelReadQuestionnairesFile.SIASScore(sub_idx);
        PSS14.(sex_nm)(iS) = excelReadQuestionnairesFile.PSS_14Score(sub_idx);
        
        % CTQ
        CTQ_emotionalA.(sex_nm)(iS) = excelReadQuestionnairesFile.CTQEmotionalAbuseScore(sub_idx);
        CTQ_physicalA.(sex_nm)(iS) = excelReadQuestionnairesFile.CTQPhysicalAbuseScore(sub_idx);
        CTQ_sexA.(sex_nm)(iS) = excelReadQuestionnairesFile.CTQSexualAbuseScore(sub_idx);
        CTQ_minDenial.(sex_nm)(iS) = excelReadQuestionnairesFile.CTQMinimizationDenialScore(sub_idx);
        CTQ_emotionalN.(sex_nm)(iS) = excelReadQuestionnairesFile.CTQEmotionalNeglectScore(sub_idx);
        CTQ_physicalN.(sex_nm)(iS) = excelReadQuestionnairesFile.CTQPhysicalNeglectScore(sub_idx);
        
        % motivation questionnaires
        JPIR.(sex_nm)(iS) = excelReadQuestionnairesFile.JPI_RScore(sub_idx);
        MPSTEFS_physical.(sex_nm)(iS) = excelReadQuestionnairesFile.MPSTEFSPhysicalTraitScore(sub_idx);
        MPSTEFS_mental.(sex_nm)(iS) = excelReadQuestionnairesFile.MPSTEFSMentalTraitScore(sub_idx);
        MADRS_S.(sex_nm)(iS) = excelReadQuestionnairesFile.MADRS_SCorrected(sub_idx);
        SPSRQ_R.(sex_nm)(iS) = excelReadQuestionnairesFile.RewardScore(sub_idx);
        SPSRQ_P.(sex_nm)(iS) = excelReadQuestionnairesFile.PunishmentScore(sub_idx);
        BIS_AI.(sex_nm)(iS) = excelReadQuestionnairesFile.AttentionalImpulsivenessScore(sub_idx);
        BIS_MI.(sex_nm)(iS) = excelReadQuestionnairesFile.MotorImpulsivenessScore(sub_idx);
        BIS_NPI.(sex_nm)(iS) = excelReadQuestionnairesFile.Non_PlanningImpulsivenessScore(sub_idx);
        Lars_e_AI_EverydayProd.(sex_nm)(iS) = excelReadQuestionnairesFile.AI_EP(sub_idx);
        Lars_e_AI_Init.(sex_nm)(iS) = excelReadQuestionnairesFile.AI_I(sub_idx);
        Lars_e_ActionInit.(sex_nm)(iS) = Lars_e_AI_Init.(sex_nm)(iS) + Lars_e_AI_EverydayProd.(sex_nm)(iS);
        Lars_e_IC_Interest.(sex_nm)(iS) = excelReadQuestionnairesFile.IC_I(sub_idx);
        Lars_e_IC_Novelty.(sex_nm)(iS) = excelReadQuestionnairesFile.IC_N(sub_idx);
        Lars_e_IC_Motiv.(sex_nm)(iS) = excelReadQuestionnairesFile.IC_M(sub_idx);
        Lars_e_IC_Social.(sex_nm)(iS) = excelReadQuestionnairesFile.IC_S(sub_idx);
        Lars_e_IntellectCuriosity.(sex_nm)(iS) = Lars_e_IC_Interest.(sex_nm)(iS) + Lars_e_IC_Novelty.(sex_nm)(iS) +...
            Lars_e_IC_Motiv.(sex_nm)(iS) + Lars_e_IC_Social.(sex_nm)(iS);
        Lars_e_EmotResp.(sex_nm)(iS) = excelReadQuestionnairesFile.ER(sub_idx);
        Lars_e_SelfAwareness.(sex_nm)(iS) = excelReadQuestionnairesFile.SA(sub_idx);
        SHAP.(sex_nm)(iS) = excelReadQuestionnairesFile.SHAPScore(sub_idx);
        IPAQ.(sex_nm)(iS) = cleaned_IPAQ(sub_idx);
        IPAQinactivity.(sex_nm)(iS) = excelReadQuestionnairesFile.IPAQInactivity(sub_idx);
        
        % dominance/competition
        PRF_D.(sex_nm)(iS) = excelReadQuestionnairesFile.PRF_DScore(sub_idx);
        CI_contentiousness.(sex_nm)(iS) = excelReadQuestionnairesFile.Contentiousness(sub_idx);
        CI_enjCompet.(sex_nm)(iS) = excelReadQuestionnairesFile.EnjoymentOfCompetition(sub_idx);
        
        % Hexaco
        hexaco1_HonestyHumility.(sex_nm)(iS) = excelReadQuestionnairesFile.Honesty_Humility(sub_idx);
        hexaco2_emotion.(sex_nm)(iS) = excelReadQuestionnairesFile.Emotionality(sub_idx);
        hexaco3_extraversion.(sex_nm)(iS) = excelReadQuestionnairesFile.Extraversion(sub_idx);
        hexaco4_agreeableness.(sex_nm)(iS) = excelReadQuestionnairesFile.Agreeableness(sub_idx);
        hexaco5_consciousness.(sex_nm)(iS) = excelReadQuestionnairesFile.Conscientiousness(sub_idx);
        hexaco6_openness.(sex_nm)(iS) = excelReadQuestionnairesFile.OpennessToExperience(sub_idx);
    end % subject loop
end % loop over sex

%% regroup questionnaires by category
% general
questionnaires.general.n_covid = n_covid;
questionnaires.general.age = age;
questionnaires.general.weight = weight;
questionnaires.general.height = height;
questionnaires.general.BMI = BMI;

% socio-economic status (education, money, etc.)
questionnaires.SES.ISCE = ISCE;
questionnaires.SES.ladder = socialLadder;
questionnaires.SES.money = money;

% sleep
questionnaires.sleep.avg = avg_sleep;
questionnaires.sleep.prevDay = prevDay_sleep;
questionnaires.sleep.avg_min_prevDay = avg_min_prevDay_sleep;

% stress/anxiety
questionnaires.stress_anxiety.STAI_T = STAI_T;
questionnaires.stress_anxiety.SIAS = SIAS;
questionnaires.stress_anxiety.PSS14 = PSS14;

% CTQ
questionnaires.CTQ.emotionalAbuse = CTQ_emotionalA;
questionnaires.CTQ.physicalAbuse = CTQ_physicalA;
questionnaires.CTQ.sexualAbuse = CTQ_sexA;
questionnaires.CTQ.minimizationDenial = CTQ_minDenial;
questionnaires.CTQ.emotionalNeglect = CTQ_emotionalN;
questionnaires.CTQ.physicalNeglect = CTQ_physicalN;

% motivation
questionnaires.Motivation.JPIR = JPIR;
questionnaires.Motivation.MPSTEFS_physical = MPSTEFS_physical;
questionnaires.Motivation.MPSTEFS_mental = MPSTEFS_mental;
questionnaires.Motivation.MADRS_S = MADRS_S;
questionnaires.Motivation.SPSRQ_R = SPSRQ_R;
questionnaires.Motivation.SPSRQ_P = SPSRQ_P;
questionnaires.Motivation.BIS_AI = BIS_AI;
questionnaires.Motivation.BIS_MI = BIS_MI;
questionnaires.Motivation.BIS_NPI = BIS_NPI;
questionnaires.Motivation.Lars_e_AI_EverydayProd = Lars_e_AI_EverydayProd;
questionnaires.Motivation.Lars_e_AI_Init = Lars_e_AI_Init;
questionnaires.Motivation.Lars_e_ActionInit = Lars_e_ActionInit;
questionnaires.Motivation.Lars_e_IC_Interest = Lars_e_IC_Interest;
questionnaires.Motivation.Lars_e_IC_Novelty = Lars_e_IC_Novelty;
questionnaires.Motivation.Lars_e_IC_Motiv = Lars_e_IC_Motiv;
questionnaires.Motivation.Lars_e_IC_Social = Lars_e_IC_Social;
questionnaires.Motivation.Lars_e_IntellectCuriosity = Lars_e_IntellectCuriosity;
questionnaires.Motivation.Lars_e_EmotResp = Lars_e_EmotResp;
questionnaires.Motivation.Lars_e_SelfAwareness = Lars_e_SelfAwareness;
questionnaires.Motivation.SHAP = SHAP;
% questionnaires.Motivation.IPAQ = IPAQ; % remove IPAQ as too noisy (+ too
% high range compared to others)
questionnaires.Motivation.IPAQinactivity = IPAQinactivity;

% dominance/competition
questionnaires.dominance_compet.PRF_D = PRF_D;
questionnaires.dominance_compet.CI_contentiousness = CI_contentiousness;
questionnaires.dominance_compet.CI_enjCompet = CI_enjCompet;

% hexaco
questionnaires.hexaco.honestyHumility = hexaco1_HonestyHumility;
questionnaires.hexaco.emotion = hexaco2_emotion;
questionnaires.hexaco.extraversion = hexaco3_extraversion;
questionnaires.hexaco.agreeableness = hexaco4_agreeableness;
questionnaires.hexaco.consciousness = hexaco5_consciousness;
questionnaires.hexaco.openness = hexaco6_openness;

%% extract number of categories + name of all questionnaires in each category
categ_quests = fieldnames(questionnaires);
n_categ = length(categ_quests);

%% method for multiple comparisons correction
corr_method = 'bonferroni';
corr_method_nm = ['corr_',corr_method];

%% colour code
[pSize, lW, col, mSize] = general_fig_prm;
female_col = col.red;
male_col = col.blue_dark;

%% perform test for each questionnaire
for iCateg = 1:n_categ
    categ_nm = categ_quests{iCateg};
    
    % category subfields
    quest_names = fieldnames(questionnaires.(categ_nm));
    n_quests = length(quest_names);
    
    pval_unc = [];
    for iQ = 1:n_quests
        quest_nm = quest_names{iQ};
        
        % perform comparison on raw data first
        [~,pval.uncorrected.(categ_nm).(quest_nm)] = ttest2(...
            questionnaires.(categ_nm).(quest_nm).male,...
            questionnaires.(categ_nm).(quest_nm).female);
        pval_unc = [pval_unc, pval.uncorrected.(categ_nm).(quest_nm)];
        
        % store significant results in a different subfield
        if pval.uncorrected.(categ_nm).(quest_nm) < 0.05
            % store p.value
            signif.unc.(categ_nm).(quest_nm).pval = pval.uncorrected.(categ_nm).(quest_nm);
            % store mean and SEM males
            [signif.unc.(categ_nm).(quest_nm).male_m,...
                signif.unc.(categ_nm).(quest_nm).male_sem] = mean_sem_sd(questionnaires.(categ_nm).(quest_nm).male,2);
            % store mean and SEM females
            [signif.unc.(categ_nm).(quest_nm).female_m,...
                signif.unc.(categ_nm).(quest_nm).female_sem] = mean_sem_sd(questionnaires.(categ_nm).(quest_nm).female,2);
        end
    end % loop over questionnaires
    
    %% attempt to correct for multiple comparisons
    % use pval_adjust if you want to correct for multiple comparisons
    % temporarily decided to correct for the number of dimensions
    % tested (considering that the tests should not be considered to be independent within the same dimension)
    % but note this is quite arbitrary... Another way would be to
    % correct for the number of questionnaires or for the number of
    % scores extracted, but again the independence between tests is
    % highly debatable...
    [pval_corr] = pval_adjust(pval_unc,corr_method);
    % extract individual p.values
    for iQ = 1:n_quests
        quest_nm = quest_names{iQ};
        pval.(corr_method_nm).(categ_nm).(quest_nm) = pval_corr(iQ);
        
        % store significant results in a different subfield
        if pval.(corr_method_nm).(categ_nm).(quest_nm) < 0.05
            % store p.value
            signif.(corr_method_nm).(categ_nm).(quest_nm).pval = pval.(corr_method_nm).(categ_nm).(quest_nm);
            % store mean and SEM males
            [signif.(corr_method_nm).(categ_nm).(quest_nm).male_m,...
                signif.(corr_method_nm).(categ_nm).(quest_nm).male_sem] = mean_sem_sd(questionnaires.(categ_nm).(quest_nm).male,2);
            % store mean and SEM females
            [signif.(corr_method_nm).(categ_nm).(quest_nm).female_m,...
                signif.(corr_method_nm).(categ_nm).(quest_nm).female_sem] = mean_sem_sd(questionnaires.(categ_nm).(quest_nm).female,2);
        end
    end % loop over questionnaires
    
    %% show results for current in a graph (+ if results are significant or not)
    fig;
    for iQ = 1:n_quests
        quest_nm = quest_names{iQ};
        jPos_male = 1 + 2*(iQ - 1);
        jPos_female = 2 + 2*(iQ - 1);
        
        if strcmp(categ_nm,'SES')
            switch quest_nm
                case 'money'
                    yyaxis right;
            end
        end
        
        % show male vs female data
        ok_males = ~isnan(questionnaires.(categ_nm).(quest_nm).male);
        male_violin = Violin({questionnaires.(categ_nm).(quest_nm).male(ok_males)},jPos_male,...
            'ViolinColor',{male_col});
        ok_females = ~isnan(questionnaires.(categ_nm).(quest_nm).female);
        female_violin = Violin({questionnaires.(categ_nm).(quest_nm).female(ok_females)},jPos_female,...
            'ViolinColor',{female_col});
        
        % add p.value indication if difference is significant        
        [l_hdl, star_hdl] = add_pval_comparison(questionnaires.(categ_nm).(quest_nm).male,...
            questionnaires.(categ_nm).(quest_nm).female,...
            pval.(corr_method_nm).(categ_nm).(quest_nm), jPos_male, jPos_female, '');
    end % questionnaire loop
    switch categ_nm
        case 'sleep'
            ylabel('Time (min)');
        case 'SES'
            switch quest_nm
                case 'money'
                    yyaxis right;
                    ylabel('Money (kCHF/year)');
                otherwise
                    ylabel('Scores');
            end
        otherwise
            ylabel('Scores');
    end
    xticks(1.5:2:n_quests*2);
    xticklabels(quest_names);
    switch categ_nm
        case 'dominance_compet'
            title('Dominance/Competitiveness');
        case 'stress_anxiety'
            title('Stress/Anxiety');
        otherwise
            title(categ_nm);
    end
end % loop over questionnaire categories
% end % function