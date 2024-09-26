function[questionnaires, categ_quests, n_categ, subject_id] = extract_questionnaires(study_nm, subject_id, NS)
% [questionnaires, categ_quests, n_categ, subject_id] = extract_questionnaires(study_nm, subject_id, NS)
% extract_questionnaires will extract all the questionnaires used in the
% study classified by general categories (stress/anxiety, motivation,
% early-life stress (CTQ), personality (hexaco), dominance, etc.).
%
% INPUTS
% study_nm: study name
%
% subject_id: list of subjects
%
% NS: number of subjects
%
% OUTPUTS
% questionnaires: structure with one subfield for each questionnaire
% category
%
% categ_quests: name of all questionnaire categories (fields of
% questionnaires variable)
%
% n_categ: number of categories in categ_quests
%
% subject_id: list of subjects used for the extraction
%

%% define inputs manually if not entered
% study name
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
% subject list
if ~exist('subject_id','var') || isempty(subject_id)
    condition = subject_condition;
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
end
% NS
if exist('subject_id','var') && (~exist('NS','var') || isempty(NS))
    NS = length(subject_id);
end

%% extract questionnaires
[excelReadQuestionnairesFile, sub_CID_list] = load_questionnaires_data;
[excelReadQuestionnairesFile2] = load_gal_data_bis(study_nm);

%% prepare variables of interest
[n_covid, ISCE,...
    money, age, sex, weight, height, BMI,...
    STAI_T, SIAS, PSS14,...
    CTQ_emotionalA, CTQ_physicalA,...
    CTQ_sexA, CTQ_minDenial,...
    CTQ_emotionalN, CTQ_physicalN,...
    MADRS_S, BIS_NPI, BIS_MI, BIS_AI,...
    Lars_e_AI_EverydayProd, Lars_e_AI_Init,...
    Lars_e_IC_Interest, Lars_e_IC_Novelty, Lars_e_IC_Motiv, Lars_e_IC_Social,...
    Lars_e_ActionInit, Lars_e_IntellectCuriosity, Lars_e_EmotResp, Lars_e_SelfAwareness,...
    MPSTEFS_physical_total, MPSTEFS_mental_total,...
    MPSTEFS_physical_energy, MPSTEFS_mental_energy,...
    MPSTEFS_physical_fatigue, MPSTEFS_mental_fatigue,...
    JPIR, SHAP, SPSRQ_R, SPSRQ_P,...
    IPAQ, IPAQinactivity,...
    PRF_D, CI_enjCompet, CI_contentiousness,...
    socialLadder,...
    hexaco1_HonestyHumility, hexaco2_emotion, hexaco3_extraversion,...
    hexaco4_agreeableness, hexaco5_consciousness, hexaco6_openness] = deal(NaN(1,NS));

% extract sleep
[avg_sleep,...
    prevDay_sleep,...
    prevDay_min_avg_sleep] = extract_sleep(study_nm, subject_id, NS);

%% extract data for each subject selected
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = strcmp(sub_CID_list,sub_nm);
    sub_idx2 = strcmp(excelReadQuestionnairesFile2.CID, ['CID',sub_nm]);
    
    % general
    n_covid(iS) = excelReadQuestionnairesFile.NombreD_infectionsAuCOVID(sub_idx);
    age(iS) = excelReadQuestionnairesFile2.Age_yearsOld_(sub_idx2);
    sex(iS) = strcmp(excelReadQuestionnairesFile2.Sexe_femaleF_maleM_(sub_idx2),'F');
    weight(iS) = excelReadQuestionnairesFile2.Weight(sub_idx2);
    height(iS) = excelReadQuestionnairesFile2.Height(sub_idx2).*100; % convert in cm (instead of m)
    BMI(iS) = excelReadQuestionnairesFile2.BMI(sub_idx2);
    
    % socio-economic status (education, money, etc.)
    ISCE(iS) = excelReadQuestionnairesFile.niveauD__ducation_bas_SurISCE_(sub_idx);
    money(iS) = str2double(excelReadQuestionnairesFile.moyenneRevenuBrutAnnuel(sub_idx))./1000; % convert in thousand CHF (not just CHF)
    socialLadder(iS) = excelReadQuestionnairesFile2.SocialComparisonLadder(sub_idx2);
    
    % stress/anxiety questionnaires
    STAI_T(iS) = excelReadQuestionnairesFile.STAITraitScore(sub_idx);
    SIAS(iS) = excelReadQuestionnairesFile.SIASScore(sub_idx);
    PSS14(iS) = excelReadQuestionnairesFile.PSS_14Score(sub_idx);
    
    % CTQ
    CTQ_emotionalA(iS) = excelReadQuestionnairesFile.CTQEmotionalAbuseScore(sub_idx);
    CTQ_physicalA(iS) = excelReadQuestionnairesFile.CTQPhysicalAbuseScore(sub_idx);
    CTQ_sexA(iS) = excelReadQuestionnairesFile.CTQSexualAbuseScore(sub_idx);
    CTQ_minDenial(iS) = excelReadQuestionnairesFile.CTQMinimizationDenialScore(sub_idx);
    CTQ_emotionalN(iS) = excelReadQuestionnairesFile.CTQEmotionalNeglectScore(sub_idx);
    CTQ_physicalN(iS) = excelReadQuestionnairesFile.CTQPhysicalNeglectScore(sub_idx);
    
    % motivation questionnaires
    JPIR(iS) = excelReadQuestionnairesFile.JPI_RScore(sub_idx);
    MPSTEFS_physical_total(iS) = excelReadQuestionnairesFile.MPSTEFSPhysicalTraitScore(sub_idx);
    MPSTEFS_mental_total(iS) = excelReadQuestionnairesFile.MPSTEFSMentalTraitScore(sub_idx);
    MPSTEFS_physical_energy(iS) = excelReadQuestionnairesFile.MPSTEFSEnergyPhysical(sub_idx);
    MPSTEFS_mental_energy(iS) = excelReadQuestionnairesFile.MPSTEFSEnergyMental(sub_idx);
    MPSTEFS_physical_fatigue(iS) = excelReadQuestionnairesFile.MPSTEFSFatiguePhysical(sub_idx);
    MPSTEFS_mental_fatigue(iS) = excelReadQuestionnairesFile.MPSTEFSFatigueMental(sub_idx);
    MADRS_S(iS) = excelReadQuestionnairesFile.MADRS_SCorrected(sub_idx);
    SPSRQ_R(iS) = excelReadQuestionnairesFile.RewardScore(sub_idx);
    SPSRQ_P(iS) = excelReadQuestionnairesFile.PunishmentScore(sub_idx);
    BIS_AI(iS) = excelReadQuestionnairesFile.AttentionalImpulsivenessScore(sub_idx);
    BIS_MI(iS) = excelReadQuestionnairesFile.MotorImpulsivenessScore(sub_idx);
    BIS_NPI(iS) = excelReadQuestionnairesFile.Non_PlanningImpulsivenessScore(sub_idx);
    Lars_e_AI_EverydayProd(iS) = excelReadQuestionnairesFile.AI_EP(sub_idx);
    Lars_e_AI_Init(iS) = excelReadQuestionnairesFile.AI_I(sub_idx);
    Lars_e_ActionInit(iS) = Lars_e_AI_Init(iS) + Lars_e_AI_EverydayProd(iS);
    Lars_e_IC_Interest(iS) = excelReadQuestionnairesFile.IC_I(sub_idx);
    Lars_e_IC_Novelty(iS) = excelReadQuestionnairesFile.IC_N(sub_idx);
    Lars_e_IC_Motiv(iS) = excelReadQuestionnairesFile.IC_M(sub_idx);
    Lars_e_IC_Social(iS) = excelReadQuestionnairesFile.IC_S(sub_idx);
    Lars_e_IntellectCuriosity(iS) = Lars_e_IC_Interest(iS) + Lars_e_IC_Novelty(iS) +...
        Lars_e_IC_Motiv(iS) + Lars_e_IC_Social(iS);
    Lars_e_EmotResp(iS) = excelReadQuestionnairesFile.ER(sub_idx);
    Lars_e_SelfAwareness(iS) = excelReadQuestionnairesFile.SA(sub_idx);
    SHAP(iS) = excelReadQuestionnairesFile.SHAPScore(sub_idx);
    IPAQinactivity(iS) = excelReadQuestionnairesFile.IPAQInactivity(sub_idx);
    
    % dominance/competition
    PRF_D(iS) = excelReadQuestionnairesFile.PRF_DScore(sub_idx);
    CI_contentiousness(iS) = excelReadQuestionnairesFile.Contentiousness(sub_idx);
    CI_enjCompet(iS) = excelReadQuestionnairesFile.EnjoymentOfCompetition(sub_idx);
    
    % Hexaco
    hexaco1_HonestyHumility(iS) = excelReadQuestionnairesFile.Honesty_Humility(sub_idx);
    hexaco2_emotion(iS) = excelReadQuestionnairesFile.Emotionality(sub_idx);
    hexaco3_extraversion(iS) = excelReadQuestionnairesFile.Extraversion(sub_idx);
    hexaco4_agreeableness(iS) = excelReadQuestionnairesFile.Agreeableness(sub_idx);
    hexaco5_consciousness(iS) = excelReadQuestionnairesFile.Conscientiousness(sub_idx);
    hexaco6_openness(iS) = excelReadQuestionnairesFile.OpennessToExperience(sub_idx);
end % subject loop

%% load clean IPAQ values
IPAQ(:) = IPAQ_rescoring(study_nm, subject_id, NS);

%% regroup questionnaires by category
% general
questionnaires.general.n_covid = n_covid;
questionnaires.general.age = age;
questionnaires.general.sex = sex;
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
questionnaires.sleep.prevDay_min_avg_sleep = prevDay_min_avg_sleep;

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
questionnaires.Motivation.MPSTEFS_physical = MPSTEFS_physical_total;
questionnaires.Motivation.MPSTEFS_mental = MPSTEFS_mental_total;
questionnaires.Motivation.MPSTEFS_physical_energy = MPSTEFS_physical_energy;
questionnaires.Motivation.MPSTEFS_mental_energy = MPSTEFS_mental_energy;
questionnaires.Motivation.MPSTEFS_physical_fatigue = MPSTEFS_physical_fatigue;
questionnaires.Motivation.MPSTEFS_mental_fatigue = MPSTEFS_mental_fatigue;
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
questionnaires.Motivation.IPAQ = IPAQ;
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

end % function