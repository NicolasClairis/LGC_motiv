% script to create a mega-pool of all the information available for all
% subjects and create a huge table accordingly.
%
% Will also perform a circos plot to show correlations between different
% categories (requires circularGraph toolbox https://github.com/paul-kassebaum-mathworks/circularGraph
% for this graph).

%% extract list of all subjects
study_nm = 'study1';
condition = 'fullList';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% define all variables of interest (except 138 FFQ variables)
[general.age, general.sex, general.ISCE,...
    general.weight, general.height, general.BMI,...
    general.n_covid, general.money,...
    general.avg_sleep, general.prevDay_sleep, general.prevDay_min_avg_sleep,...
    general.hexaco.honestyHumility, general.hexaco.emotion, general.hexaco.extraversion,...
    general.hexaco.agreeableness, general.hexaco.consciousness, general.hexaco.openness,...
    salivary.TESTO1, salivary.TESTO4,...
    salivary.CORT1, salivary.CORT2, salivary.CORT3, salivary.CORT4,...
    salivary.TESTO_div_CORT1, salivary.TESTO_div_CORT4,...
    salivary.IL1b, salivary.IL6, salivary.IL18,...
    whole_blood.NAM, whole_blood.NMN, whole_blood.NR,...
    whole_blood.NAD, whole_blood.NADH,...
    whole_blood.NADP, whole_blood.NADPH,...
    whole_blood.MeNam, whole_blood.MeXPY,...
    whole_blood.NAD_div_NADH, whole_blood.NADP_div_NADPH,...
    whole_blood.tNAD, whole_blood.tNAD_with_precursors,...
    whole_blood.tNAD_with_byproducts, whole_blood.NAD_byproducts,...
    plasma.aa.Ala, plasma.aa.Arg, plasma.aa.Asp, plasma.aa.Asn, plasma.aa.Glu,...
    plasma.aa.Gln, plasma.aa.Gly, plasma.aa.His, plasma.aa.Ile, plasma.aa.Leu,...
    plasma.aa.Lys, plasma.aa.Met, plasma.aa.Phe, plasma.aa.Pro, plasma.aa.Ser,...
    plasma.aa.Thr, plasma.aa.Tyr, plasma.aa.Val, plasma.aa.Tau,...
    plasma.aa.trans4hydroxyproline, plasma.aa.x3MethylHis, plasma.aa.x1MethylHis,...
    plasma.Lac, plasma.fa.AcetoAcetate, plasma.fa.BHB, plasma.fa.Acetate,...
    plasma.fa.Propionate, plasma.fa.IsoButyrate, plasma.fa.Butyrate,...
    plasma.fa.x2MethylButyrate, plasma.fa.isovalerate, plasma.fa.valerate,...
    plasma.fa.hexanoicA, plasma.fa.octanoicA, plasma.fa.decanoicA,...
    plasma.fa.AABA, plasma.fa.L_Citrulline, plasma.fa.SDMA, plasma.fa.ADMA,...
    plasma.fa.sarcosine,...
    brainM.dmPFC.Asp, brainM.dmPFC.GABA, brainM.dmPFC.Gln, brainM.dmPFC.Glu,...
    brainM.dmPFC.Glx, brainM.dmPFC.Gln_div_Glu, brainM.dmPFC.GSH, brainM.dmPFC.Gly,...
    brainM.dmPFC.Ins, brainM.dmPFC.Lac, brainM.dmPFC.Tau, brainM.dmPFC.NAA,...
    brainM.dmPFC.NAAG, brainM.dmPFC.GPC_PCho, brainM.dmPFC.Cr_PCr,...
    brainM.aIns.Asp, brainM.aIns.GABA, brainM.aIns.Gln, brainM.aIns.Glu,...
    brainM.aIns.Glx, brainM.aIns.Gln_div_Glu, brainM.aIns.GSH, brainM.aIns.Gly,...
    brainM.aIns.Ins, brainM.aIns.Lac, brainM.aIns.Tau, brainM.aIns.NAA,...
    brainM.aIns.NAAG, brainM.aIns.GPC_PCho, brainM.aIns.Cr_PCr,...
    behavior.questionnaires.stress_anxiety.PSS14,...
    behavior.questionnaires.stress_anxiety.STAI_T,...
    behavior.questionnaires.stress_anxiety.SIAS,...
    behavior.questionnaires.stress_anxiety.CTQ_emotionalA,...
    behavior.questionnaires.stress_anxiety.CTQ_physicalA,...
    behavior.questionnaires.stress_anxiety.CTQ_sexA,...
    behavior.questionnaires.stress_anxiety.CTQ_minDenial,...
    behavior.questionnaires.stress_anxiety.CTQ_emotionalN,...
    behavior.questionnaires.stress_anxiety.CTQ_physicalN,...
    behavior.questionnaires.motivation.MADRS_S,...
    behavior.questionnaires.motivation.BIS_NPI,...
    behavior.questionnaires.motivation.BIS_MI,...
    behavior.questionnaires.motivation.BIS_AI,...
    behavior.questionnaires.motivation.Lars_e_AI_EverydayProd,...
    behavior.questionnaires.motivation.Lars_e_AI_Init,...
    behavior.questionnaires.motivation.Lars_e_IC_Interest,...
    behavior.questionnaires.motivation.Lars_e_IC_Novelty,...
    behavior.questionnaires.motivation.Lars_e_IC_Motiv,...
    behavior.questionnaires.motivation.Lars_e_IC_Social,...
    behavior.questionnaires.motivation.Lars_e_ActionInit,...
    behavior.questionnaires.motivation.Lars_e_IntellectCuriosity,...
    behavior.questionnaires.motivation.Lars_e_EmotResp,...
    behavior.questionnaires.motivation.Lars_e_SelfAwareness,...
    behavior.questionnaires.motivation.MPSTEFS_physical,...
    behavior.questionnaires.motivation.MPSTEFS_mental,...
    behavior.questionnaires.motivation.JPIR,...
    behavior.questionnaires.motivation.SHAP,...
    behavior.questionnaires.motivation.SPSRQ_R,...
    behavior.questionnaires.motivation.SPSRQ_P,...
    behavior.questionnaires.motivation.IPAQ,...
    behavior.questionnaires.motivation.IPAQinactivity,...
    behavior.questionnaires.dominance.PRF_D,...
    behavior.questionnaires.dominance.CI_enjCompet,...
    behavior.questionnaires.dominance.CI_contentiousness,...
    behavior.questionnaires.dominance.socialLadder,...
    behavior.stress_ratings.S1, behavior.stress_ratings.S2,...
    behavior.stress_ratings.S3, behavior.stress_ratings.S4,...
    behavior.stress_ratings.S4_min_S1, behavior.stress_ratings.S4_min_S3,...
    behavior.fatigue_ratings.F1, behavior.fatigue_ratings.F2,...
    behavior.fatigue_ratings.F3, behavior.fatigue_ratings.F4,...
    behavior.fatigue_ratings.F4_min_F1, behavior.fatigue_ratings.F4_min_F3,...
    behavior.task.choices_HE, behavior.task.choices_HPE, behavior.task.choices_HME,...
    behavior.task.prm_kR, behavior.task.prm_kP,...
    behavior.task.prm_kEp, behavior.task.prm_kEm,...
    behavior.task.prm_kFp, behavior.task.prm_kLm,...
    behavior.task.prm_kBias] = deal(NaN(1,NS));

%% extract all the relevant data

% salivary hormones + interleukins
[TESTO_data] = load_TESTO(study_nm, subject_id);
[CORT_data] = load_CORT(study_nm, subject_id);
[IL_data] = load_IL(study_nm, subject_id);

% whole-blood NAD-omics
[wholeBlood, sub_List] = load_blood_NAD(study_nm, subject_id);

% plasma metabolites (amino-acids + acids)
[plasmaM, mb_names, n_mb] = load_plasma_metabolites(subject_id);

% brain metabolites extracted with 1H-MRS
brainMetabolites = metabolite_load(subject_id);

% load FFQ
[nutrition, nutri_categories_names, n_nutri_categ, nutri_names, n_nutri] = load_FFQ(study_nm, subject_id, NS);

% load behavior (questionnaires, % choices, behavioral parameters)
[questionnaires, categ_quests, n_quest_categ, subject_id] = extract_questionnaires(study_nm, subject_id, NS);
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, 0);
prm = prm_extraction(study_nm, subject_id, 'bayesian', '5');
% load also stress and fatigue questionnaires on the day of the experiment
[~,...
    preMRS_stress, postMRS_stress,...
    prefMRI_stress, postfMRI_stress,...
    deltaStressPrePostExp, deltaStressPrePostfMRI] = extract_subjective_stress_ratings(study_nm, subject_id, NS);
[~,...
    preMRS_fatigue, postMRS_fatigue,...
    prefMRI_fatigue, postfMRI_fatigue,...
    deltaFatiguePrePostExp, deltaFatiguePrePostfMRI] = extract_subjective_fatigue_ratings(study_nm, subject_id, NS);

%% load all in a matrix/table
mega_mtrx_names = {'gal_age','gal_sex','gal_ISCE','gal_weight','gal_height','gal_BMI',...
    'gal_n_covid','gal_money','gal_avg_sleep','gal_prevDay_sleep','gal_prevDay_min_avg_sleep',...
    'gal_hexaco_hh','gal_hexaco_emotion','gal_hexaco_extraversion',...
    'gal_hexaco_agreeableness','gal_hexaco_consciousness','gal_hexaco_openness',...
    'salivary_TESTO1','salivary_TESTO4',...
    'salivary_CORT1','salivary_CORT2','salivary_CORT3','salivary_CORT4',...
    'salivary_TESTO_div_CORT1','salivary_TESTO_div_CORT4',...
    'salivary_IL1b','salivary_IL6','salivary_IL18',...
    'wholeB_NAM','wholeB_NMN','wholeB_NR','wholeB_NAD','wholeB_NADH',...
    'wholeB_NADP','wholeB_NADPH','wholeB_MeNam','wholeB_MeXPY',...
    'wholeB_NAD_div_NADH','wholeB_NADP_div_NADPH','wholeB_tNAD',...
    'wholeB_tNAD_with_precursors','wholeB_tNAD_with_byproducts',...
    'wholeB_NAD_byproducts',...
    'plasma_aa_Ala','plasma_aa_Arg','plasma_aa_Asp','plasma_aa_Asn',...
    'plasma_aa_Glu','plasma_aa_Gln','plasma_aa_Gly','plasma_aa_His',...
    'plasma_aa_Ile','plasma_aa_Leu','plasma_aa_Lys','plasma_aa_Met',...
    'plasma_aa_Phe','plasma_aa_Pro','plasma_aa_Ser','plasma_aa_Thr',...
    'plasma_aa_Tyr','plasma_aa_Val','plasma_aa_Tau','plasma_aa_T4HPro',...
    'plasma_aa_x3MethylHis','plasma_aa_x1MethylHis',...
    'plasma_Lac','plasma_fa_AcetoAcetate','plasma_fa_BHB','plasma_fa_Acetate',...
    'plasma_fa_Propionate', 'plasma_fa_IsoButyrate', 'plasma_fa_Butyrate',...
    'plasma_fa_x2MethylButyrate', 'plasma_fa_isovalerate', 'plasma_fa_valerate',...
    'plasma_fa_hexanoicA', 'plasma_fa_octanoicA', 'plasma_fa_decanoicA',...
    'plasma_fa_AABA', 'plasma_fa_L_Citrulline', 'plasma_fa_SDMA', 'plasma_fa_ADMA',...
    'plasma_fa_sarcosine',...
    'brainM_dmPFC_Asp', 'brainM_dmPFC_GABA', 'brainM_dmPFC_Gln', 'brainM_dmPFC_Glu',...
    'brainM_dmPFC_Glx', 'brainM_dmPFC_Gln_div_Glu', 'brainM_dmPFC_GSH', 'brainM_dmPFC_Gly',...
    'brainM_dmPFC_Ins', 'brainM_dmPFC_Lac', 'brainM_dmPFC_Tau', 'brainM_dmPFC_NAA',...
    'brainM_dmPFC_NAAG', 'brainM_dmPFC_GPC_PCho', 'brainM_dmPFC_Cr_PCr',...
    'brainM_aIns_Asp', 'brainM_aIns_GABA', 'brainM_aIns_Gln', 'brainM_aIns_Glu',...
    'brainM_aIns_Glx', 'brainM_aIns_Gln_div_Glu', 'brainM_aIns_GSH', 'brainM_aIns_Gly',...
    'brainM_aIns_Ins', 'brainM_aIns_Lac', 'brainM_aIns_Tau', 'brainM_aIns_NAA',...
    'brainM_aIns_NAAG', 'brainM_aIns_GPC_PCho', 'brainM_aIns_Cr_PCr',...
    'behavior_questionnaires_stress_anxiety_PSS14',...
    'behavior_questionnaires_stress_anxiety_STAI_T',...
    'behavior_questionnaires_stress_anxiety_SIAS',...
    'behavior_questionnaires_stress_anxiety_CTQ_emotionalA',...
    'behavior_questionnaires_stress_anxiety_CTQ_physicalA',...
    'behavior_questionnaires_stress_anxiety_CTQ_sexA',...
    'behavior_questionnaires_stress_anxiety_CTQ_minDenial',...
    'behavior_questionnaires_stress_anxiety_CTQ_emotionalN',...
    'behavior_questionnaires_stress_anxiety_CTQ_physicalN',...
    'behavior_questionnaires_motivation_MADRS_S',...
    'behavior_questionnaires_motivation_BIS_NPI',...
    'behavior_questionnaires_motivation_BIS_MI',...
    'behavior_questionnaires_motivation_BIS_AI',...
    'behavior_questionnaires_motivation_Lars_e_AI_EverydayProd',...
    'behavior_questionnaires_motivation_Lars_e_AI_Init',...
    'behavior_questionnaires_motivation_Lars_e_IC_Interest',...
    'behavior_questionnaires_motivation_Lars_e_IC_Novelty',...
    'behavior_questionnaires_motivation_Lars_e_IC_Motiv',...
    'behavior_questionnaires_motivation_Lars_e_IC_Social',...
    'behavior_questionnaires_motivation_Lars_e_ActionInit',...
    'behavior_questionnaires_motivation_Lars_e_IntellectCuriosity',...
    'behavior_questionnaires_motivation_Lars_e_EmotResp',...
    'behavior_questionnaires_motivation_Lars_e_SelfAwareness',...
    'behavior_questionnaires_motivation_MPSTEFS_physical',...
    'behavior_questionnaires_motivation_MPSTEFS_mental',...
    'behavior_questionnaires_motivation_JPIR',...
    'behavior_questionnaires_motivation_SHAP',...
    'behavior_questionnaires_motivation_SPSRQ_R',...
    'behavior_questionnaires_motivation_SPSRQ_P',...
    'behavior_questionnaires_motivation_IPAQ',...
    'behavior_questionnaires_motivation_IPAQinactivity',...
    'behavior_questionnaires_dominance_PRF_D',...
    'behavior_questionnaires_dominance_CI_enjCompet',...
    'behavior_questionnaires_dominance_CI_contentiousness',...
    'behavior_questionnaires_dominance_socialLadder',...
    
    'FFQ'
    
    'behavior_stress_ratings_S1', 'behavior_stress_ratings_S2',...
    'behavior_stress_ratings_S3', 'behavior_stress_ratings_S4',...
    'behavior_stress_ratings_S4_min_S1', 'behavior_stress_ratings_S4_min_S3',...
    'behavior_fatigue_ratings_F1', 'behavior_fatigue_ratings_F2',...
    'behavior_fatigue_ratings_F3', 'behavior_fatigue_ratings_F4',...
    'behavior_fatigue_ratings_F4_min_F1', 'behavior_fatigue_ratings_F4_min_F3',...
    'behavior_task_choices_HE', 'behavior_task_choices_HPE', 'behavior_task_choices_HME',...
    'behavior_task_prm_kR', 'behavior_task_prm_kP',...
    'behavior_task_prm_kEp', 'behavior_task_prm_kEm',...
    'behavior_task_prm_kFp', 'behavior_task_prm_kLm',...
    'behavior_task_prm_kBias'};
n_mega_mtrx_vars = length(mega_mtrx_names);
mega_mtrx = deal(NaN(n_mega_mtrx_vars, NS));
for iVar = 1:n_mega_mtrx_vars
    var_nm = mega_mtrx_names{iVar};
    switch var_nm
        case 'gal_age'
            mega_mtrx(iVar, :) = questionnaires.general.age;
        case 'gal_sex'
            mega_mtrx(iVar, :) = ;
        case 'gal_ISCE'
            mega_mtrx(iVar, :) = questionnaires.SES.ISCE;
        case 'gal_weight'
            mega_mtrx(iVar, :) = questionnaires.general.weight;
        case 'gal_height'
            mega_mtrx(iVar, :) = questionnaires.general.height;
        case 'gal_BMI'
            mega_mtrx(iVar, :) = questionnaires.general.BMI;
        case 'gal_n_covid'
            mega_mtrx(iVar, :) = questionnaires.general.n_covid;
        case 'gal_money'
            mega_mtrx(iVar, :) = questionnaires.SES.money;
        case 'gal_avg_sleep'
            mega_mtrx(iVar, :) = questionnaires.sleep.avg;
        case 'gal_prevDay_sleep'
            mega_mtrx(iVar, :) = questionnaires.sleep.prevDay;
        case 'gal_prevDay_min_avg_sleep'
            mega_mtrx(iVar, :) = questionnaires.sleep.prevDay_min_avg_sleep;
        case 'gal_hexaco_hh'
            mega_mtrx(iVar, :) = questionnaires.hexaco.honestyHumility;
        case 'gal_hexaco_emotion'
            mega_mtrx(iVar, :) = questionnaires.hexaco.emotion;
        case 'gal_hexaco_extraversion'
            mega_mtrx(iVar, :) = questionnaires.hexaco.extraversion;
        case 'gal_hexaco_agreeableness'
            mega_mtrx(iVar, :) = questionnaires.hexaco.agreeableness;
        case 'gal_hexaco_consciousness'
            mega_mtrx(iVar, :) = questionnaires.hexaco.consciousness;
        case 'gal_hexaco_openness'
            mega_mtrx(iVar, :) = questionnaires.hexaco.openness;
        case 'salivary_TESTO1'
            mega_mtrx(iVar, :) = TESTO_data.TESTO(1,:);
        case 'salivary_TESTO4'
            mega_mtrx(iVar, :) = TESTO_data.TESTO(4,:);
        case 'salivary_CORT1'
            mega_mtrx(iVar, :) = CORT_data.CORT(1,:);
        case 'salivary_CORT2'
            mega_mtrx(iVar, :) = CORT_data.CORT(2,:);
        case 'salivary_CORT3'
            mega_mtrx(iVar, :) = CORT_data.CORT(3,:);
        case 'salivary_CORT4'
            mega_mtrx(iVar, :) = CORT_data.CORT(4,:);
        case 'salivary_TESTO_div_CORT1'
            mega_mtrx(iVar, :) = TESTO_data.TESTO(1,:)./CORT_data.CORT(1,:);
        case 'salivary_TESTO_div_CORT4'
            mega_mtrx(iVar, :) = TESTO_data.TESTO(4,:)./CORT_data.CORT(4,:);
        case 'salivary_IL1b'
            mega_mtrx(iVar, :) = IL_data.IL1;
        case 'salivary_IL6'
            mega_mtrx(iVar, :) = IL_data.IL6;
        case 'salivary_IL18'
            mega_mtrx(iVar, :) = IL_data.IL18;
        case 'wholeB_NAM'
            mega_mtrx(iVar, :) = wholeBlood.Nam;
        case 'wholeB_NMN'
            mega_mtrx(iVar, :) = wholeBlood.NMN;
        case 'wholeB_NR'
            mega_mtrx(iVar, :) = wholeBlood.NR;
        case 'wholeB_NAD'
            mega_mtrx(iVar, :) = wholeBlood.NAD;
        case 'wholeB_NADH'
            mega_mtrx(iVar, :) = wholeBlood.NADH;
        case 'wholeB_NADP'
            mega_mtrx(iVar, :) = wholeBlood.NADP;
        case 'wholeB_NADPH'
            mega_mtrx(iVar, :) = wholeBlood.NADPH;
        case 'wholeB_MeNam'
            mega_mtrx(iVar, :) = wholeBlood.MeNam;
        case 'wholeB_MeXPY'
            mega_mtrx(iVar, :) = wholeBlood.MeXPY;
        case 'wholeB_NAD_div_NADH'
            mega_mtrx(iVar, :) = wholeBlood.NAD_div_NADH;
        case 'wholeB_NADP_div_NADPH'
            mega_mtrx(iVar, :) = wholeBlood.NADP_div_NADPH;
        case 'wholeB_tNAD'
            mega_mtrx(iVar, :) = wholeBlood.total_NAD;
        case 'wholeB_tNAD_with_precursors'
            mega_mtrx(iVar, :) = wholeBlood.total_NAD_with_precursors;
        case 'wholeB_tNAD_with_byproducts'
            mega_mtrx(iVar, :) = wholeBlood.total_NAD_with_byproducts;
        case 'wholeB_NAD_byproducts'
            mega_mtrx(iVar, :) = wholeBlood.total_NAD_byproducts;
        case 'plasma_aa_Ala'
            mega_mtrx(iVar, :) = plasmaM.Ala;
        case 'plasma_aa_Arg'
            mega_mtrx(iVar, :) = plasmaM.Arg;
        case 'plasma_aa_Asp'
            mega_mtrx(iVar, :) = plasmaM.Asp;
        case 'plasma_aa_Asn'
            mega_mtrx(iVar, :) = plasmaM.Asn;
        case 'plasma_aa_Glu'
            mega_mtrx(iVar, :) = plasmaM.Glu;
        case 'plasma_aa_Gln'
            mega_mtrx(iVar, :) = plasmaM.Gln;
        case 'plasma_aa_Gly'
            mega_mtrx(iVar, :) = plasmaM.Gly;
        case 'plasma_aa_His'
            mega_mtrx(iVar, :) = plasmaM.His;
        case 'plasma_aa_Ile'
            mega_mtrx(iVar, :) = plasmaM.Ile;
        case 'plasma_aa_Leu'
            mega_mtrx(iVar, :) = plasmaM.Leu;
        case 'plasma_aa_Lys'
            mega_mtrx(iVar, :) = plasmaM.Lys;
        case 'plasma_aa_Met'
            mega_mtrx(iVar, :) = plasmaM.Met;
        case 'plasma_aa_Phe'
            mega_mtrx(iVar, :) = plasmaM.Phe;
        case 'plasma_aa_Pro'
            mega_mtrx(iVar, :) = plasmaM.Pro;
        case 'plasma_aa_Ser'
            mega_mtrx(iVar, :) = plasmaM.Ser;
        case 'plasma_aa_Thr'
            mega_mtrx(iVar, :) = plasmaM.Thr;
        case 'plasma_aa_Tyr'
            mega_mtrx(iVar, :) = plasmaM.Tyr;
        case 'plasma_aa_Val'
            mega_mtrx(iVar, :) = plasmaM.Val;
        case 'plasma_aa_Tau'
            mega_mtrx(iVar, :) = plasmaM.Tau;
        case 'plasma_aa_T4HPro'
            mega_mtrx(iVar, :) = plasmaM.Trans_4_Hydroxyproline;
        case 'plasma_aa_x3MethylHis'
            mega_mtrx(iVar, :) = plasmaM.x3_MethylHistidine;
        case 'plasma_aa_x1MethylHis'
            mega_mtrx(iVar, :) = plasmaM.x1_Methylhistidine;
        case 'plasma_Lac'
            mega_mtrx(iVar, :) = plasmaM.Lac;
        case 'plasma_fa_AcetoAcetate'
            mega_mtrx(iVar, :) = plasmaM.Acetoacetate;
        case 'plasma_fa_BHB'
            mega_mtrx(iVar, :) = plasmaM.Hydroxybutyrate;
        case 'plasma_fa_Acetate'
            mega_mtrx(iVar, :) = plasmaM.Acetate;
        case 'plasma_fa_Propionate'
            mega_mtrx(iVar, :) = plasmaM.Propionate;
        case 'plasma_fa_IsoButyrate'
            mega_mtrx(iVar, :) = plasmaM.Isobutyrate;
        case 'plasma_fa_Butyrate'
            mega_mtrx(iVar, :) = plasmaM.Butyrate;
        case 'plasma_fa_x2MethylButyrate'
            mega_mtrx(iVar, :) = plasmaM.x2_Methylbutyrate;
        case 'plasma_fa_isovalerate'
            mega_mtrx(iVar, :) = plasmaM.Isovalerate;
        case 'plasma_fa_valerate'
            mega_mtrx(iVar, :) = plasmaM.Valerate;
        case 'plasma_fa_hexanoicA'
            mega_mtrx(iVar, :) = plasmaM.HexanoicAcid;
        case 'plasma_fa_octanoicA'
            mega_mtrx(iVar, :) = plasmaM.OctanoicAcid;
        case 'plasma_fa_decanoicA'
            mega_mtrx(iVar, :) = plasmaM.DecanoicAcid;
        case 'plasma_fa_AABA'
            mega_mtrx(iVar, :) = plasmaM.a_AminobutyricA_AABA_;
        case 'plasma_fa_L_Citrulline'
            mega_mtrx(iVar, :) = plasmaM.L_Citrulline;
        case 'plasma_fa_SDMA'
            mega_mtrx(iVar, :) = plasmaM.SymmetricDimethylarginine_SDMA_;
        case 'plasma_fa_ADMA'
            mega_mtrx(iVar, :) = plasmaM.AsymmetricDimethylarginine_ADMA_;
        case 'plasma_fa_sarcosine'
            mega_mtrx(iVar, :) = plasmaM.Sarcosine;
        case 'brainM_dmPFC_Asp'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Asp;
        case 'brainM_dmPFC_GABA'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.GABA;
        case 'brainM_dmPFC_Gln'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Gln;
        case 'brainM_dmPFC_Glu'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Glu;
        case 'brainM_dmPFC_Glx'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Glx;
        case 'brainM_dmPFC_Gln_div_Glu'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Gln_div_Glu;
        case 'brainM_dmPFC_GSH'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.GSH;
        case 'brainM_dmPFC_Gly'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Gly;
        case 'brainM_dmPFC_Ins'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Ins;
        case 'brainM_dmPFC_Lac'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Lac;
        case 'brainM_dmPFC_Tau'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Tau;
        case 'brainM_dmPFC_NAA'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.NAA;
        case 'brainM_dmPFC_NAAG'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.NAAG;
        case 'brainM_dmPFC_GPC_PCho'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.GPC_PCho;
        case 'brainM_dmPFC_Cr_PCr'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Cr_PCr;
        case 'brainM_aIns_Asp'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Asp;
        case 'brainM_aIns_GABA'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.GABA;
        case 'brainM_aIns_Gln'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Gln;
        case 'brainM_aIns_Glu'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Glu;
        case 'brainM_aIns_Glx'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Glx;
        case 'brainM_aIns_Gln_div_Glu'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Gln_div_Glu;
        case 'brainM_aIns_GSH'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.GSH;
        case 'brainM_aIns_Gly'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Gly;
        case 'brainM_aIns_Ins'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Ins;
        case 'brainM_aIns_Lac'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Lac;
        case 'brainM_aIns_Tau'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Tau;
        case 'brainM_aIns_NAA'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.NAA;
        case 'brainM_aIns_NAAG'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.NAAG;
        case 'brainM_aIns_GPC_PCho'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.GPC_PCho;
        case 'brainM_aIns_Cr_PCr'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Cr_PCr;
        case 'behavior_questionnaires_stress_anxiety_PSS14'
            mega_mtrx(iVar, :) = questionnaires.stress_anxiety.PSS14;
        case 'behavior_questionnaires_stress_anxiety_STAI_T'
            mega_mtrx(iVar, :) = questionnaires.stress_anxiety.STAI_T;
        case 'behavior_questionnaires_stress_anxiety_SIAS'
            mega_mtrx(iVar, :) = questionnaires.stress_anxiety.SIAS;
        case 'behavior_questionnaires_stress_anxiety_CTQ_emotionalA'
            mega_mtrx(iVar, :) = questionnaires.CTQ.emotionalAbuse;
        case 'behavior_questionnaires_stress_anxiety_CTQ_physicalA'
            mega_mtrx(iVar, :) = questionnaires.CTQ.physicalAbuse;
        case 'behavior_questionnaires_stress_anxiety_CTQ_sexA'
            mega_mtrx(iVar, :) = questionnaires.CTQ.sexualAbuse;
        case 'behavior_questionnaires_stress_anxiety_CTQ_minDenial'
            mega_mtrx(iVar, :) = questionnaires.CTQ.minimizationDenial;
        case 'behavior_questionnaires_stress_anxiety_CTQ_emotionalN'
            mega_mtrx(iVar, :) = questionnaires.CTQ.emotionalNeglect;
        case 'behavior_questionnaires_stress_anxiety_CTQ_physicalN'
            mega_mtrx(iVar, :) = questionnaires.CTQ.physicalNeglect;
        case 'behavior_questionnaires_motivation_MADRS_S'
            mega_mtrx(iVar, :) = questionnaires.Motivation.MADRS_S;
        case 'behavior_questionnaires_motivation_BIS_NPI'
            mega_mtrx(iVar, :) = questionnaires.Motivation.BIS_NPI;
        case 'behavior_questionnaires_motivation_BIS_MI'
            mega_mtrx(iVar, :) = questionnaires.Motivation.BIS_MI;
        case 'behavior_questionnaires_motivation_BIS_AI'
            mega_mtrx(iVar, :) = questionnaires.Motivation.BIS_AI;
        case 'behavior_questionnaires_motivation_Lars_e_AI_EverydayProd'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_AI_EverydayProd;
        case 'behavior_questionnaires_motivation_Lars_e_AI_Init'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_AI_Init;
        case 'behavior_questionnaires_motivation_Lars_e_IC_Interest'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IC_Interest;
        case 'behavior_questionnaires_motivation_Lars_e_IC_Novelty'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IC_Novelty;
        case 'behavior_questionnaires_motivation_Lars_e_IC_Motiv'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IC_Motiv;
        case 'behavior_questionnaires_motivation_Lars_e_IC_Social'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IC_Social;
        case 'behavior_questionnaires_motivation_Lars_e_ActionInit'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_ActionInit;
        case 'behavior_questionnaires_motivation_Lars_e_IntellectCuriosity'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IntellectCuriosity;
        case 'behavior_questionnaires_motivation_Lars_e_EmotResp'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_EmotResp;
        case 'behavior_questionnaires_motivation_Lars_e_SelfAwareness'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_SelfAwareness;
        case 'behavior_questionnaires_motivation_MPSTEFS_physical'
            mega_mtrx(iVar, :) = questionnaires.Motivation.MPSTEFS_physical;
        case 'behavior_questionnaires_motivation_MPSTEFS_mental'
            mega_mtrx(iVar, :) = questionnaires.Motivation.MPSTEFS_mental;
        case 'behavior_questionnaires_motivation_JPIR'
            mega_mtrx(iVar, :) = questionnaires.Motivation.JPIR;
        case 'behavior_questionnaires_motivation_SHAP'
            mega_mtrx(iVar, :) = questionnaires.Motivation.SHAP;
        case 'behavior_questionnaires_motivation_SPSRQ_R'
            mega_mtrx(iVar, :) = questionnaires.Motivation.SPSRQ_R;
        case 'behavior_questionnaires_motivation_SPSRQ_P'
            mega_mtrx(iVar, :) = questionnaires.Motivation.SPSRQ_P;
        case 'behavior_questionnaires_motivation_IPAQ'
            mega_mtrx(iVar, :) = questionnaires.Motivation.
        case 'behavior_questionnaires_motivation_IPAQinactivity'
            mega_mtrx(iVar, :) = questionnaires.Motivation.IPAQinactivity;
        case 'behavior_questionnaires_dominance_PRF_D'
            mega_mtrx(iVar, :) = questionnaires.dominance_compet.PRF_D;
        case 'behavior_questionnaires_dominance_CI_enjCompet'
            mega_mtrx(iVar, :) = questionnaires.dominance_compet.CI_enjCompet;
        case 'behavior_questionnaires_dominance_CI_contentiousness'
            mega_mtrx(iVar, :) = questionnaires.dominance_compet.CI_contentiousness;
        case 'behavior_questionnaires_dominance_socialLadder'
            mega_mtrx(iVar, :) = questionnaires.SES.ladder;

        case 'FFQ'
            
        case 'behavior_stress_ratings_S1'
            mega_mtrx(iVar, :) = preMRS_stress;
        case 'behavior_stress_ratings_S2'
            mega_mtrx(iVar, :) = postMRS_stress;
        case 'behavior_stress_ratings_S3'
            mega_mtrx(iVar, :) = prefMRI_stress;
        case 'behavior_stress_ratings_S4'
            mega_mtrx(iVar, :) = postfMRI_stress;
        case 'behavior_stress_ratings_S4_min_S1'
            mega_mtrx(iVar, :) = deltaStressPrePostExp;
        case 'behavior_stress_ratings_S4_min_S3'
            mega_mtrx(iVar, :) = deltaStressPrePostfMRI;
        case 'behavior_fatigue_ratings_F1'
            mega_mtrx(iVar, :) = preMRS_fatigue;
        case 'behavior_fatigue_ratings_F2'
            mega_mtrx(iVar, :) = postMRS_fatigue;
        case 'behavior_fatigue_ratings_F3'
            mega_mtrx(iVar, :) = prefMRI_fatigue;
        case 'behavior_fatigue_ratings_F4'
            mega_mtrx(iVar, :) = postfMRI_fatigue;
        case 'behavior_fatigue_ratings_F4_min_F1'
            mega_mtrx(iVar, :) = deltaFatiguePrePostExp;
        case 'behavior_fatigue_ratings_F4_min_F3'
            mega_mtrx(iVar, :) = deltaFatiguePrePostfMRI;
        case 'behavior_task_choices_HE'
            mega_mtrx(iVar, :) = choice_hE.
        case 'behavior_task_choices_HPE'
            mega_mtrx(iVar, :) = choice_hE.
        case 'behavior_task_choices_HME'
            mega_mtrx(iVar, :) = choice_hE.
        case 'behavior_task_prm_kR'
            mega_mtrx(iVar, :) = prm.kR;
        case 'behavior_task_prm_kP'
            mega_mtrx(iVar, :) = prm.kP;
        case 'behavior_task_prm_kEp'
            mega_mtrx(iVar, :) = prm.kEp;
        case 'behavior_task_prm_kEm'
            mega_mtrx(iVar, :) = prm.kEm;
        case 'behavior_task_prm_kFp'
            mega_mtrx(iVar, :) = prm.kFp;
        case 'behavior_task_prm_kLm'
            mega_mtrx(iVar, :) = prm.kLm;
        case 'behavior_task_prm_kBias'
            mega_mtrx(iVar, :) = prm.kBias;
    end
end % variable loop

%% perform correlation matrix between all variables
[corr_mtrx, pval_mtrx] = deal(NaN(n_mega_mtrx_vars, n_mega_mtrx_vars));
for iVar = 1:n_mega_mtrx_vars
    for jVar = 1:n_mega_mtrx_vars
        [corr_mtrx(iVar, jVar),...
            pval_mtrx(iVar, jVar)] = corr(mega_mtrx(iVar,:)', mega_mtrx(jVar,:)');
    end % variable loop
end % variable loop

%% perform generalized linear mixed model (GLME) + basic correlation matrix for Nestlé
tbl = ;
glme = fitglme(tbl, formula);

%% display corresponding figures
% figure parameters
[pSize, ~, col] = general_fig_prm;
% define which colormap you want to use (see full list here if you are not
% happy with the selection:
% https://ch.mathworks.com/help/matlab/ref/colormap.html)
% color_range_choices = 'hot';
% color_range_choices = 'turbo';
% color_range_choices = 'jet';
color_range_choices = redblue(45);

% correlation range
corr_range = [-0.65 0.65];

%% correlation matrix
fig;
subplot_hdl = subplot(1,1,1);
imagesc(corr_mtrx, corr_range);
colormap(subplot_hdl, color_range_choices);
cbar = colorbar;
cbar.Label.String = 'r';
xticks(1:n_mega_mtrx_vars);
xticklabels(mega_mtrx_names);
yticks(1:n_mega_mtrx_vars);
yticklabels(mega_mtrx_names);
% add stars in the graph if some correlations are significant
for iVar = 1:n_mega_mtrx_vars
    for jVar = 1:n_mega_mtrx_vars
        if pval_mtrx(iVar, jVar) <= 0.05
            if pval_mtrx(iVar, jVar) > 0.01 && pval_mtrx(iVar, jVar) <= 0.05
                pval_hdl = text(jVar, iVar, '*');
                warning('check iVar vs jVar');
            elseif pval_mtrx(iVar, jVar) > 0.001 && pval_mtrx(iVar, jVar) <= 0.01
                pval_hdl = text(jVar, iVar, '**');
            elseif pval_mtrx(iVar, jVar) <= 0.001
                pval_hdl = text(jVar, iVar, '***');
            end % p.value
            % adjust p.value parameters
            pval_hdl.Color = col.black;
            pval_hdl.FontSize = 70;
            pval_hdl.FontWeight = 'bold';
            pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
            pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
        end % when p.value is significant
    end % loop over Y variables
end % loop over X variables