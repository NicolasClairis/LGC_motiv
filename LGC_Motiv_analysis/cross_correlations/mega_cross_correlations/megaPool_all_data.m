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

%% working directory
saveFolder = fullfile('P:','boulot','postdoc_CarmenSandi',...
    'results','mega_correlation_matrix');

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
    'behavior_questionnaires_motivation_IPAQinactivity',...
    'behavior_questionnaires_dominance_PRF_D',...
    'behavior_questionnaires_dominance_CI_enjCompet',...
    'behavior_questionnaires_dominance_CI_contentiousness',...
    'behavior_questionnaires_dominance_socialLadder',...
    'FFQ_raw_energy_kCal','FFQ_raw_energy_kJ','FFQ_raw_energy_kCalWithFibers','FFQ_raw_energy_kJ_withFibers',...
    'FFQ_raw_proteins_aa_proteins','FFQ_raw_proteins_aa_essentialAA','FFQ_raw_proteins_aa_nonEssentialAA',...
    'FFQ_raw_proteins_aa_Ile','FFQ_raw_proteins_aa_Leu','FFQ_raw_proteins_aa_Lys',...
    'FFQ_raw_proteins_aa_Met','FFQ_raw_proteins_aa_Cys','FFQ_raw_proteins_aa_Phe',...
    'FFQ_raw_proteins_aa_Tyr','FFQ_raw_proteins_aa_Thr','FFQ_raw_proteins_aa_Trp',...
    'FFQ_raw_proteins_aa_Val','FFQ_raw_proteins_aa_Arg','FFQ_raw_proteins_aa_His',...
    'FFQ_raw_proteins_aa_Ala','FFQ_raw_proteins_aa_Asp','FFQ_raw_proteins_aa_Glu',...
    'FFQ_raw_proteins_aa_Gly','FFQ_raw_proteins_aa_Pro','FFQ_raw_proteins_aa_Ser',...
    'FFQ_raw_vitamins_ARE','FFQ_raw_vitamins_AR','FFQ_raw_vitamins_ABC',...
    'FFQ_raw_vitamins_D','FFQ_raw_vitamins_ETE','FFQ_raw_vitamins_EAT',...
    'FFQ_raw_vitamins_K','FFQ_raw_vitamins_B1','FFQ_raw_vitamins_B2',...
    'FFQ_raw_vitamins_B3_niacin','FFQ_raw_vitamins_B3_NE','FFQ_raw_vitamins_B5',...
    'FFQ_raw_vitamins_B6','FFQ_raw_vitamins_B7','FFQ_raw_vitamins_B9',...
    'FFQ_raw_vitamins_B12','FFQ_raw_vitamins_C',...
    'FFQ_raw_minerals_ashes','FFQ_raw_minerals_Na','FFQ_raw_minerals_K',...
    'FFQ_raw_minerals_Ca','FFQ_raw_minerals_Mg','FFQ_raw_minerals_P',...
    'FFQ_raw_minerals_S','FFQ_raw_minerals_Cl','FFQ_raw_minerals_Fe',...
    'FFQ_raw_minerals_Zn','FFQ_raw_minerals_Cu','FFQ_raw_minerals_Mn',...
    'FFQ_raw_minerals_F','FFQ_raw_minerals_I','FFQ_raw_minerals_NaCl',...
    'FFQ_raw_carbs_absorbableCarbohydrates','FFQ_raw_carbs_mannitol','FFQ_raw_carbs_sorbitol',...
    'FFQ_raw_carbs_xylitol','FFQ_raw_carbs_sugarAlcohols','FFQ_raw_carbs_Glc',...
    'FFQ_raw_carbs_fructose','FFQ_raw_carbs_galactose','FFQ_raw_carbs_monosaccharides',...
    'FFQ_raw_carbs_saccharose','FFQ_raw_carbs_maltose','FFQ_raw_carbs_lactose',...
    'FFQ_raw_carbs_disaccharides','FFQ_raw_carbs_sugar','FFQ_raw_carbs_reabsorbableOligosaccharides',...
    'FFQ_raw_carbs_nonReabsorbableOligosaccharides','FFQ_raw_carbs_glycogen_animalStarch',...
    'FFQ_raw_carbs_starch','FFQ_raw_carbs_polysaccharides','FFQ_raw_carbs_polypentoses',...
    'FFQ_raw_carbs_polyhexoses','FFQ_raw_carbs_polyuronicA',...
    'FFQ_raw_fibers_fibers','FFQ_raw_fibers_solubleFibers','FFQ_raw_fibers_insolubleFibers',...
    'FFQ_raw_fibers_cellulose','FFQ_raw_fibers_lignine',...
    'FFQ_raw_fa_fat','FFQ_raw_fa_monoInsaturated','FFQ_raw_fa_polyInsaturated',...
    'FFQ_raw_fa_shortChain','FFQ_raw_fa_mediumChain','FFQ_raw_fa_longChain',...
    'FFQ_raw_fa_butanoicA','FFQ_raw_fa_hexanoicA','FFQ_raw_fa_octanoicA','FFQ_raw_fa_decanoicA_capricA',...
    'FFQ_raw_fa_dodecanoicA','FFQ_raw_fa_tetradecanoicA','FFQ_raw_fa_pentadecanoicA',...
    'FFQ_raw_fa_hexadecanoicA','FFQ_raw_fa_heptadecanoicA','FFQ_raw_fa_octadecanoicA',...
    'FFQ_raw_fa_eicosanoicA','FFQ_raw_fa_decanoicA_behenicA','FFQ_raw_fa_tetracosanoicA','FFQ_raw_fa_saturatedFA',...
    'FFQ_raw_fa_tetradecenoicA','FFQ_raw_fa_pentadecenoicA','FFQ_raw_fa_hexadecenoicA',...
    'FFQ_raw_fa_heptadecenoicA','FFQ_raw_fa_octadecenoicA','FFQ_raw_fa_eicosenoicA',...
    'FFQ_raw_fa_decosenoicA','FFQ_raw_fa_tetracosenoicA','FFQ_raw_fa_hexadecadienoicA',...
    'FFQ_raw_fa_hexacedatetraenoicA','FFQ_raw_fa_octadecadienoicA_linoleicA','FFQ_raw_fa_octadecatrienoicA_linolenicA',...
    'FFQ_raw_fa_octradecatetraenoicA','FFQ_raw_fa_nonAdecatrienoicA','FFQ_raw_fa_eicosadienoicA',...
    'FFQ_raw_fa_eicosatrienoicA','FFQ_raw_fa_eicosatetraenoicA_arachidonicA','FFQ_raw_fa_eicosapentaenoicA',...
    'FFQ_raw_fa_docosadienoicA','FFQ_raw_fa_docosatrienoicA','FFQ_raw_fa_docosatetraenoicA',...
    'FFQ_raw_fa_docosapentaenoicA','FFQ_raw_fa_docosahexaenoicA','FFQ_raw_fa_omegas3',...
    'FFQ_raw_fa_omegas6','FFQ_raw_fa_glycerine_lipids','FFQ_raw_fa_cholesterol',...
    'FFQ_raw_fa_GFPS_PolyI_div_Sat_fa_ratio',...
    'FFQ_raw_water','FFQ_raw_organicA','FFQ_raw_alcohol','FFQ_raw_uricA',...
    'FFQ_raw_purin','FFQ_raw_bread',...
    'FFQ_norm_energy_kCal','FFQ_norm_energy_kJ','FFQ_norm_energy_kCalWithFibers','FFQ_norm_energy_kJ_withFibers',...
    'FFQ_norm_proteins_aa_proteins','FFQ_norm_proteins_aa_essentialAA','FFQ_norm_proteins_aa_nonEssentialAA',...
    'FFQ_norm_proteins_aa_Ile','FFQ_norm_proteins_aa_Leu','FFQ_norm_proteins_aa_Lys',...
    'FFQ_norm_proteins_aa_Met','FFQ_norm_proteins_aa_Cys','FFQ_norm_proteins_aa_Phe',...
    'FFQ_norm_proteins_aa_Tyr','FFQ_norm_proteins_aa_Thr','FFQ_norm_proteins_aa_Trp',...
    'FFQ_norm_proteins_aa_Val','FFQ_norm_proteins_aa_Arg','FFQ_norm_proteins_aa_His',...
    'FFQ_norm_proteins_aa_Ala','FFQ_norm_proteins_aa_Asp','FFQ_norm_proteins_aa_Glu',...
    'FFQ_norm_proteins_aa_Gly','FFQ_norm_proteins_aa_Pro','FFQ_norm_proteins_aa_Ser',...
    'FFQ_norm_vitamins_ARE','FFQ_norm_vitamins_AR','FFQ_norm_vitamins_ABC',...
    'FFQ_norm_vitamins_D','FFQ_norm_vitamins_ETE','FFQ_norm_vitamins_EAT',...
    'FFQ_norm_vitamins_K','FFQ_norm_vitamins_B1','FFQ_norm_vitamins_B2',...
    'FFQ_norm_vitamins_B3_niacin','FFQ_norm_vitamins_B3_NE','FFQ_norm_vitamins_B5',...
    'FFQ_norm_vitamins_B6','FFQ_norm_vitamins_B7','FFQ_norm_vitamins_B9',...
    'FFQ_norm_vitamins_B12','FFQ_norm_vitamins_C',...
    'FFQ_norm_minerals_ashes','FFQ_norm_minerals_Na','FFQ_norm_minerals_K',...
    'FFQ_norm_minerals_Ca','FFQ_norm_minerals_Mg','FFQ_norm_minerals_P',...
    'FFQ_norm_minerals_S','FFQ_norm_minerals_Cl','FFQ_norm_minerals_Fe',...
    'FFQ_norm_minerals_Zn','FFQ_norm_minerals_Cu','FFQ_norm_minerals_Mn',...
    'FFQ_norm_minerals_F','FFQ_norm_minerals_I','FFQ_norm_minerals_NaCl',...
    'FFQ_norm_carbs_absorbableCarbohydrates','FFQ_norm_carbs_mannitol','FFQ_norm_carbs_sorbitol',...
    'FFQ_norm_carbs_xylitol','FFQ_norm_carbs_sugarAlcohols','FFQ_norm_carbs_Glc',...
    'FFQ_norm_carbs_fructose','FFQ_norm_carbs_galactose','FFQ_norm_carbs_monosaccharides',...
    'FFQ_norm_carbs_saccharose','FFQ_norm_carbs_maltose','FFQ_norm_carbs_lactose',...
    'FFQ_norm_carbs_disaccharides','FFQ_norm_carbs_sugar','FFQ_norm_carbs_reabsorbableOligosaccharides',...
    'FFQ_norm_carbs_nonReabsorbableOligosaccharides','FFQ_norm_carbs_glycogen_animalStarch',...
    'FFQ_norm_carbs_starch','FFQ_norm_carbs_polysaccharides','FFQ_norm_carbs_polypentoses',...
    'FFQ_norm_carbs_polyhexoses','FFQ_norm_carbs_polyuronicA',...
    'FFQ_norm_fibers_fibers','FFQ_norm_fibers_solubleFibers','FFQ_norm_fibers_insolubleFibers',...
    'FFQ_norm_fibers_cellulose','FFQ_norm_fibers_lignine',...
    'FFQ_norm_fa_fat','FFQ_norm_fa_monoInsaturated','FFQ_norm_fa_polyInsaturated',...
    'FFQ_norm_fa_shortChain','FFQ_norm_fa_mediumChain','FFQ_norm_fa_longChain',...
    'FFQ_norm_fa_butanoicA','FFQ_norm_fa_hexanoicA','FFQ_norm_fa_octanoicA','FFQ_norm_fa_decanoicA_capricA',...
    'FFQ_norm_fa_dodecanoicA','FFQ_norm_fa_tetradecanoicA','FFQ_norm_fa_pentadecanoicA',...
    'FFQ_norm_fa_hexadecanoicA','FFQ_norm_fa_heptadecanoicA','FFQ_norm_fa_octadecanoicA',...
    'FFQ_norm_fa_eicosanoicA','FFQ_norm_fa_decanoicA_behenicA','FFQ_norm_fa_tetracosanoicA','FFQ_norm_fa_saturatedFA',...
    'FFQ_norm_fa_tetradecenoicA','FFQ_norm_fa_pentadecenoicA','FFQ_norm_fa_hexadecenoicA',...
    'FFQ_norm_fa_heptadecenoicA','FFQ_norm_fa_octadecenoicA','FFQ_norm_fa_eicosenoicA',...
    'FFQ_norm_fa_decosenoicA','FFQ_norm_fa_tetracosenoicA','FFQ_norm_fa_hexadecadienoicA',...
    'FFQ_norm_fa_hexacedatetraenoicA','FFQ_norm_fa_octadecadienoicA_linoleicA','FFQ_norm_fa_octadecatrienoicA_linolenicA',...
    'FFQ_norm_fa_octradecatetraenoicA','FFQ_norm_fa_nonAdecatrienoicA','FFQ_norm_fa_eicosadienoicA',...
    'FFQ_norm_fa_eicosatrienoicA','FFQ_norm_fa_eicosatetraenoicA_arachidonicA','FFQ_norm_fa_eicosapentaenoicA',...
    'FFQ_norm_fa_docosadienoicA','FFQ_norm_fa_docosatrienoicA','FFQ_norm_fa_docosatetraenoicA',...
    'FFQ_norm_fa_docosapentaenoicA','FFQ_norm_fa_docosahexaenoicA','FFQ_norm_fa_omegas3',...
    'FFQ_norm_fa_omegas6','FFQ_norm_fa_glycerine_lipids','FFQ_norm_fa_cholesterol',...
    'FFQ_norm_fa_GFPS_PolyI_div_Sat_fa_ratio',...
    'FFQ_norm_water','FFQ_norm_organicA','FFQ_norm_alcohol','FFQ_norm_uricA',...
    'FFQ_norm_purin','FFQ_norm_bread',...
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
category = cell(1,n_mega_mtrx_vars);
for iVar = 1:n_mega_mtrx_vars
    var_nm = mega_mtrx_names{iVar};
    switch var_nm
        case 'gal_age'
            mega_mtrx(iVar, :) = questionnaires.general.age;
            category{iVar} = 'age';
        case 'gal_sex'
            mega_mtrx(iVar, :) = questionnaires.general.sex;
            category{iVar} = 'sex';
        case 'gal_ISCE'
            mega_mtrx(iVar, :) = questionnaires.SES.ISCE;
            category{iVar} = 'ISCE';
        case 'gal_weight'
            mega_mtrx(iVar, :) = questionnaires.general.weight;
            category{iVar} = 'weight';
        case 'gal_height'
            mega_mtrx(iVar, :) = questionnaires.general.height;
            category{iVar} = 'height';
        case 'gal_BMI'
            mega_mtrx(iVar, :) = questionnaires.general.BMI;
            category{iVar} = 'BMI';
        case 'gal_n_covid'
            mega_mtrx(iVar, :) = questionnaires.general.n_covid;
            category{iVar} = 'n_covid';
        case 'gal_money'
            mega_mtrx(iVar, :) = questionnaires.SES.money;
            category{iVar} = 'money';
        case 'gal_avg_sleep'
            mega_mtrx(iVar, :) = questionnaires.sleep.avg;
            category{iVar} = 'sleep';
        case 'gal_prevDay_sleep'
            mega_mtrx(iVar, :) = questionnaires.sleep.prevDay;
            category{iVar} = 'sleep';
        case 'gal_prevDay_min_avg_sleep'
            mega_mtrx(iVar, :) = questionnaires.sleep.prevDay_min_avg_sleep;
            category{iVar} = 'sleep';
        case 'gal_hexaco_hh'
            mega_mtrx(iVar, :) = questionnaires.hexaco.honestyHumility;
            category{iVar} = 'hexaco';
        case 'gal_hexaco_emotion'
            mega_mtrx(iVar, :) = questionnaires.hexaco.emotion;
            category{iVar} = 'hexaco';
        case 'gal_hexaco_extraversion'
            mega_mtrx(iVar, :) = questionnaires.hexaco.extraversion;
            category{iVar} = 'hexaco';
        case 'gal_hexaco_agreeableness'
            mega_mtrx(iVar, :) = questionnaires.hexaco.agreeableness;
            category{iVar} = 'hexaco';
        case 'gal_hexaco_consciousness'
            mega_mtrx(iVar, :) = questionnaires.hexaco.consciousness;
            category{iVar} = 'hexaco';
        case 'gal_hexaco_openness'
            mega_mtrx(iVar, :) = questionnaires.hexaco.openness;
            category{iVar} = 'hexaco';
        case 'salivary_TESTO1'
            mega_mtrx(iVar, :) = TESTO_data.TESTO(1,:);
            category{iVar} = 'salivary_Testosterone';
        case 'salivary_TESTO4'
            mega_mtrx(iVar, :) = TESTO_data.TESTO(4,:);
            category{iVar} = 'salivary_Testosterone';
        case 'salivary_CORT1'
            mega_mtrx(iVar, :) = CORT_data.CORT(1,:);
            category{iVar} = 'salivary_Cortisol';
        case 'salivary_CORT2'
            mega_mtrx(iVar, :) = CORT_data.CORT(2,:);
            category{iVar} = 'salivary_Cortisol';
        case 'salivary_CORT3'
            mega_mtrx(iVar, :) = CORT_data.CORT(3,:);
            category{iVar} = 'salivary_Cortisol';
        case 'salivary_CORT4'
            mega_mtrx(iVar, :) = CORT_data.CORT(4,:);
            category{iVar} = 'salivary_Cortisol';
        case 'salivary_TESTO_div_CORT1'
            mega_mtrx(iVar, :) = TESTO_data.TESTO(1,:)./CORT_data.CORT(1,:);
            category{iVar} = 'salivary_Testo_div_Cort';
        case 'salivary_TESTO_div_CORT4'
            mega_mtrx(iVar, :) = TESTO_data.TESTO(4,:)./CORT_data.CORT(4,:);
            category{iVar} = 'salivary_Testo_div_Cort';
        case 'salivary_IL1b'
            mega_mtrx(iVar, :) = IL_data.IL1b;
            category{iVar} = 'salivary_IL';
        case 'salivary_IL6'
            mega_mtrx(iVar, :) = IL_data.IL6;
            category{iVar} = 'salivary_IL';
        case 'salivary_IL18'
            mega_mtrx(iVar, :) = IL_data.IL18;
            category{iVar} = 'salivary_IL';
        case 'wholeB_NAM'
            mega_mtrx(iVar, :) = wholeBlood.Nam;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_NMN'
            mega_mtrx(iVar, :) = wholeBlood.NMN;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_NR'
            mega_mtrx(iVar, :) = wholeBlood.NR;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_NAD'
            mega_mtrx(iVar, :) = wholeBlood.NAD;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_NADH'
            mega_mtrx(iVar, :) = wholeBlood.NADH;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_NADP'
            mega_mtrx(iVar, :) = wholeBlood.NADP;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_NADPH'
            mega_mtrx(iVar, :) = wholeBlood.NADPH;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_MeNam'
            mega_mtrx(iVar, :) = wholeBlood.MeNam;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_MeXPY'
            mega_mtrx(iVar, :) = wholeBlood.MeXPY;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_NAD_div_NADH'
            mega_mtrx(iVar, :) = wholeBlood.NAD_div_NADH;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_NADP_div_NADPH'
            mega_mtrx(iVar, :) = wholeBlood.NADP_div_NADPH;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_tNAD'
            mega_mtrx(iVar, :) = wholeBlood.total_NAD;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_tNAD_with_precursors'
            mega_mtrx(iVar, :) = wholeBlood.total_NAD_with_precursors;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_tNAD_with_byproducts'
            mega_mtrx(iVar, :) = wholeBlood.total_NAD_with_byproducts;
            category{iVar} = 'wholeB_NADomics';
        case 'wholeB_NAD_byproducts'
            mega_mtrx(iVar, :) = wholeBlood.total_NAD_byproducts;
            category{iVar} = 'wholeB_NADomics';
        case 'plasma_aa_Ala'
            mega_mtrx(iVar, :) = plasmaM.Ala;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Arg'
            mega_mtrx(iVar, :) = plasmaM.Arg;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Asp'
            mega_mtrx(iVar, :) = plasmaM.Asp;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Asn'
            mega_mtrx(iVar, :) = plasmaM.Asn;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Glu'
            mega_mtrx(iVar, :) = plasmaM.Glu;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Gln'
            mega_mtrx(iVar, :) = plasmaM.Gln;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Gly'
            mega_mtrx(iVar, :) = plasmaM.Gly;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_His'
            mega_mtrx(iVar, :) = plasmaM.His;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Ile'
            mega_mtrx(iVar, :) = plasmaM.Ile;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Leu'
            mega_mtrx(iVar, :) = plasmaM.Leu;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Lys'
            mega_mtrx(iVar, :) = plasmaM.Lys;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Met'
            mega_mtrx(iVar, :) = plasmaM.Met;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Phe'
            mega_mtrx(iVar, :) = plasmaM.Phe;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Pro'
            mega_mtrx(iVar, :) = plasmaM.Pro;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Ser'
            mega_mtrx(iVar, :) = plasmaM.Ser;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Thr'
            mega_mtrx(iVar, :) = plasmaM.Thr;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Tyr'
            mega_mtrx(iVar, :) = plasmaM.Tyr;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Val'
            mega_mtrx(iVar, :) = plasmaM.Val;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_Tau'
            mega_mtrx(iVar, :) = plasmaM.Tau;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_T4HPro'
            mega_mtrx(iVar, :) = plasmaM.Trans_4_Hydroxyproline;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_x3MethylHis'
            mega_mtrx(iVar, :) = plasmaM.x3_MethylHistidine;
            category{iVar} = 'plasma_aa';
        case 'plasma_aa_x1MethylHis'
            mega_mtrx(iVar, :) = plasmaM.x1_Methylhistidine;
            category{iVar} = 'plasma_aa';
        case 'plasma_Lac'
            mega_mtrx(iVar, :) = plasmaM.Lac;
            category{iVar} = 'plasma_Lac';
        case 'plasma_fa_AcetoAcetate'
            mega_mtrx(iVar, :) = plasmaM.Acetoacetate;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_BHB'
            mega_mtrx(iVar, :) = plasmaM.Hydroxybutyrate;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_Acetate'
            mega_mtrx(iVar, :) = plasmaM.Acetate;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_Propionate'
            mega_mtrx(iVar, :) = plasmaM.Propionate;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_IsoButyrate'
            mega_mtrx(iVar, :) = plasmaM.Isobutyrate;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_Butyrate'
            mega_mtrx(iVar, :) = plasmaM.Butyrate;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_x2MethylButyrate'
            mega_mtrx(iVar, :) = plasmaM.x2_Methylbutyrate;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_isovalerate'
            mega_mtrx(iVar, :) = plasmaM.Isovalerate;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_valerate'
            mega_mtrx(iVar, :) = plasmaM.Valerate;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_hexanoicA'
            mega_mtrx(iVar, :) = plasmaM.HexanoicAcid;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_octanoicA'
            mega_mtrx(iVar, :) = plasmaM.OctanoicAcid;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_decanoicA'
            mega_mtrx(iVar, :) = plasmaM.DecanoicAcid;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_AABA'
            mega_mtrx(iVar, :) = plasmaM.a_AminobutyricA_AABA_;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_L_Citrulline'
            mega_mtrx(iVar, :) = plasmaM.L_Citrulline;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_SDMA'
            mega_mtrx(iVar, :) = plasmaM.SymmetricDimethylarginine_SDMA_;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_ADMA'
            mega_mtrx(iVar, :) = plasmaM.AsymmetricDimethylarginine_ADMA_;
            category{iVar} = 'plasma_fa';
        case 'plasma_fa_sarcosine'
            mega_mtrx(iVar, :) = plasmaM.Sarcosine;
            category{iVar} = 'plasma_fa';
        case 'brainM_dmPFC_Asp'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Asp;
            category{iVar} = 'dmPFC_aa';
        case 'brainM_dmPFC_GABA'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.GABA;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_Gln'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Gln;
            category{iVar} = 'dmPFC_aa';
        case 'brainM_dmPFC_Glu'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Glu;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_Glx'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Glx;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_Gln_div_Glu'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Gln_div_Glu;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_GSH'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.GSH;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_Gly'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Gly;
            category{iVar} = 'dmPFC_aa';
        case 'brainM_dmPFC_Ins'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Ins;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_Lac'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Lac;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_Tau'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Tau;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_NAA'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.NAA;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_NAAG'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.NAAG;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_GPC_PCho'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.GPC_PCho;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_dmPFC_Cr_PCr'
            mega_mtrx(iVar, :) = brainMetabolites.dmPFC.Cr_PCr;
            category{iVar} = 'dmPFC_mb';
        case 'brainM_aIns_Asp'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Asp;
            category{iVar} = 'aIns_aa';
        case 'brainM_aIns_GABA'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.GABA;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_Gln'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Gln;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_Glu'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Glu;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_Glx'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Glx;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_Gln_div_Glu'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Gln_div_Glu;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_GSH'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.GSH;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_Gly'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Gly;
            category{iVar} = 'aIns_aa';
        case 'brainM_aIns_Ins'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Ins;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_Lac'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Lac;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_Tau'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Tau;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_NAA'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.NAA;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_NAAG'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.NAAG;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_GPC_PCho'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.GPC_PCho;
            category{iVar} = 'aIns_mb';
        case 'brainM_aIns_Cr_PCr'
            mega_mtrx(iVar, :) = brainMetabolites.aIns.Cr_PCr;
            category{iVar} = 'aIns_mb';
        case 'behavior_questionnaires_stress_anxiety_PSS14'
            mega_mtrx(iVar, :) = questionnaires.stress_anxiety.PSS14;
            category{iVar} = 'stress_anxiety';
        case 'behavior_questionnaires_stress_anxiety_STAI_T'
            mega_mtrx(iVar, :) = questionnaires.stress_anxiety.STAI_T;
            category{iVar} = 'stress_anxiety';
        case 'behavior_questionnaires_stress_anxiety_SIAS'
            mega_mtrx(iVar, :) = questionnaires.stress_anxiety.SIAS;
            category{iVar} = 'stress_anxiety';
        case 'behavior_questionnaires_stress_anxiety_CTQ_emotionalA'
            mega_mtrx(iVar, :) = questionnaires.CTQ.emotionalAbuse;
            category{iVar} = 'stress_anxiety';
        case 'behavior_questionnaires_stress_anxiety_CTQ_physicalA'
            mega_mtrx(iVar, :) = questionnaires.CTQ.physicalAbuse;
            category{iVar} = 'stress_anxiety';
        case 'behavior_questionnaires_stress_anxiety_CTQ_sexA'
            mega_mtrx(iVar, :) = questionnaires.CTQ.sexualAbuse;
            category{iVar} = 'stress_anxiety';
        case 'behavior_questionnaires_stress_anxiety_CTQ_minDenial'
            mega_mtrx(iVar, :) = questionnaires.CTQ.minimizationDenial;
            category{iVar} = 'stress_anxiety';
        case 'behavior_questionnaires_stress_anxiety_CTQ_emotionalN'
            mega_mtrx(iVar, :) = questionnaires.CTQ.emotionalNeglect;
            category{iVar} = 'stress_anxiety';
        case 'behavior_questionnaires_stress_anxiety_CTQ_physicalN'
            mega_mtrx(iVar, :) = questionnaires.CTQ.physicalNeglect;
            category{iVar} = 'stress_anxiety';
        case 'behavior_questionnaires_motivation_MADRS_S'
            mega_mtrx(iVar, :) = questionnaires.Motivation.MADRS_S;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_BIS_NPI'
            mega_mtrx(iVar, :) = questionnaires.Motivation.BIS_NPI;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_BIS_MI'
            mega_mtrx(iVar, :) = questionnaires.Motivation.BIS_MI;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_BIS_AI'
            mega_mtrx(iVar, :) = questionnaires.Motivation.BIS_AI;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_AI_EverydayProd'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_AI_EverydayProd;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_AI_Init'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_AI_Init;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_IC_Interest'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IC_Interest;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_IC_Novelty'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IC_Novelty;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_IC_Motiv'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IC_Motiv;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_IC_Social'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IC_Social;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_ActionInit'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_ActionInit;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_IntellectCuriosity'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_IntellectCuriosity;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_EmotResp'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_EmotResp;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_Lars_e_SelfAwareness'
            mega_mtrx(iVar, :) = questionnaires.Motivation.Lars_e_SelfAwareness;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_MPSTEFS_physical'
            mega_mtrx(iVar, :) = questionnaires.Motivation.MPSTEFS_physical;
            category{iVar} = 'fatigue';
        case 'behavior_questionnaires_motivation_MPSTEFS_mental'
            mega_mtrx(iVar, :) = questionnaires.Motivation.MPSTEFS_mental;
            category{iVar} = 'fatigue';
        case 'behavior_questionnaires_motivation_JPIR'
            mega_mtrx(iVar, :) = questionnaires.Motivation.JPIR;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_SHAP'
            mega_mtrx(iVar, :) = questionnaires.Motivation.SHAP;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_SPSRQ_R'
            mega_mtrx(iVar, :) = questionnaires.Motivation.SPSRQ_R;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_SPSRQ_P'
            mega_mtrx(iVar, :) = questionnaires.Motivation.SPSRQ_P;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_motivation_IPAQinactivity'
            mega_mtrx(iVar, :) = questionnaires.Motivation.IPAQinactivity;
            category{iVar} = 'motivation';
        case 'behavior_questionnaires_dominance_PRF_D'
            mega_mtrx(iVar, :) = questionnaires.dominance_compet.PRF_D;
            category{iVar} = 'dominance';
        case 'behavior_questionnaires_dominance_CI_enjCompet'
            mega_mtrx(iVar, :) = questionnaires.dominance_compet.CI_enjCompet;
            category{iVar} = 'dominance';
        case 'behavior_questionnaires_dominance_CI_contentiousness'
            mega_mtrx(iVar, :) = questionnaires.dominance_compet.CI_contentiousness;
            category{iVar} = 'dominance';
        case 'behavior_questionnaires_dominance_socialLadder'
            mega_mtrx(iVar, :) = questionnaires.SES.ladder;
            category{iVar} = 'dominance';
        case 'FFQ_raw_energy_kCal'
            mega_mtrx(iVar, :) = nutrition.raw.energy.kCal;
            category{iVar} = 'FFQ_raw_energy';
        case 'FFQ_raw_energy_kJ'
            mega_mtrx(iVar, :) = nutrition.raw.energy.kJ;
            category{iVar} = 'FFQ_raw_energy';
        case 'FFQ_raw_energy_kCalWithFibers'
            mega_mtrx(iVar, :) = nutrition.raw.energy.kCal_withFibers;
            category{iVar} = 'FFQ_raw_energy';
        case 'FFQ_raw_energy_kJ_withFibers'
            mega_mtrx(iVar, :) = nutrition.raw.energy.kJ_withFibers;
            category{iVar} = 'FFQ_raw_energy';
        case 'FFQ_raw_proteins_aa_proteins'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.proteins;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_essentialAA'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.essentialAA;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_nonEssentialAA'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.nonEssentialAA;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Ile'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Ile;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Leu'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Leu;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Lys'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Lys;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Met'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Met;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Cys'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Cys;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Phe'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Phe;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Tyr'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Tyr;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Thr'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Thr;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Trp'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Trp;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Val'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Val;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Arg'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Arg;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_His'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.His;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Ala'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Ala;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Asp'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Asp;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Glu'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Glu;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Gly'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Gly;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Pro'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Pro;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_proteins_aa_Ser'
            mega_mtrx(iVar, :) = nutrition.raw.proteins_aa.Ser;
            category{iVar} = 'FFQ_raw_proteins_aa';
        case 'FFQ_raw_vitamins_ARE'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.A_retinolEquivalent;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_AR'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.A_retinol;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_ABC'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.A_bCarotene;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_D'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.D_calciferole;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_ETE'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.E_tocopherolEquivalent;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_EAT'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.E_alphaTocopherol;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_K'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.K_phylloquinone;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_B1'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.B1_thiamin;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_B2'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.B2_riboflavin;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_B3_niacin'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.B3_niacin;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_B3_NE'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.B3_niacinEquivalents;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_B5'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.B5_pantothenicA;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_B6'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.B6_pyridoxin;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_B7'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.B7_biotin;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_B9'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.B9_folicA;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_B12'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.B12_cobalamin;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_vitamins_C'
            mega_mtrx(iVar, :) = nutrition.raw.vitamins.C_ascorbicA;
            category{iVar} = 'FFQ_raw_vitamins';
        case 'FFQ_raw_minerals_ashes'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.minerals_ashes;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_Na'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.Na;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_K'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.K;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_Ca'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.Ca;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_Mg'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.Mg;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_P'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.P;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_S'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.S;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_Cl'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.Cl;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_Fe'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.Fe;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_Zn'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.Zn;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_Cu'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.Cu;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_Mn'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.Mn;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_F'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.F;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_I'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.I;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_minerals_NaCl'
            mega_mtrx(iVar, :) = nutrition.raw.minerals.NaCl_salt;
            category{iVar} = 'FFQ_raw_minerals';
        case 'FFQ_raw_carbs_absorbableCarbohydrates'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.absorbableCarbohydrates;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_mannitol'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.mannitol;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_sorbitol'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.sorbitol;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_xylitol'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.xylitol;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_sugarAlcohols'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.sugarAlcohols;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_Glc'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.Glc;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_fructose'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.fructose;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_galactose'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.galactose;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_monosaccharides'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.monosaccharides;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_saccharose'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.saccharose;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_maltose'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.maltose;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_lactose'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.lactose;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_disaccharides'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.disaccharides;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_sugar'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.sugar;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_reabsorbableOligosaccharides'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.reabsorbableOligosaccharides;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_nonReabsorbableOligosaccharides'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.nonReabsorbableOligosaccharides;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_glycogen_animalStarch'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.glycogen_animalStarch;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_starch'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.starch;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_polysaccharides'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.polysaccharides;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_polypentoses'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.polypentoses;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_polyhexoses'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.polyhexoses;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_carbs_polyuronicA'
            mega_mtrx(iVar, :) = nutrition.raw.carbohydrates.polyuronicA;
            category{iVar} = 'FFQ_raw_carbs';
        case 'FFQ_raw_fibers_fibers'
            mega_mtrx(iVar, :) = nutrition.raw.fibers.fibers;
            category{iVar} = 'FFQ_raw_fibers';
        case 'FFQ_raw_fibers_solubleFibers'
            mega_mtrx(iVar, :) = nutrition.raw.fibers.solubleFibers;
            category{iVar} = 'FFQ_raw_fibers';
        case 'FFQ_raw_fibers_insolubleFibers'
            mega_mtrx(iVar, :) = nutrition.raw.fibers.insolubleFibers;
            category{iVar} = 'FFQ_raw_fibers';
        case 'FFQ_raw_fibers_cellulose'
            mega_mtrx(iVar, :) = nutrition.raw.fibers.cellulose;
            category{iVar} = 'FFQ_raw_fibers';
        case 'FFQ_raw_fibers_lignine'
            mega_mtrx(iVar, :) = nutrition.raw.fibers.lignine;
            category{iVar} = 'FFQ_raw_fibers';
        case 'FFQ_raw_fa_fat'
            mega_mtrx(iVar, :) = nutrition.raw.fa.fat;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_monoInsaturated'
            mega_mtrx(iVar, :) = nutrition.raw.fa.fa_monoInsaturated;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_polyInsaturated'
            mega_mtrx(iVar, :) = nutrition.raw.fa.fa_polyInsaturated;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_shortChain'
            mega_mtrx(iVar, :) = nutrition.raw.fa.fa_shortChain;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_mediumChain'
            mega_mtrx(iVar, :) = nutrition.raw.fa.fa_mediumChain;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_longChain'
            mega_mtrx(iVar, :) = nutrition.raw.fa.fa_longChain;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_butanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.butanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_hexanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.hexanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_octanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.octanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_decanoicA_capricA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.decanoicA_capricA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_dodecanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.dodecanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_tetradecanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.tetradecanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_pentadecanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.pentadecanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_hexadecanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.hexadecanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_heptadecanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.heptadecanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_octadecanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.octadecanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_eicosanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.eicosanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_decanoicA_behenicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.decanoicA_behenicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_tetracosanoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.tetracosanoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_saturatedFA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.saturatedFA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_tetradecenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.tetradecenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_pentadecenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.pentadecenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_hexadecenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.hexadecenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_heptadecenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.heptadecenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_octadecenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.octadecenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_eicosenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.eicosenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_decosenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.decosenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_tetracosenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.tetracosenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_hexadecadienoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.hexadecadienoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_hexacedatetraenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.hexacedatetraenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_octadecadienoicA_linoleicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.octadecadienoicA_linoleicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_octadecatrienoicA_linolenicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.octadecatrienoicA_linolenicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_octradecatetraenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.octradecatetraenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_nonAdecatrienoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.nonAdecatrienoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_eicosadienoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.eicosadienoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_eicosatrienoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.eicosatrienoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_eicosatetraenoicA_arachidonicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.eicosatetraenoicA_arachidonicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_eicosapentaenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.eicosapentaenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_docosadienoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.docosadienoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_docosatrienoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.docosatrienoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_docosatetraenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.docosatetraenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_docosapentaenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.docosapentaenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_docosahexaenoicA'
            mega_mtrx(iVar, :) = nutrition.raw.fa.docosahexaenoicA;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_omegas3'
            mega_mtrx(iVar, :) = nutrition.raw.fa.omegas3;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_omegas6'
            mega_mtrx(iVar, :) = nutrition.raw.fa.omegas6;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_glycerine_lipids'
            mega_mtrx(iVar, :) = nutrition.raw.fa.glycerine_lipids;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_cholesterol'
            mega_mtrx(iVar, :) = nutrition.raw.fa.cholesterol;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_fa_GFPS_PolyI_div_Sat_fa_ratio'
            mega_mtrx(iVar, :) = nutrition.raw.fa.GFPS_PolyI_div_Sat_fa_ratio;
            category{iVar} = 'FFQ_raw_fa';
        case 'FFQ_raw_water'
            mega_mtrx(iVar, :) = nutrition.raw.other.water;
            category{iVar} = 'FFQ_raw_water';
        case 'FFQ_raw_organicA'
            mega_mtrx(iVar, :) = nutrition.raw.other.organicAcids;
            category{iVar} = 'FFQ_raw_organicA';
        case 'FFQ_raw_alcohol'
            mega_mtrx(iVar, :) = nutrition.raw.other.alcohol;
            category{iVar} = 'FFQ_raw_alcohol';
        case 'FFQ_raw_uricA'
            mega_mtrx(iVar, :) = nutrition.raw.other.uricA;
            category{iVar} = 'FFQ_raw_uricA';
        case 'FFQ_raw_purin'
            mega_mtrx(iVar, :) = nutrition.raw.other.purin;
            category{iVar} = 'FFQ_raw_purin';
        case 'FFQ_raw_bread'
            mega_mtrx(iVar, :) = nutrition.raw.other.bread;
            category{iVar} = 'FFQ_raw_bread';
        case 'FFQ_norm_energy_kCal'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.energy.kCal;
            category{iVar} = 'FFQ_energy';
        case 'FFQ_norm_energy_kJ'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.energy.kJ;
            category{iVar} = 'FFQ_energy';
        case 'FFQ_norm_energy_kCalWithFibers'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.energy.kCal_withFibers;
            category{iVar} = 'FFQ_energy';
        case 'FFQ_norm_energy_kJ_withFibers'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.energy.kJ_withFibers;
            category{iVar} = 'FFQ_energy';
        case 'FFQ_norm_proteins_aa_proteins'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.proteins;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_essentialAA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.essentialAA;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_nonEssentialAA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.nonEssentialAA;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Ile'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Ile;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Leu'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Leu;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Lys'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Lys;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Met'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Met;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Cys'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Cys;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Phe'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Phe;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Tyr'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Tyr;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Thr'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Thr;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Trp'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Trp;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Val'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Val;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Arg'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Arg;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_His'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.His;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Ala'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Ala;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Asp'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Asp;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Glu'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Glu;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Gly'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Gly;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Pro'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Pro;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_proteins_aa_Ser'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.proteins_aa.Ser;
            category{iVar} = 'FFQ_proteins_aa';
        case 'FFQ_norm_vitamins_ARE'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.A_retinolEquivalent;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_AR'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.A_retinol;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_ABC'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.A_bCarotene;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_D'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.D_calciferole;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_ETE'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.E_tocopherolEquivalent;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_EAT'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.E_alphaTocopherol;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_K'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.K_phylloquinone;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_B1'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.B1_thiamin;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_B2'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.B2_riboflavin;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_B3_niacin'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.B3_niacin;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_B3_NE'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.B3_niacinEquivalents;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_B5'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.B5_pantothenicA;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_B6'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.B6_pyridoxin;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_B7'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.B7_biotin;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_B9'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.B9_folicA;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_B12'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.B12_cobalamin;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_vitamins_C'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.vitamins.C_ascorbicA;
            category{iVar} = 'FFQ_vitamins';
        case 'FFQ_norm_minerals_ashes'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.minerals_ashes;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_Na'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.Na;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_K'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.K;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_Ca'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.Ca;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_Mg'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.Mg;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_P'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.P;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_S'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.S;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_Cl'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.Cl;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_Fe'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.Fe;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_Zn'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.Zn;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_Cu'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.Cu;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_Mn'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.Mn;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_F'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.F;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_I'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.I;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_minerals_NaCl'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.minerals.NaCl_salt;
            category{iVar} = 'FFQ_minerals';
        case 'FFQ_norm_carbs_absorbableCarbohydrates'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.absorbableCarbohydrates;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_mannitol'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.mannitol;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_sorbitol'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.sorbitol;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_xylitol'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.xylitol;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_sugarAlcohols'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.sugarAlcohols;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_Glc'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.Glc;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_fructose'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.fructose;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_galactose'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.galactose;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_monosaccharides'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.monosaccharides;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_saccharose'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.saccharose;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_maltose'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.maltose;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_lactose'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.lactose;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_disaccharides'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.disaccharides;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_sugar'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.sugar;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_reabsorbableOligosaccharides'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.reabsorbableOligosaccharides;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_nonReabsorbableOligosaccharides'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.nonReabsorbableOligosaccharides;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_glycogen_animalStarch'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.glycogen_animalStarch;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_starch'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.starch;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_polysaccharides'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.polysaccharides;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_polypentoses'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.polypentoses;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_polyhexoses'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.polyhexoses;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_carbs_polyuronicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.carbohydrates.polyuronicA;
            category{iVar} = 'FFQ_carbs';
        case 'FFQ_norm_fibers_fibers'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fibers.fibers;
            category{iVar} = 'FFQ_fibers';
        case 'FFQ_norm_fibers_solubleFibers'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fibers.solubleFibers;
            category{iVar} = 'FFQ_fibers';
        case 'FFQ_norm_fibers_insolubleFibers'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fibers.insolubleFibers;
            category{iVar} = 'FFQ_fibers';
        case 'FFQ_norm_fibers_cellulose'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fibers.cellulose;
            category{iVar} = 'FFQ_fibers';
        case 'FFQ_norm_fibers_lignine'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fibers.lignine;
            category{iVar} = 'FFQ_fibers';
        case 'FFQ_norm_fa_fat'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.fat;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_monoInsaturated'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.fa_monoInsaturated;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_polyInsaturated'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.fa_polyInsaturated;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_shortChain'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.fa_shortChain;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_mediumChain'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.fa_mediumChain;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_longChain'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.fa_longChain;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_butanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.butanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_hexanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.hexanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_octanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.octanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_decanoicA_capricA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.decanoicA_capricA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_dodecanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.dodecanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_tetradecanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.tetradecanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_pentadecanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.pentadecanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_hexadecanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.hexadecanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_heptadecanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.heptadecanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_octadecanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.octadecanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_eicosanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.eicosanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_decanoicA_behenicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.decanoicA_behenicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_tetracosanoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.tetracosanoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_saturatedFA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.saturatedFA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_tetradecenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.tetradecenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_pentadecenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.pentadecenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_hexadecenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.hexadecenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_heptadecenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.heptadecenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_octadecenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.octadecenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_eicosenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.eicosenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_decosenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.decosenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_tetracosenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.tetracosenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_hexadecadienoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.hexadecadienoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_hexacedatetraenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.hexacedatetraenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_octadecadienoicA_linoleicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.octadecadienoicA_linoleicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_octadecatrienoicA_linolenicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.octadecatrienoicA_linolenicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_octradecatetraenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.octradecatetraenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_nonAdecatrienoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.nonAdecatrienoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_eicosadienoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.eicosadienoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_eicosatrienoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.eicosatrienoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_eicosatetraenoicA_arachidonicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.eicosatetraenoicA_arachidonicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_eicosapentaenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.eicosapentaenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_docosadienoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.docosadienoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_docosatrienoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.docosatrienoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_docosatetraenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.docosatetraenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_docosapentaenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.docosapentaenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_docosahexaenoicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.docosahexaenoicA;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_omegas3'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.omegas3;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_omegas6'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.omegas6;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_glycerine_lipids'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.glycerine_lipids;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_cholesterol'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.cholesterol;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_fa_GFPS_PolyI_div_Sat_fa_ratio'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.fa.GFPS_PolyI_div_Sat_fa_ratio;
            category{iVar} = 'FFQ_fa';
        case 'FFQ_norm_water'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.other.water;
            category{iVar} = 'FFQ_water';
        case 'FFQ_norm_organicA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.other.organicAcids;
            category{iVar} = 'FFQ_organicA';
        case 'FFQ_norm_alcohol'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.other.alcohol;
            category{iVar} = 'FFQ_alcohol';
        case 'FFQ_norm_uricA'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.other.uricA;
            category{iVar} = 'FFQ_uricA';
        case 'FFQ_norm_purin'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.other.purin;
            category{iVar} = 'FFQ_purin';
        case 'FFQ_norm_bread'
            mega_mtrx(iVar, :) = nutrition.norm_by_totalCal.other.bread;
            category{iVar} = 'FFQ_bread';
        case 'behavior_stress_ratings_S1'
            mega_mtrx(iVar, :) = preMRS_stress;
            category{iVar} = 'stress_anxiety';
        case 'behavior_stress_ratings_S2'
            mega_mtrx(iVar, :) = postMRS_stress;
            category{iVar} = 'stress_anxiety';
        case 'behavior_stress_ratings_S3'
            mega_mtrx(iVar, :) = prefMRI_stress;
            category{iVar} = 'stress_anxiety';
        case 'behavior_stress_ratings_S4'
            mega_mtrx(iVar, :) = postfMRI_stress;
            category{iVar} = 'stress_anxiety';
        case 'behavior_stress_ratings_S4_min_S1'
            mega_mtrx(iVar, :) = deltaStressPrePostExp;
            category{iVar} = 'stress_anxiety';
        case 'behavior_stress_ratings_S4_min_S3'
            mega_mtrx(iVar, :) = deltaStressPrePostfMRI;
            category{iVar} = 'stress_anxiety';
        case 'behavior_fatigue_ratings_F1'
            mega_mtrx(iVar, :) = preMRS_fatigue;
            category{iVar} = 'fatigue';
        case 'behavior_fatigue_ratings_F2'
            mega_mtrx(iVar, :) = postMRS_fatigue;
            category{iVar} = 'fatigue';
        case 'behavior_fatigue_ratings_F3'
            mega_mtrx(iVar, :) = prefMRI_fatigue;
            category{iVar} = 'fatigue';
        case 'behavior_fatigue_ratings_F4'
            mega_mtrx(iVar, :) = postfMRI_fatigue;
            category{iVar} = 'fatigue';
        case 'behavior_fatigue_ratings_F4_min_F1'
            mega_mtrx(iVar, :) = deltaFatiguePrePostExp;
            category{iVar} = 'fatigue';
        case 'behavior_fatigue_ratings_F4_min_F3'
            mega_mtrx(iVar, :) = deltaFatiguePrePostfMRI;
            category{iVar} = 'fatigue';
        case 'behavior_task_choices_HE'
            mega_mtrx(iVar, :) = choice_hE.EpEm;
            category{iVar} = 'motivation';
        case 'behavior_task_choices_HPE'
            mega_mtrx(iVar, :) = choice_hE.Ep;
            category{iVar} = 'motivation';
        case 'behavior_task_choices_HME'
            mega_mtrx(iVar, :) = choice_hE.Em;
            category{iVar} = 'motivation';
        case 'behavior_task_prm_kR'
            mega_mtrx(iVar, :) = prm.kR;
            category{iVar} = 'motivation';
        case 'behavior_task_prm_kP'
            mega_mtrx(iVar, :) = prm.kP;
            category{iVar} = 'motivation';
        case 'behavior_task_prm_kEp'
            mega_mtrx(iVar, :) = prm.kEp;
            category{iVar} = 'motivation';
        case 'behavior_task_prm_kEm'
            mega_mtrx(iVar, :) = prm.kEm;
            category{iVar} = 'motivation';
        case 'behavior_task_prm_kFp'
            mega_mtrx(iVar, :) = prm.kFp;
            category{iVar} = 'motivation';
        case 'behavior_task_prm_kLm'
            mega_mtrx(iVar, :) = prm.kLm;
            category{iVar} = 'motivation';
        case 'behavior_task_prm_kBias'
            mega_mtrx(iVar, :) = prm.kBias;
            category{iVar} = 'motivation';
    end
end % variable loop

%% create alternative category variable with very short name and even more general
category_bis = category;
category_bis = strrep(category_bis,'salivary_Testosterone','St');
category_bis = strrep(category_bis,'salivary_Cortisol','Sc');
category_bis = strrep(category_bis,'salivary_Testo_div_Cort','Stc');
category_bis = strrep(category_bis,'salivary_IL','Sil');
category_bis = strrep(category_bis,'wholeB_NADomics','WBn');
category_bis = strrep(category_bis,'plasma_aa','Paa');
category_bis = strrep(category_bis,'plasma_Lac','Plac');
category_bis = strrep(category_bis,'plasma_fa','Pfa');
category_bis = strrep(category_bis,'dmPFC_aa','D');
category_bis = strrep(category_bis,'dmPFC_mb','D');
category_bis = strrep(category_bis,'aIns_aa','A');
category_bis = strrep(category_bis,'aIns_mb','A');
category_bis = strrep(category_bis,'FFQ_raw_energy','FFQre');
category_bis = strrep(category_bis,'FFQ_raw_proteins_aa','FFQraa');
category_bis = strrep(category_bis,'FFQ_raw_vitamins','FFQrv');
category_bis = strrep(category_bis,'FFQ_raw_minerals','FFQrm');
category_bis = strrep(category_bis,'FFQ_raw_carbs','FFQrc');
category_bis = strrep(category_bis,'FFQ_raw_fibers','FFQrf');
category_bis = strrep(category_bis,'FFQ_raw_fa','FFQrfa');
category_bis = strrep(category_bis,'FFQ_energy','FFQe');
category_bis = strrep(category_bis,'FFQ_proteins_aa','FFQaa');
category_bis = strrep(category_bis,'FFQ_vitamins','FFQv');
category_bis = strrep(category_bis,'FFQ_minerals','FFQm');
category_bis = strrep(category_bis,'FFQ_carbs','FFQc');
category_bis = strrep(category_bis,'FFQ_fibers','FFQf');
category_bis = strrep(category_bis,'FFQ_fa','FFQfa');
category_bis = strrep(category_bis,'stress_anxiety','Bs');
category_bis = strrep(category_bis,'motivation','Bm');
category_bis = strrep(category_bis,'fatigue','Bf');
category_bis = strrep(category_bis,'dominance','Bd');

%% perform correlation matrix between all variables
[corr_mtrx, pval_mtrx] = deal(NaN(n_mega_mtrx_vars, n_mega_mtrx_vars));
[mtrx_var_nm1, mtrx_var_nm2,...
    mtrx_categ_nm1, mtrx_categ_nm2] =  deal(cell(n_mega_mtrx_vars, n_mega_mtrx_vars));
for iVar = 1:n_mega_mtrx_vars % loop through rows (all variables)
    for jVar = 1:n_mega_mtrx_vars % loop through columns (all variables)
        okSubs = ~isnan(mega_mtrx(iVar,:).*mega_mtrx(jVar,:));
        
        % extract variable name
        mtrx_var_nm1{iVar, jVar} = mega_mtrx_names{iVar};
        mtrx_var_nm2{iVar, jVar} = mega_mtrx_names{jVar};
        % extract corresponding category of variable
        mtrx_categ_nm1{iVar, jVar} = category{iVar};
        mtrx_categ_nm2{iVar, jVar} = category{jVar};
        
        % perform correlation
        if sum(okSubs) > 0 % at least some subjects ok => matlab can perform the correlation
            [corr_mtrx(iVar, jVar),...
                pval_mtrx(iVar, jVar)] = corr(mega_mtrx(iVar,okSubs)', mega_mtrx(jVar,okSubs)');
        else % no good subjects
            % manually force correlation at 0, otherwise Matlab puts a -1
            % correlation through imagesc when given NaN as input which is
            % completely misleading
            corr_mtrx(iVar, jVar) = 0;
            pval_mtrx(iVar, jVar) = 1;
        end
    end % variable loop
end % variable loop

%% extract corrected p.value for the number of tests
% first need to reshape the p.values
pval_linearized = reshape(pval_mtrx,[n_mega_mtrx_vars*n_mega_mtrx_vars,1]);
% perform the correction
pval_corr_linear = pval_adjust(pval_linearized, 'bonferroni');
% reshape matrix to have the same shape as for the other matrices
% However, note that the correction here is overly stringent as all
% correlations are performed twice
pval_mtrx_Bonferroni_corr = reshape(pval_corr_linear,[n_mega_mtrx_vars,n_mega_mtrx_vars]);

%% perform GLM including age and sex
[GLM_b_mtrx_age, GLM_pval_mtrx_age,...
    GLM_b_mtrx_sex, GLM_pval_mtrx_sex,...
    GLM_b_mtrx_vars, GLM_pval_mtrx_vars] = deal(NaN(n_mega_mtrx_vars, n_mega_mtrx_vars));
% extract age and sex for each subject (co-variates)
age = nanzscore(mega_mtrx(strcmp(mega_mtrx_names,'gal_age'),:))'; % zscore age to make beta comparable
binary_var_sex = mega_mtrx(strcmp(mega_mtrx_names,'gal_sex'),:)';

% perform LME test on each variable with each variable as requested...
for iVar = 3:n_mega_mtrx_vars % loop through explained variables (except age and sex)
    explained_var_nm = mega_mtrx_names{iVar};
    for jVar = 3:n_mega_mtrx_vars % loop through explaining variables (except age and sex)
        explaining_var_nm = mega_mtrx_names{jVar};
        
        % name for storage
        curr_test_nm = [explained_var_nm,'_f_',explaining_var_nm];
        
        % zscore input and output to make betas comparable (unless data is
        % binary)
        % explaining variable
        switch explaining_var_nm
            case {'gal_sex'} % binary variables (no zscore)
                explaining_var = mega_mtrx(jVar,:)';
            otherwise % continuous variables (zscore)
                explaining_var = nanzscore(mega_mtrx(jVar,:))';
        end
        % explained variable
        switch explained_var_nm
            case {'gal_sex'} % binary variables (no zscore)
                explained_var = mega_mtrx(iVar,:)';
            otherwise % continuous variables (zscore)
                explained_var = nanzscore(mega_mtrx(iVar,:))';
        end
        % perform GLM
        [betas.(curr_test_nm),~,stats_tmp] = glmfit([age, binary_var_sex, explaining_var], explained_var,'normal');
        pval.(curr_test_nm).age = stats_tmp.p(2);
        pval.(curr_test_nm).sex = stats_tmp.p(3);
        pval.(curr_test_nm).(explaining_var_nm) = stats_tmp.p(4);
        
        % store data in big matrices (with explained variable varying in
        % the Y.axis (i.e. rows) and explaining variable varying in the
        % X.axis (i.e. columns)
        % age
        GLM_b_mtrx_age(iVar, jVar) = betas.(curr_test_nm)(2);
        GLM_pval_mtrx_age(iVar, jVar) = stats_tmp.p(2);
        % sex
        GLM_b_mtrx_sex(iVar, jVar) = betas.(curr_test_nm)(3);
        GLM_pval_mtrx_sex(iVar, jVar) = stats_tmp.p(3);
        % variable
        GLM_b_mtrx_vars(iVar, jVar) = betas.(curr_test_nm)(4);
        GLM_pval_mtrx_vars(iVar, jVar) = stats_tmp.p(4);
    end % loop through explaining variables
end % loop through explained variable


%% extract smaller correlation matrices

% select relevant categories
circulatory_mb_idx = ismember(category,{'salivary_Testosterone','salivary_Cortisol',...
    'salivary_Testo_div_Cort','salivary_IL','wholeB_NADomics',...
    'plasma_aa','plasma_Lac','plasma_fa'});
brain_mb_idx = ismember(category,{'dmPFC_aa','dmPFC_mb','aIns_aa','aIns_mb'});
bhv_idx = ismember(category,{'stress_anxiety','motivation','fatigue','dominance'});
FFQ_raw_idx = ismember(category,{'FFQ_raw_energy','FFQ_raw_proteins_aa',...
    'FFQ_raw_vitamins','FFQ_raw_minerals','FFQ_raw_carbs','FFQ_raw_fibers',...
    'FFQ_raw_fa'});
FFQ_norm_idx = ismember(category,{'FFQ_energy','FFQ_proteins_aa',...
    'FFQ_vitamins','FFQ_minerals','FFQ_carbs','FFQ_fibers',...
    'FFQ_fa'});
% combinations
circulatory_brain_mb_idx = circulatory_mb_idx | brain_mb_idx;
% circulatory_brain_mb_FFQ_raw_idx = circulatory_brain_mb_idx | FFQ_raw_idx;
% circulatory_brain_mb_FFQ_norm_idx = circulatory_brain_mb_idx | FFQ_norm_idx;
% circulatory_brain_mb_bhv_idx = circulatory_brain_mb_idx | bhv_idx;
circulatory_mb_bhv_idx = circulatory_mb_idx | bhv_idx;
brain_mb_bhv_idx = brain_mb_idx | bhv_idx;

% brain metabolites + circulatory metabolites (whole-blood + plasma + salivary)
n_circulatory_brain_mb = sum(circulatory_brain_mb_idx);
[brain_mb_circulatory_mb_corr_mtrx,...
    brain_mb_circulatory_mb_pval_mtrx] = deal(NaN(n_circulatory_brain_mb));
brain_mb_circulatory_mb_corr_mtrx(:,:) = corr_mtrx(circulatory_brain_mb_idx, circulatory_brain_mb_idx);
brain_mb_circulatory_mb_pval_mtrx(:,:) = pval_mtrx(circulatory_brain_mb_idx, circulatory_brain_mb_idx);
brain_mb_circulatory_mb_names = mega_mtrx_names(circulatory_brain_mb_idx);
% create label for circular plot
brain_mb_circulatory_mb_categ = category_bis(circulatory_brain_mb_idx);

% behavior + circulatory metabolites (whole-blood + plasma + salivary))
n_circulatory_mb_bhv = sum(circulatory_mb_bhv_idx);
[circulatory_mb_bhv_corr_mtrx,...
    circulatory_mb_bhv_pval_mtrx] = deal(NaN(n_circulatory_mb_bhv));
circulatory_mb_bhv_corr_mtrx(:,:) = corr_mtrx(circulatory_mb_bhv_idx, circulatory_mb_bhv_idx);
circulatory_mb_bhv_pval_mtrx(:,:) = pval_mtrx(circulatory_mb_bhv_idx, circulatory_mb_bhv_idx);
circulatory_mb_bhv_names = mega_mtrx_names(circulatory_mb_bhv_idx);
% create label for circular plot
circulatory_mb_bhv_categ = category_bis(circulatory_mb_bhv_idx);

% behavior + brain metabolites
n_brain_mb_bhv = sum(brain_mb_bhv_idx);
[brain_mb_bhv_corr_mtrx,...
    brain_mb_bhv_pval_mtrx] = deal(NaN(n_brain_mb_bhv));
brain_mb_bhv_corr_mtrx(:,:) = corr_mtrx(brain_mb_bhv_idx, brain_mb_bhv_idx);
brain_mb_bhv_pval_mtrx(:,:) = pval_mtrx(brain_mb_bhv_idx, brain_mb_bhv_idx);
brain_mb_bhv_names = mega_mtrx_names(brain_mb_bhv_idx);
% create label for circular plot
brain_mb_bhv_categ = category_bis(brain_mb_bhv_idx);

%% save all the data to avoid relaunching everything every single time
save([saveFolder,filesep,'mega_correlation_mtrx_data.mat']);

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
corr_range = [-1 1];

%% correlation matrix
% display a first correlation matrix with all tests
apply_pval_threshold = false;
pval_threshold = [];
disp_signif_stars = false;
corr_plot(corr_mtrx, pval_mtrx,...
    corr_range, [], [], [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);

% display a correlation matrix filtering based on the corrected p.values
apply_pval_threshold = true;
pval_threshold = 0.05;
disp_signif_stars = false;
corr_plot(corr_mtrx, pval_mtrx_Bonferroni_corr,...
    corr_range, [], [], [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);

%% GLM results after correcting for sex and age
GLM_range = [-2 2];
% show two figures, one with all tests and one with filtering on p.value
for iFig = 1:2
    switch iFig
        case 1 % no filter
            apply_pval_threshold = false;
            pval_threshold = [];
        case 2 % filter on p.value
            apply_pval_threshold = true;
            pval_threshold = 0.05;
    end
    disp_signif_stars = false; % do not show stars
    fig_correl = fig;
    
    % age variable
%     GLM_age_range = [min(min(GLM_b_mtrx_age,[],'omitnan'),[],'omitnan'),...
%         max(max(GLM_b_mtrx_age,[],'omitnan'),[],'omitnan')];
    subplot_age_hdl = subplot(1,3,1);
    title('Age');
    corr_plot(GLM_b_mtrx_age, GLM_pval_mtrx_age,...
        GLM_range, [], [], [], [],...
        apply_pval_threshold, pval_threshold, disp_signif_stars,...
        fig_correl, subplot_age_hdl);
    
    % sex variable
%     GLM_sex_range = [min(min(GLM_b_mtrx_sex,[],'omitnan'),[],'omitnan'),...
%         max(max(GLM_b_mtrx_sex,[],'omitnan'),[],'omitnan')];
    subplot_sex_hdl = subplot(1,3,2);
    title('Sex');
    corr_plot(GLM_b_mtrx_sex, GLM_pval_mtrx_sex,...
        GLM_range, [], [], [], [],...
        apply_pval_threshold, pval_threshold, disp_signif_stars,...
        fig_correl, subplot_sex_hdl);
    
    % variable correlation
%     GLM_var_range = [min(min(GLM_b_mtrx_vars,[],'omitnan'),[],'omitnan'),...
%         max(max(GLM_b_mtrx_vars,[],'omitnan'),[],'omitnan')];
    subplot_var_hdl = subplot(1,3,3);
    title('Variables of interest');
    corr_plot(GLM_b_mtrx_vars, GLM_pval_mtrx_vars,...
        GLM_range, [], [], [], [],...
        apply_pval_threshold, pval_threshold, disp_signif_stars,...
        fig_correl, subplot_var_hdl);
end % loop through figures (with/without filtering on p.value)

%% perform circular correlation matrices between specific parts
pval_thresh = 0.001;
% brain metabolites + circulatory metabolites (whole-blood + plasma +
% salivary)
circ_disp(brain_mb_circulatory_mb_corr_mtrx,...
    brain_mb_circulatory_mb_pval_mtrx, pval_thresh, brain_mb_circulatory_mb_categ);

% behavior + circulatory metabolites (whole-blood + plasma + salivary)
circ_disp(circulatory_mb_bhv_corr_mtrx,...
    circulatory_mb_bhv_pval_mtrx, pval_thresh, circulatory_mb_bhv_categ);

% behavior + brain metabolites
circ_disp(brain_mb_bhv_corr_mtrx,...
    brain_mb_bhv_pval_mtrx, pval_thresh, brain_mb_bhv_categ);

%% save excel file for Gephi
% create one table for nodes
table_nodes = table;
table_nodes.Id = mega_mtrx_names';
table_nodes.Label = mega_mtrx_names';
% table_nodes.interval = repmat([],n_mega_mtrx_vars,1);
% save node tables
writetable(table_nodes,[saveFolder,filesep,'Nodes_megamtrx_',num2str(NS),'subs.xlsx']);

% create one table for edges (ie correlations)
n_corrs = n_mega_mtrx_vars*n_mega_mtrx_vars;
table_edges = table;
table_edges.Source = reshape(mtrx_var_nm1,[n_corrs,1]);
table_edges.Target = reshape(mtrx_var_nm2,[n_corrs,1]);
table_edges.weight = reshape(abs(corr_mtrx),[n_corrs,1]); % use |r| as visualization software Gephi does not deal with negative values 
table_edges.weight2 = reshape(pval_mtrx,[n_corrs,1]);
table_edges.Type = repmat('undirected',n_corrs,1);
% save edges table
writetable(table_edges,[saveFolder,filesep,'Edges_megamtrx_',num2str(NS),'subs.xlsx']);

% create one table for Nestlé with main information of interest for
% correlations (using absolute correlation coefficient as some softwares
% cannot process negative values)
n_corrs = n_mega_mtrx_vars*n_mega_mtrx_vars;
table_correl = table;
% variable names
table_correl.var1 = reshape(mtrx_var_nm1,[n_corrs,1]);
table_correl.category_var1 = reshape(mtrx_categ_nm1,[n_corrs,1]);
table_correl.var2 = reshape(mtrx_var_nm2,[n_corrs,1]);
table_correl.category_var2 = reshape(mtrx_categ_nm2,[n_corrs,1]);
% regressors
table_correl.r_corr = reshape(abs(corr_mtrx),[n_corrs,1]); % use |r|
table_correl.pvalue = reshape(pval_mtrx,[n_corrs,1]);
% save edges table
writetable(table_correl,[saveFolder,filesep,'crosscorrel_abs_r_table_',num2str(NS),'subs.xlsx']);

% create one table for Nestlé with main information of interest for
% correlations (using signed correlation coefficients)
n_corrs = n_mega_mtrx_vars*n_mega_mtrx_vars;
table_correl = table;
% variable names
table_correl.var1 = reshape(mtrx_var_nm1,[n_corrs,1]);
table_correl.category_var1 = reshape(mtrx_categ_nm1,[n_corrs,1]);
table_correl.var2 = reshape(mtrx_var_nm2,[n_corrs,1]);
table_correl.category_var2 = reshape(mtrx_categ_nm2,[n_corrs,1]);
% regressors
table_correl.r_corr = reshape(corr_mtrx,[n_corrs,1]); % use r
table_correl.pvalue = reshape(pval_mtrx,[n_corrs,1]);
% save edges table
writetable(table_correl,[saveFolder,filesep,'crosscorrel_signed_r_table_',num2str(NS),'subs.xlsx']);

% create second table for Nestlé with main information of interest for
% GLM
n_corrs = n_mega_mtrx_vars*n_mega_mtrx_vars;
table_GLM = table;
% variable names
table_GLM.var1 = reshape(mtrx_var_nm1,[n_corrs,1]);
table_GLM.category_var1 = reshape(mtrx_categ_nm1,[n_corrs,1]);
table_GLM.var2 = reshape(mtrx_var_nm2,[n_corrs,1]);
table_GLM.category_var2 = reshape(mtrx_categ_nm2,[n_corrs,1]);
% regressors
table_GLM.b_age = reshape(abs(GLM_b_mtrx_age),[n_corrs,1]);
table_GLM.pvalue_age = reshape(GLM_pval_mtrx_age,[n_corrs,1]);
table_GLM.b_sex = reshape(abs(GLM_b_mtrx_sex),[n_corrs,1]);
table_GLM.pvalue_sex = reshape(GLM_pval_mtrx_sex,[n_corrs,1]);
table_GLM.b_var = reshape(abs(GLM_b_mtrx_vars),[n_corrs,1]);
table_GLM.pvalue_var = reshape(GLM_pval_mtrx_vars,[n_corrs,1]);
% save edges table
writetable(table_GLM,[saveFolder,filesep,'GLM_table_abs_betas_',num2str(NS),'subs.xlsx']);

% create second table for Nestlé with main information of interest for
% GLM
n_corrs = n_mega_mtrx_vars*n_mega_mtrx_vars;
table_GLM = table;
% variable names
table_GLM.var1 = reshape(mtrx_var_nm1,[n_corrs,1]);
table_GLM.category_var1 = reshape(mtrx_categ_nm1,[n_corrs,1]);
table_GLM.var2 = reshape(mtrx_var_nm2,[n_corrs,1]);
table_GLM.category_var2 = reshape(mtrx_categ_nm2,[n_corrs,1]);
% regressors
table_GLM.b_age = reshape(GLM_b_mtrx_age,[n_corrs,1]);
table_GLM.pvalue_age = reshape(GLM_pval_mtrx_age,[n_corrs,1]);
table_GLM.b_sex = reshape(GLM_b_mtrx_sex,[n_corrs,1]);
table_GLM.pvalue_sex = reshape(GLM_pval_mtrx_sex,[n_corrs,1]);
table_GLM.b_var = reshape(GLM_b_mtrx_vars,[n_corrs,1]);
table_GLM.pvalue_var = reshape(GLM_pval_mtrx_vars,[n_corrs,1]);
% save edges table
writetable(table_GLM,[saveFolder,filesep,'GLM_table_signed_betas_',num2str(NS),'subs.xlsx']);

%% selection of a few variables of interest and show corresponding correlation matrix
vars_of_interest = listdlg('PromptString','Select the vars of interest',...
    'ListString',mega_mtrx_names);
n_vars_of_interest = length(vars_of_interest);
corr_selecta = corr_mtrx(vars_of_interest, :);
pval_corr_selecta = pval_mtrx(vars_of_interest, :);
for iV = 1:n_vars_of_interest
    signif_corr_selecta.(mega_mtrx_names{vars_of_interest(iV)}) = mega_mtrx_names(pval_corr_selecta(iV,:) < 0.05);
end

apply_pval_threshold = false;
pval_threshold = [];
disp_signif_stars = false;
corr_plot(corr_selecta, pval_corr_selecta,...
    corr_range, [], mega_mtrx_names(vars_of_interest), [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);