function[nutrition, nutri_categories_names, n_categ, nutri_names, n_nutri] = load_FFQ(study_nm, subject_id, NS)
% [nutrition, nutri_categories_names, n_categ, nutri_names, n_nutri] = load_FFQ(study_nm, subject_id, NS)
% load_FFQ will extract the proportion of nutrients derived from the FFQ
% questionnaire with and without normalization by total caloric intake.
% Nutrients are classified by nutrient type to facilitate further
% processing (glucids, fat, fibers, proteins/amino-acids, etc.).
%
% INPUTS
% study_nm: study name ('study1' by default)
%
% subject_id: list of subjects
%
% NS: number of subjects
%
% OUTPUTS
% nutrition: structure containing all the nutrient values with and without
% caloric normalization
%
% nutri_categories_names: nutrition categories
%
% n_categ: number of categories
%
% nutri_names: cell with list of nutrient names
%
% n_nutri: number of nutrients considered

%% define inputs
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
if ~exist('subject_id','var') || isempty(subject_id) || ~exist('NS','var') || isempty(NS)
    condition = subject_condition;
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
end

%% working directory
root = LGCM_root_paths;
switch root
    case ['E:',filesep]
        gitPath = fullfile('C:','Users','clairis','Desktop');
    case {[filesep,filesep,fullfile('sv-nas1.rcp.epfl.ch','Sandi-lab',...
            'human_data_private','raw_data_subject'),filesep],...
            [fullfile('L:','human_data_private','raw_data_subject'),filesep]}
        gitPath = fullfile('C:','Users','clairis','Documents');
    otherwise
        error('case not ready yet');
end
nutritionPath = [fullfile(gitPath,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'nutrition'),filesep];

%% initialize variables of interest
% normalization type
nutri_norm_names = {'raw','norm_by_totalCal'};
n_nutri_norm = length(nutri_norm_names);
% nutrition categories
nutri_categories_names = {'energy','proteins_aa','vitamins',...
    'minerals','carbohydrates','fibers','fa','other'};
n_categ = length(nutri_categories_names);
% nutrition names for each category
nutri_names.energy = {'kCal','kJ','kCal_withFibers','kJ_withFibers'};
n_nutri.energy = length(nutri_names.energy);
nutri_names.proteins_aa = {'proteins','essentialAA','nonEssentialAA',...
    'Ile','Leu','Lys','Met','Cys','Phe','Tyr','Thr','Trp','Val','Arg',...
    'His','Ala','Asp','Glu','Gly','Pro','Ser'};
n_nutri.proteins_aa = length(nutri_names.proteins_aa);
nutri_names.vitamins = {'A_retinolEquivalent','A_retinol','A_bCarotene',...
    'D_calciferole','E_tocopherolEquivalent','E_alphaTocopherol',...
    'K_phylloquinone','B1_thiamin','B2_riboflavin','B3_niacin','B3_niacinEquivalents',...
    'B5_pantothenicA','B6_pyridoxin','B7_biotin','B9_folicA','B12_cobalamin',...
    'C_ascorbicA'};
n_nutri.vitamins = length(nutri_names.vitamins);
nutri_names.minerals = {'minerals_ashes','Na','K','Ca','Mg','P','S','Cl',...
    'Fe','Zn','Cu','Mn','F','I','NaCl_salt'};
n_nutri.minerals = length(nutri_names.minerals);
nutri_names.carbohydrates = {'absorbableCarbohydrates','mannitol','sorbitol',...
    'xylitol','sugarAlcohols','Glc','fructose','galactose',...
    'monosaccharides','saccharose','maltose','lactose','disaccharides',...
    'sugar','reabsorbableOligosaccharides','nonReabsorbableOligosaccharides',...
    'glycogen_animalStarch','starch','polysaccharides','polypentoses',...
    'polyhexoses','polyuronicA'};
n_nutri.carbohydrates = length(nutri_names.carbohydrates);
nutri_names.fibers = {'fibers','solubleFibers','insolubleFibers',...
    'cellulose','lignine'};
n_nutri.fibers = length(nutri_names.fibers);
nutri_names.fa = {'fat','fa_monoInsaturated','fa_polyInsaturated',...
    'fa_shortChain','fa_mediumChain','fa_longChain',...
    'butanoicA','hexanoicA','octanoicA','decanoicA_capricA','dodecanoicA',...
    'tetradecanoicA','pentadecanoicA','hexadecanoicA','heptadecanoicA',...
    'octadecanoicA','eicosanoicA','decanoicA_behenicA','tetracosanoicA',...
    'saturatedFA','tetradecenoicA','pentadecenoicA','hexadecenoicA',...
    'heptadecenoicA','octadecenoicA','eicosenoicA','decosenoicA',...
    'tetracosenoicA','hexadecadienoicA','hexacedatetraenoicA',...
    'octadecadienoicA_linoleicA','octadecatrienoicA_linolenicA',...
    'octradecatetraenoicA','nonAdecatrienoicA','eicosadienoicA',...
    'eicosatrienoicA','eicosatetraenoicA_arachidonicA','eicosapentaenoicA',...
    'docosadienoicA','docosatrienoicA','docosatetraenoicA','docosapentaenoicA',...
    'docosahexaenoicA','omegas3','omegas6','glycerine_lipids','cholesterol',...
    'GFPS_PolyI_div_Sat_fa_ratio'};
n_nutri.fa = length(nutri_names.fa);
nutri_names.other = {'water','organicAcids','alcohol','uricA','purin','bread'};
n_nutri.other = length(nutri_names.other);

for iNutriNorm = 1:n_nutri_norm
    nutri_norm_nm = nutri_norm_names{iNutriNorm};
    for iCat = 1:n_categ
        categ_nm = nutri_categories_names{iCat};
        for iNutrient = 1:n_nutri.(categ_nm)
            nutri_nm = nutri_names.(categ_nm){iNutrient};
            nutrition.(nutri_norm_nm).(categ_nm).(nutri_nm) = NaN(1,NS);
        end
    end % category
end % nutrition norm

%% load data
FFQ_struct = load([nutritionPath,'FFQ_metabolites_extracted_results.mat']);

%% extract relevant data
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = strcmp(sub_nm, FFQ_struct.subject_id);
    
    % energy
    nutrition.raw.energy.kCal(iS) = FFQ_struct.FFQ_results.calories_kcalPerWeek(sub_idx); % (GCAL)
    nutrition.raw.energy.kJ(iS) = FFQ_struct.FFQ_results.GJEnergie_kJ_portion__PerWeek(sub_idx); % (GJ)
    nutrition.raw.energy.kCal_withFibers(iS) = FFQ_struct.FFQ_results.GCALZBEnergie_kcalAvecFibres_portion__PerWeek(sub_idx); % (GCALZB)
    nutrition.raw.energy.kJ_withFibers(iS) = FFQ_struct.FFQ_results.GJZBEnergie_kJAvecFibres_portion__PerWeek(sub_idx); % (GJZB)
    
    % proteins/amino-acids
    nutrition.raw.proteins_aa.proteins(iS) = FFQ_struct.FFQ_results.proteins_mgPerWeek(sub_idx); % (ZE)
    nutrition.raw.proteins_aa.essentialAA(iS) = FFQ_struct.FFQ_results.EEA_mg_portion__PerWeek(sub_idx); % (EEA)
    nutrition.raw.proteins_aa.nonEssentialAA(iS) = FFQ_struct.FFQ_results.ENA_mg_portion__PerWeek(sub_idx); % (ENA)
    nutrition.raw.proteins_aa.Ile(iS) = FFQ_struct.FFQ_results.EILE_mg_portion__PerWeek(sub_idx); % (EILE)
    nutrition.raw.proteins_aa.Leu(iS) = FFQ_struct.FFQ_results.ELEU_mg_portion__PerWeek(sub_idx); % (ELEU)
    nutrition.raw.proteins_aa.Lys(iS) = FFQ_struct.FFQ_results.ELYS_mg_portion__PerWeek(sub_idx); % (ELYS)
    nutrition.raw.proteins_aa.Met(iS) = FFQ_struct.FFQ_results.EMET_mg_portion__PerWeek(sub_idx); % (EMET)
    nutrition.raw.proteins_aa.Cys(iS) = FFQ_struct.FFQ_results.Cys_mgPerWeek(sub_idx); % (ECYS)
    nutrition.raw.proteins_aa.Phe(iS) = FFQ_struct.FFQ_results.EPHE_mg_portion__PerWeek(sub_idx); % (EPHE)
    nutrition.raw.proteins_aa.Tyr(iS) = FFQ_struct.FFQ_results.ETYR_mg_portion__PerWeek(sub_idx); % (ETYR)
    nutrition.raw.proteins_aa.Thr(iS) = FFQ_struct.FFQ_results.ETHR_mg_portion__PerWeek(sub_idx); % (ETHR)
    nutrition.raw.proteins_aa.Trp(iS) = FFQ_struct.FFQ_results.Trp_mgPerWeek(sub_idx); % (ETRP)
    nutrition.raw.proteins_aa.Val(iS) = FFQ_struct.FFQ_results.EVAL_mg_portion__PerWeek(sub_idx); % (EVAL)
    nutrition.raw.proteins_aa.Arg(iS) = FFQ_struct.FFQ_results.EARG_mg_portion__PerWeek(sub_idx); % (EARG)
    nutrition.raw.proteins_aa.His(iS) = FFQ_struct.FFQ_results.EHIS_mg_portion__PerWeek(sub_idx); % (EHIS)
    nutrition.raw.proteins_aa.Ala(iS) = FFQ_struct.FFQ_results.EALA_mg_portion__PerWeek(sub_idx); % (EALA)
    nutrition.raw.proteins_aa.Asp(iS) = FFQ_struct.FFQ_results.Asp_mgPerWeek(sub_idx); % (EASP)
    nutrition.raw.proteins_aa.Glu(iS) = FFQ_struct.FFQ_results.Glu_mgPerWeek(sub_idx); % (EGLU)
    nutrition.raw.proteins_aa.Gly(iS) = FFQ_struct.FFQ_results.Gly_mgPerWeek(sub_idx); % (EGLY)
    nutrition.raw.proteins_aa.Pro(iS) = FFQ_struct.FFQ_results.EPRO_mg_portion__PerWeek(sub_idx); % (EPRO)
    nutrition.raw.proteins_aa.Ser(iS) = FFQ_struct.FFQ_results.ESER_mg_portion__PerWeek(sub_idx); % (ESER)
    
    % vitamins
    nutrition.raw.vitamins.A_retinolEquivalent(iS) = FFQ_struct.FFQ_results.VA_vitamineA_quivalentsR_tinol___g_portion__PerWeek(sub_idx); % (VA)
    nutrition.raw.vitamins.A_retinol(iS) = FFQ_struct.FFQ_results.VAR_vitamineAR_tinol___g_portion__PerWeek(sub_idx); % (VAR)
    nutrition.raw.vitamins.A_bCarotene(iS) = FFQ_struct.FFQ_results.VAC_vitamineABeta_carot_nes___g_portion__PerWeek(sub_idx); % (VAC)
    nutrition.raw.vitamins.D_calciferole(iS) = FFQ_struct.FFQ_results.VD_vitamineD___g_portion__PerWeek(sub_idx); % (VD)
    nutrition.raw.vitamins.E_tocopherolEquivalent(iS) = FFQ_struct.FFQ_results.VE_vitamineE_quivalentsTocoph_rol___g_portion__PerWeek(sub_idx); % (VE)
    nutrition.raw.vitamins.E_alphaTocopherol(iS) = FFQ_struct.FFQ_results.VEAT_VitamineE_AlphaTocoph_rol___g_portion__PerWeek(sub_idx); % (VEAT)
    nutrition.raw.vitamins.K_phylloquinone(iS) = FFQ_struct.FFQ_results.VK_VitamineK_Phylloquinone___g_portion__PerWeek(sub_idx); % (VK)
    nutrition.raw.vitamins.B1_thiamin(iS) = FFQ_struct.FFQ_results.VB1_VitamineB1_Thiamine___g_portion__PerWeek(sub_idx); % (VB1)
    nutrition.raw.vitamins.B2_riboflavin(iS) = FFQ_struct.FFQ_results.VB2_VitamineB2_Riboflavine___g_portion__PerWeek(sub_idx); % (VB2)
    nutrition.raw.vitamins.B3_niacin(iS) = FFQ_struct.FFQ_results.NA_ugPerWeek(sub_idx); % (VB3)
    nutrition.raw.vitamins.B3_niacinEquivalents(iS) = FFQ_struct.FFQ_results.NE_ugPerWeek(sub_idx); % (VB3A)
    nutrition.raw.vitamins.B5_pantothenicA(iS) = FFQ_struct.FFQ_results.VB5__g_portion__PerWeek(sub_idx); % (VB5)
    nutrition.raw.vitamins.B6_pyridoxin(iS) = FFQ_struct.FFQ_results.VB6__g_portion__PerWeek(sub_idx); % (VB6)
    nutrition.raw.vitamins.B7_biotin(iS) = FFQ_struct.FFQ_results.VB7__g_portion__PerWeek(sub_idx); % (VB7)
    nutrition.raw.vitamins.B9_folicA(iS) = FFQ_struct.FFQ_results.FolicAcid_mgPerWeek(sub_idx); % (VB9G)
    nutrition.raw.vitamins.B12_cobalamin(iS) = FFQ_struct.FFQ_results.VB12__g_portion__PerWeek(sub_idx); % (VB12)
    nutrition.raw.vitamins.C_ascorbicA(iS) = FFQ_struct.FFQ_results.VC_vitamineC___g_portion__PerWeek(sub_idx); % (VC)
    
    % minerals
    nutrition.raw.minerals.minerals_ashes(iS) = FFQ_struct.FFQ_results.ZM_min_raux__mg_portion__PerWeek(sub_idx); % (ZM) minerals and ashes
    nutrition.raw.minerals.Na(iS) = FFQ_struct.FFQ_results.MNA_sodium__mg_portion__PerWeek(sub_idx); % (MNA) sodium
    nutrition.raw.minerals.K(iS) = FFQ_struct.FFQ_results.MK_potassium__mg_portion__PerWeek(sub_idx); % (MK) potassium
    nutrition.raw.minerals.Ca(iS) = FFQ_struct.FFQ_results.MCA_calcium__mg_portion__PerWeek(sub_idx); % (MCA) calcium
    nutrition.raw.minerals.Mg(iS) = FFQ_struct.FFQ_results.Mg_mgPerWeek(sub_idx); % (MMG) magnesium
    nutrition.raw.minerals.P(iS) = FFQ_struct.FFQ_results.MP_phosphore__mg_portion__PerWeek(sub_idx); % (MP) phosphorus
    nutrition.raw.minerals.S(iS) = FFQ_struct.FFQ_results.MS_mg_portion__PerWeek(sub_idx); % (MS) sulfur
    nutrition.raw.minerals.Cl(iS) = FFQ_struct.FFQ_results.MCL_mg_portion__PerWeek(sub_idx); % (MCL) chlorure
    nutrition.raw.minerals.Fe(iS) = FFQ_struct.FFQ_results.MFE_fer___g_portion__PerWeek(sub_idx); % (MFE)
    nutrition.raw.minerals.Zn(iS) = FFQ_struct.FFQ_results.MZN__g_portion__PerWeek(sub_idx); % (MZN) zinc
    nutrition.raw.minerals.Cu(iS) = FFQ_struct.FFQ_results.MCU_cuivre___g_portion__PerWeek(sub_idx); % (MCU)
    nutrition.raw.minerals.Mn(iS) = FFQ_struct.FFQ_results.MMN__g_portion__PerWeek(sub_idx); % (MMN) manganese
    nutrition.raw.minerals.F(iS) = FFQ_struct.FFQ_results.MF__g_portion__PerWeek(sub_idx); % (MF) fluorure
    nutrition.raw.minerals.I(iS) = FFQ_struct.FFQ_results.MJ__g_portion__PerWeek(sub_idx); % (MJ) iodure
    nutrition.raw.minerals.NaCl_salt(iS) = FFQ_struct.FFQ_results.GMKO_sel__mg_portion__PerWeek(sub_idx); % (GMKO)
    
    % glucids
    nutrition.raw.carbohydrates.absorbableCarbohydrates(iS) = FFQ_struct.FFQ_results.glucids_mgPerWeek(sub_idx); % (ZK)
    nutrition.raw.carbohydrates.mannitol(iS) = FFQ_struct.FFQ_results.KAM_mg_portion__PerWeek(sub_idx); % (KAM)
    nutrition.raw.carbohydrates.sorbitol(iS) = FFQ_struct.FFQ_results.KAS_mg_portion__PerWeek(sub_idx); % (KAS)
    nutrition.raw.carbohydrates.xylitol(iS) = FFQ_struct.FFQ_results.KAX_mg_portion__PerWeek(sub_idx); % (KAX)
    nutrition.raw.carbohydrates.sugarAlcohols(iS) = FFQ_struct.FFQ_results.KA_mg_portion__PerWeek(sub_idx); % (KA)
    nutrition.raw.carbohydrates.Glc(iS) = FFQ_struct.FFQ_results.Glc_mgPerWeek(sub_idx); % (KMT)
    nutrition.raw.carbohydrates.fructose(iS) = FFQ_struct.FFQ_results.KMF_fructose__mg_portion__PerWeek(sub_idx); % (KMF)
    nutrition.raw.carbohydrates.galactose(iS) = FFQ_struct.FFQ_results.KMG_galactose__mg_portion__PerWeek(sub_idx); % (KMG)
    nutrition.raw.carbohydrates.monosaccharides(iS) = FFQ_struct.FFQ_results.KM_monosaccharides__mg_portion__PerWeek(sub_idx); % (KM)
    nutrition.raw.carbohydrates.saccharose(iS) = FFQ_struct.FFQ_results.KDS_saccharose__mg_portion__PerWeek(sub_idx); % (KDS)
    nutrition.raw.carbohydrates.maltose(iS) = FFQ_struct.FFQ_results.KDM_maltose__mg_portion__PerWeek(sub_idx); % (KDM)
    nutrition.raw.carbohydrates.lactose(iS) = FFQ_struct.FFQ_results.Lactose_mgPerWeek(sub_idx); % (KDL)
    nutrition.raw.carbohydrates.disaccharides(iS) = FFQ_struct.FFQ_results.KD_mg_portion__PerWeek(sub_idx); % (KD)
    nutrition.raw.carbohydrates.sugar(iS) = FFQ_struct.FFQ_results.sugars_mgPerWeek(sub_idx); % (KMD)
    nutrition.raw.carbohydrates.reabsorbableOligosaccharides(iS) = FFQ_struct.FFQ_results.KPOR_mg_portion__PerWeek(sub_idx); % (KPOR)
    nutrition.raw.carbohydrates.nonReabsorbableOligosaccharides(iS) = FFQ_struct.FFQ_results.KPON_mg_portion__PerWeek(sub_idx); % (KPON)
    nutrition.raw.carbohydrates.glycogen_animalStarch(iS) = FFQ_struct.FFQ_results.KPG_glycog_ne_amidonAnimal___mg_portion__PerWeek(sub_idx); % (KPG)
    nutrition.raw.carbohydrates.starch(iS) = FFQ_struct.FFQ_results.KPS_amidon__mg_portion__PerWeek(sub_idx); % (KPS)
    nutrition.raw.carbohydrates.polysaccharides(iS) = FFQ_struct.FFQ_results.KP_mg_portion__PerWeek(sub_idx); % (KP)
    nutrition.raw.carbohydrates.polypentoses(iS) = FFQ_struct.FFQ_results.KBP_mg_portion__PerWeek(sub_idx); % (KBP)
    nutrition.raw.carbohydrates.polyhexoses(iS) = FFQ_struct.FFQ_results.KBH_mg_portion__PerWeek(sub_idx); % (KBH)
    nutrition.raw.carbohydrates.polyuronicA(iS) = FFQ_struct.FFQ_results.KBU_mg_portion__PerWeek(sub_idx); % (KBU)
    
    % fibers (kind of carbohydrates in fact, but with very different
    % nutritional properties)
    nutrition.raw.fibers.fibers(iS) = FFQ_struct.FFQ_results.ZB_fibres__mg_portion__PerWeek(sub_idx); % (ZB)
    nutrition.raw.carbohydrates.solubleFibers(iS) = FFQ_struct.FFQ_results.KBW_fibresSolublesDansL_eau__mg_portion__PerWeek(sub_idx); % (KBW)
    nutrition.raw.carbohydrates.insolubleFibers(iS) = FFQ_struct.FFQ_results.KBN_fibresInsolublesDansL_eau__mg_portion__PerWeek(sub_idx); % (KBN)
    nutrition.raw.fibers.cellulose(iS) = FFQ_struct.FFQ_results.KBC_cellulose__mg_portion__PerWeek(sub_idx); % (KBC)
    nutrition.raw.fibers.lignine(iS) = FFQ_struct.FFQ_results.KBL_mg_portion__PerWeek(sub_idx); % (KBL)
    
    % fatty acids
    nutrition.raw.fa.fat(iS) = FFQ_struct.FFQ_results.fat_mgPerWeek(sub_idx); % (ZF)
    nutrition.raw.fa.fa_monoInsaturated(iS) = FFQ_struct.FFQ_results.FU_acidesGrasMonoinsatur_s__mg_portion__PerWeek(sub_idx); % (FU)
    nutrition.raw.fa.fa_polyInsaturated(iS) = FFQ_struct.FFQ_results.FP_acidesGrasPolyinsatur_s__mg_portion__PerWeek(sub_idx); % (FP)
    nutrition.raw.fa.fa_shortChain(iS) = FFQ_struct.FFQ_results.FK_acidesGras_Cha_neCourte__mg_portion__PerWeek(sub_idx); % (FK)
    nutrition.raw.fa.fa_mediumChain(iS) = FFQ_struct.FFQ_results.FM_acidesGras_Cha_neMoyenne__mg_portion__PerWeek(sub_idx); % (FM)
    nutrition.raw.fa.fa_longChain(iS) = FFQ_struct.FFQ_results.FL_acidesGras_Cha_neLongue__mg_portion__PerWeek(sub_idx); % (FL)
    nutrition.raw.fa.butanoicA(iS) = FFQ_struct.FFQ_results.F40_mg_portion__PerWeek(sub_idx); % (F40)
    nutrition.raw.fa.hexanoicA(iS) = FFQ_struct.FFQ_results.F60_mg_portion__PerWeek(sub_idx); % (F60)
    nutrition.raw.fa.octanoicA(iS) = FFQ_struct.FFQ_results.F80_mg_portion__PerWeek(sub_idx); % (F80)
    nutrition.raw.fa.decanoicA_capricA(iS) = FFQ_struct.FFQ_results.F100_mg_portion__PerWeek(sub_idx); % (F100)
    nutrition.raw.fa.dodecanoicA(iS) = FFQ_struct.FFQ_results.F120_mg_portion__PerWeek(sub_idx); % (F120)
    nutrition.raw.fa.tetradecanoicA(iS) = FFQ_struct.FFQ_results.F140_mg_portion__PerWeek(sub_idx); % (F140)
    nutrition.raw.fa.pentadecanoicA(iS) = FFQ_struct.FFQ_results.F150_mg_portion__PerWeek(sub_idx); % (F150)
    nutrition.raw.fa.hexadecanoicA(iS) = FFQ_struct.FFQ_results.F160_mg_portion__PerWeek(sub_idx); % (F160)
    nutrition.raw.fa.heptadecanoicA(iS) = FFQ_struct.FFQ_results.F170_mg_portion__PerWeek(sub_idx); % (F170)
    nutrition.raw.fa.octadecanoicA(iS) = FFQ_struct.FFQ_results.F180_mg_portion__PerWeek(sub_idx); % (F180)
    nutrition.raw.fa.eicosanoicA(iS) = FFQ_struct.FFQ_results.F200_mg_portion__PerWeek(sub_idx); % (F200)
    nutrition.raw.fa.decanoicA_behenicA(iS) = FFQ_struct.FFQ_results.F220_mg_portion__PerWeek(sub_idx); % (F220)
    nutrition.raw.fa.tetracosanoicA(iS) = FFQ_struct.FFQ_results.F240_mg_portion__PerWeek(sub_idx); % (F240)
    nutrition.raw.fa.saturatedFA(iS) = FFQ_struct.FFQ_results.FS_acidesGrasSatur_s__mg_portion__PerWeek(sub_idx); % (FS)
    nutrition.raw.fa.tetradecenoicA(iS) = FFQ_struct.FFQ_results.F141_mg_portion__PerWeek(sub_idx); % (F141)
    nutrition.raw.fa.pentadecenoicA(iS) = FFQ_struct.FFQ_results.F151_mg_portion__PerWeek(sub_idx); % (F151)
    nutrition.raw.fa.hexadecenoicA(iS) = FFQ_struct.FFQ_results.F161_mg_portion__PerWeek(sub_idx); % (F161)
    nutrition.raw.fa.heptadecenoicA(iS) = FFQ_struct.FFQ_results.F171_mg_portion__PerWeek(sub_idx); % (F171)
    nutrition.raw.fa.octadecenoicA(iS) = FFQ_struct.FFQ_results.F181_mg_portion__PerWeek(sub_idx); % (F181)
    nutrition.raw.fa.eicosenoicA(iS) = FFQ_struct.FFQ_results.F201_mg_portion__PerWeek(sub_idx); % (F201)
    nutrition.raw.fa.decosenoicA(iS) = FFQ_struct.FFQ_results.F221_mg_portion__PerWeek(sub_idx); % (F221)
    nutrition.raw.fa.tetracosenoicA(iS) = FFQ_struct.FFQ_results.F241_mg_portion__PerWeek(sub_idx); % (F241)
    nutrition.raw.fa.hexadecadienoicA(iS) = FFQ_struct.FFQ_results.F162_mg_portion__PerWeek(sub_idx); % (F162)
    nutrition.raw.fa.hexacedatetraenoicA(iS) = FFQ_struct.FFQ_results.F164_mg_portion__PerWeek(sub_idx); % (F164)
    nutrition.raw.fa.octadecadienoicA_linoleicA(iS) = FFQ_struct.FFQ_results.F182_mg_portion__PerWeek(sub_idx); % (F182)
    nutrition.raw.fa.octadecatrienoicA_linolenicA(iS) = FFQ_struct.FFQ_results.F183_mg_portion__PerWeek(sub_idx); % (F183)
    nutrition.raw.fa.octradecatetraenoicA(iS) = FFQ_struct.FFQ_results.F184_mg_portion__PerWeek(sub_idx); % (F184)
    nutrition.raw.fa.nonAdecatrienoicA(iS) = FFQ_struct.FFQ_results.F193_mg_portion__PerWeek(sub_idx); % (F193)
    nutrition.raw.fa.eicosadienoicA(iS) = FFQ_struct.FFQ_results.F202_mg_portion__PerWeek(sub_idx); % (F202)
    nutrition.raw.fa.eicosatrienoicA(iS) = FFQ_struct.FFQ_results.F203_mg_portion__PerWeek(sub_idx); % (F203)
    nutrition.raw.fa.eicosatetraenoicA_arachidonicA(iS) = FFQ_struct.FFQ_results.F204_mg_portion__PerWeek(sub_idx); % (F204)
    nutrition.raw.fa.eicosapentaenoicA(iS) = FFQ_struct.FFQ_results.F205_mg_portion__PerWeek(sub_idx); % (F205)
    nutrition.raw.fa.docosadienoicA(iS) = FFQ_struct.FFQ_results.F222_mg_portion__PerWeek(sub_idx); % (F222)
    nutrition.raw.fa.docosatrienoicA(iS) = FFQ_struct.FFQ_results.F223_mg_portion__PerWeek(sub_idx); % (F223)
    nutrition.raw.fa.docosatetraenoicA(iS) = FFQ_struct.FFQ_results.F224_mg_portion__PerWeek(sub_idx); % (F224)
    nutrition.raw.fa.docosapentaenoicA(iS) = FFQ_struct.FFQ_results.F225_mg_portion__PerWeek(sub_idx); % (F225)
    nutrition.raw.fa.docosahexaenoicA(iS) = FFQ_struct.FFQ_results.F226_mg_portion__PerWeek(sub_idx); % (F226)
    nutrition.raw.fa.omegas3(iS) = FFQ_struct.FFQ_results.omega3_mgPerWeek(sub_idx); % (FO3)
    nutrition.raw.fa.omegas6(iS) = FFQ_struct.FFQ_results.FO6_omega6__mg_portion__PerWeek(sub_idx); % (FO6)
    nutrition.raw.fa.glycerine_lipids(iS) = FFQ_struct.FFQ_results.FG_mg_portion__PerWeek(sub_idx); % (FG)
    nutrition.raw.fa.cholesterol(iS) = FFQ_struct.FFQ_results.FC_Chol_st_rol__mg_portion__PerWeek(sub_idx); % (FC)
    nutrition.raw.fa.GFPS_PolyI_div_Sat_fa_ratio(iS) = FFQ_struct.FFQ_results.GFPSPolyins_saturatedFattyAcidsRatio_mg_portion__PerWeek(sub_idx); % (GFPS)
    
    % general/unclassifiable
    nutrition.raw.other.water(iS) = FFQ_struct.FFQ_results.ZW_eau__mg_portion__PerWeek(sub_idx); % (ZW)
    nutrition.raw.other.organicAcids(iS) = FFQ_struct.FFQ_results.ZO_acidesOrganiques__mg_portion__PerWeek(sub_idx); % (ZO)
    nutrition.raw.other.alcohol(iS) = FFQ_struct.FFQ_results.ZA_alcool__mg_portion__PerWeek(sub_idx); % (ZA)
    nutrition.raw.other.uricA(iS) = FFQ_struct.FFQ_results.EH_mg_portion__PerWeek(sub_idx); % (EH)
    nutrition.raw.other.purin(iS) = FFQ_struct.FFQ_results.EP_mg_portion__PerWeek(sub_idx); % (EP)
    nutrition.raw.bread(iS) = FFQ_struct.FFQ_results.GKB_pain__BE__PerWeek(sub_idx); % (GKB)
    
end % subject loop

%% normalize data by total caloric intake
for iCat = 1:n_categ
    categ_nm = nutri_categories_names{iCat};
    for iNutrient = 1:n_nutri.(categ_nm)
        nutri_nm = nutri_names.(categ_nm){iNutrient};
        nutrition.norm_by_totalCal.(categ_nm).(nutri_nm) = nutrition.raw.(categ_nm).(nutri_nm)./nutrition.raw.energy.kCal;
    end
end % category

end % function