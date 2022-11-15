% script to correlate niacin nutritional intake derived from the FFQ to
% blood measurements of the NAD metabolome.
%
%

%% working directory
pcRoot = LGCM_root_paths;
switch pcRoot
    case 'E:\'
        gitPath = fullfile('C:','Users','clairis','Desktop','GitHub',...
            'LGC_motiv');
    case 'L:\human_data_private\raw_data_subject\'
        gitPath = fullfile('C:','Users','Loco','Documents','GitHub',...
            'LGC_motiv');
end

%% select study and participants
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);
dataResultsPath = fullfile(gitPath, 'LGC_Motiv_results',study_nm);

%% extract all blood data
bloodFilePath = fullfile(dataResultsPath,...
    'blood_NAD','20221114_CS_NAD_whole_blood_data.xlsx');
NAD_blood_table = readtable(bloodFilePath,...
    'Sheet','clean');
%% nutrition path
nutritionPath = [dataResultsPath, filesep,'nutrition', filesep];
%% extract all niacin data
niacineFilePath = [nutritionPath,'niacine_scoring.xlsx'];
niacin_table = readtable(niacineFilePath,...
    'Sheet','Sheet1');
%% extract all Tryptophane data
TrpFilePath = [nutritionPath,'Tryptophan_scoring.xlsx'];
Trp_table = readtable(TrpFilePath,...
    'Sheet','Sheet1');
%% extract calories
calFilePath = [nutritionPath,'calories_scoring.xlsx'];
calories_table = readtable(calFilePath,...
    'Sheet','Sheet1');
%% extract general infos of subjects
galFilePath = [dataResultsPath, filesep,...
    'summary_participants_infos.xlsx'];
sub_infos = readtable(galFilePath,...
    'Sheet','Feuil1');
%% loop through subjects to perform the correlation
blood_metabolites = {'Nam','MeXPY','NAD','NADP','NMN','MeNam','NR',...
    'NADH','NADPH','NAD_div_NADH','NADP_div_NADPH'};
nNADHmetab = length(blood_metabolites);
[niacin,...
    Trp,...
    weight,...
    BMI,...
    calories] = deal(NaN(1,NS));
for iBloodM = 1:nNADHmetab
    [blood.(blood_metabolites{iBloodM}),...
        blood.([blood_metabolites{iBloodM},'_div_weight']),...
        blood.([blood_metabolites{iBloodM},'_div_BMI']),...
        blood.([blood_metabolites{iBloodM},'_div_calories'])]= deal(NaN(1,NS));
end
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    % identify index for the current subject on the different files
    sub_blood_idx = find(strcmp(NAD_blood_table.ID, sub_nm_bis));
    sub_niacin_idx = find(strcmp(niacin_table.CID, sub_nm));
    sub_Trp_idx = find(strcmp(Trp_table.CID, sub_nm));
    sub_calories_idx = find(strcmp(calories_table.CID, sub_nm));
    sub_infos_idx = find(strcmp(sub_infos.CID, sub_nm_bis));
    % extract niacin values
    if ~isempty(sub_niacin_idx) && size(sub_niacin_idx,1) == 1
        niacin(iS) = niacin_table.NiacineParSemaine_ug_(sub_niacin_idx);
    end
    % extract Tryptophan values
    if ~isempty(sub_Trp_idx) && size(sub_Trp_idx,1) == 1
        Trp(iS) = Trp_table.TryptophaneParSemaine_ug_(sub_Trp_idx);
    end
    % extract calories values
    if ~isempty(sub_calories_idx) && size(sub_calories_idx,1) == 1
        calories(iS) = calories_table.CaloriesParSemaine_kcal_(sub_calories_idx);
    end
    % extract blood values
    for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        if ~isempty(sub_blood_idx)
            blood.(blood_m_nm)(iS) = NAD_blood_table.(blood_m_nm)(sub_blood_idx);
        end
    end % metabolite loop
    % extract subject infos
    if ~isempty(sub_infos_idx) && size(sub_infos_idx,1) == 1
        weight(iS) = sub_infos.Weight(sub_infos_idx);
        BMI(iS) = sub_infos.BMI(sub_infos_idx);
    end
end % subject loop

% normalize niacin values by caloric intake-related variables
niacin_div_weight = niacin./weight;
niacin_div_BMI = niacin./BMI;
niacin_div_calories = niacin./calories;
Trp_div_weight = Trp./weight;
Trp_div_BMI = Trp./BMI;
Trp_div_calories = Trp./calories;
niacinPlusTrp = niacin + Trp;
niacinPlusTrp_div_weight = niacinPlusTrp./weight;
niacinPlusTrp_div_BMI = niacinPlusTrp./BMI;
niacinPlusTrp_div_calories = niacinPlusTrp./calories;
normTypes = {'niacin',...
    'niacin_div_weight','niacin_div_BMI',...
    'niacin_div_calories',...
    'Trp',...
    'Trp_div_weight','Trp_div_BMI',...
    'Trp_div_calories',...
    'niacin_plus_Trp',...
    'niacin_plus_Trp_div_weight','niacin_plus_Trp_div_BMI',...
    'niacin_plus_Trp_div_calories'};
nNormTypes = length(normTypes);

%% test the correlations
for iNormType = 1:nNormTypes
    normType = normTypes{iNormType};
    switch normType
        case 'niacin'
            nutrition_var = niacin;
        case 'niacin_div_weight'
            nutrition_var = niacin_div_weight;
        case 'niacin_div_BMI'
            nutrition_var = niacin_div_BMI;
        case 'niacin_div_calories'
            nutrition_var = niacin_div_calories;
        case 'Trp'
            nutrition_var = Trp;
        case 'Trp_div_weight'
            nutrition_var = Trp_div_weight;
        case 'Trp_div_BMI'
            nutrition_var = Trp_div_BMI;
        case 'Trp_div_calories'
            nutrition_var = Trp_div_calories;
        case 'niacin_plus_Trp'
            nutrition_var = niacinPlusTrp;
        case 'niacin_plus_Trp_div_weight'
            nutrition_var = niacinPlusTrp_div_weight;
        case 'niacin_plus_Trp_div_BMI'
            nutrition_var = niacinPlusTrp_div_BMI;
        case 'niacin_plus_Trp_div_calories'
            nutrition_var = niacinPlusTrp_div_calories;
    end
    for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        % check good subjects for each measurement
        goodSubs = (~isnan(blood.(blood_m_nm)).*(~isnan(nutrition_var))) == 1;
        blood_f_niacin_nm = [blood_m_nm,'_f_',normType];
        [b_tmp,~,stats_tmp] = glmfit(...
            nutrition_var(goodSubs), blood.(blood_m_nm)(goodSubs),'normal');
        corr_vals.b.(blood_f_niacin_nm) = b_tmp(2);
        corr_vals.p.(blood_f_niacin_nm) = stats_tmp.p(2);
        fit_blood.(blood_m_nm).(normType) = NaN(1,NS);
        fit_blood.(blood_m_nm).(normType)(goodSubs) = glmval(b_tmp, nutrition_var(goodSubs), 'identity');
    end % blood metabolite
end % normalization types loop
%% display the correlation matrix
lWidth = 3;
pSize = 30;

% display each correlation
for iNormType = 1:nNormTypes
    normType = normTypes{iNormType};
    switch normType
        case 'niacin'
            nutrition_var = niacin;
        case 'niacin_div_weight'
            nutrition_var = niacin_div_weight;
        case 'niacin_div_BMI'
            nutrition_var = niacin_div_BMI;
        case 'niacin_div_calories'
            nutrition_var = niacin_div_calories;
        case 'Trp'
            nutrition_var = Trp;
        case 'Trp_div_weight'
            nutrition_var = Trp_div_weight;
        case 'Trp_div_BMI'
            nutrition_var = Trp_div_BMI;
        case 'Trp_div_calories'
            nutrition_var = Trp_div_calories;
        case 'niacin_plus_Trp'
            nutrition_var = niacinPlusTrp;
        case 'niacin_plus_Trp_div_weight'
            nutrition_var = niacinPlusTrp_div_weight;
        case 'niacin_plus_Trp_div_BMI'
            nutrition_var = niacinPlusTrp_div_BMI;
        case 'niacin_plus_Trp_div_calories'
            nutrition_var = niacinPlusTrp_div_calories;
    end
    fig;
    for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        % check good subjects for each measurement
        goodSubs = (~isnan(blood.(blood_m_nm)).*(~isnan(nutrition_var))) == 1;
        
        subplot(3,4,iBloodM);
        hold on;
        hdl = scatter(nutrition_var, blood.(blood_m_nm));
        [nutrition_sorted, nutrition_sort_idx] = sort(nutrition_var(goodSubs));
        fit_blood_tmp = fit_blood.(blood_m_nm).(normType)(goodSubs);
        fit_blood_sorted = fit_blood_tmp(nutrition_sort_idx);
        fit_hdl = plot(nutrition_sorted, fit_blood_sorted);
        switch normType
            case 'niacin'
                xlabel('niacin (μg/week)');
            case 'niacin_div_weight'
                xlabel('niacin/weight');
            case 'niacin_div_BMI'
                xlabel('niacin/BMI');
            case 'niacin_div_calories'
                xlabel('niacin/calories');
            case 'Trp'
                xlabel('Trp (μg/week)');
            case 'Trp_div_weight'
                xlabel('Trp/weight');
            case 'Trp_div_BMI'
                xlabel('Trp/BMI');
            case 'Trp_div_calories'
                xlabel('Trp/calories');
            case 'niacin_plus_Trp'
                xlabel('niacin+Trp (μg/week)');
            case 'niacin_plus_Trp_div_weight'
                xlabel('(niacin+Trp)/weight');
            case 'niacin_plus_Trp_div_BMI'
                xlabel('(niacin+Trp)/BMI');
            case 'niacin_plus_Trp_div_calories'
                xlabel('(niacin+Trp)/calories');
        end
        
        switch blood_m_nm
            case 'NAD_div_NADH'
                ylabel('NAD/NADH (μM)');
            case 'NADP_div_NADPH'
                ylabel('NADP/NADPH (μM)');
            otherwise
                ylabel([blood_m_nm,' (μM)']);
        end
        % esthetic improvements
        hdl.LineWidth = lWidth;
        hdl.MarkerEdgeColor = [0 0 0];
        fit_hdl.LineWidth = lWidth;
        fit_hdl.Color = [143 143 143]./255;
        fit_hdl.LineStyle = '--';
        legend_size(pSize);
    end % blood metabolite
end % normalization types loop

%% check which values are significant
correl_names = fieldnames(corr_vals.p);
for iCorr = 1:length(correl_names)
    correl_nm = correl_names{iCorr};
    if corr_vals.p.(correl_nm) < 0.05
        signif_correl.p005.(correl_nm) = corr_vals.p.(correl_nm);
    end
end % loop through tests