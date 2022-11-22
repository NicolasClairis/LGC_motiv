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
[bloodTable, blood_NAD_sub_List] = load_blood_NAD(study_nm);
%% nutrition path
nutritionPath = [dataResultsPath, filesep,'nutrition', filesep];
%% extract all niacin data
niacineFilePath = [nutritionPath,'niacine_scoring.xlsx'];
niacin_table = readtable(niacineFilePath,...
    'Sheet','Sheet1');
%% extract niacin equivalents data
niacineEquivFilePath = [nutritionPath,'niacine_equivalents_scoring.xlsx'];
niacin_equivalents_table = readtable(niacineEquivFilePath,...
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
blood_metabolites = fieldnames(bloodTable);
nNADHmetab = length(blood_metabolites);
[niacin,...
    niacin_equivalents,...
    Trp,...
    weight,...
    BMI,...
    calories] = deal(NaN(1,NS));
for iBloodM = 1:nNADHmetab
    [blood.(blood_metabolites{iBloodM}),...
        blood.([blood_metabolites{iBloodM},'_div_weight']),...
        blood.([blood_metabolites{iBloodM},'_div_BMI']),...
        blood.([blood_metabolites{iBloodM},'_div_calories_per_weight']),...
        blood.([blood_metabolites{iBloodM},'_div_calories'])]= deal(NaN(1,NS));
end
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    % identify index for the current subject on the different files
    sub_blood_idx = find(strcmp(blood_NAD_sub_List, sub_nm_bis));
    sub_niacin_idx = find(strcmp(niacin_table.CID, sub_nm));
    sub_niacinEquiv_idx = find(strcmp(niacin_equivalents_table.CID, sub_nm));
    sub_Trp_idx = find(strcmp(Trp_table.CID, sub_nm));
    sub_calories_idx = find(strcmp(calories_table.CID, sub_nm));
    sub_infos_idx = find(strcmp(sub_infos.CID, sub_nm_bis));
    % extract niacin values
    if ~isempty(sub_niacin_idx) && size(sub_niacin_idx,1) == 1
        niacin(iS) = niacin_table.NiacineParSemaine__g_(sub_niacin_idx)./1000;
    end
    % extract niacin equivalents values
    if ~isempty(sub_niacinEquiv_idx) && size(sub_niacinEquiv_idx,1) == 1
        niacin_equivalents(iS) = niacin_equivalents_table.x_quivalentsDeNiacineParSemaine__g_(sub_niacinEquiv_idx)./1000;
    end
    % extract Tryptophan values
    if ~isempty(sub_Trp_idx) && size(sub_Trp_idx,1) == 1
        Trp(iS) = Trp_table.TryptophaneParSemaine_mg_(sub_Trp_idx);
    end
    % extract calories values
    if ~isempty(sub_calories_idx) && size(sub_calories_idx,1) == 1
        calories(iS) = calories_table.CaloriesParSemaine_kcal_(sub_calories_idx);
    end
    % extract blood values
    for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        if ~isempty(sub_blood_idx)
            blood.(blood_m_nm)(iS) = bloodTable.(blood_m_nm)(sub_blood_idx);
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
niacin_div_calories_per_BMI = (niacin./calories)./BMI;
Trp_div_weight = Trp./weight;
Trp_div_BMI = Trp./BMI;
Trp_div_calories = Trp./calories;
Trp_div_calories_per_BMI = (Trp./calories)./BMI;
niacinEquiv_div_weight = niacin_equivalents./weight;
niacinEquiv_div_BMI = niacin_equivalents./BMI;
niacinEquiv_div_calories = niacin_equivalents./calories;
niacinEquiv_div_calories_per_BMI = (niacin_equivalents./calories)./BMI;
normTypes = {'niacin',...
    'niacin_div_weight','niacin_div_BMI',...
    'niacin_div_calories','niacin_div_calories_per_BMI',...
    'Trp',...
    'Trp_div_weight','Trp_div_BMI',...
    'Trp_div_calories','Trp_div_calories_per_BMI',...
    'niacinEquiv',...
    'niacinEquiv_div_weight','niacinEquiv_div_BMI',...
    'niacinEquiv_div_calories','niacinEquiv_div_calories_per_BMI'};
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
        case 'niacin_div_calories_per_BMI'
            nutrition_var = niacin_div_calories_per_BMI;
        case 'Trp'
            nutrition_var = Trp;
        case 'Trp_div_weight'
            nutrition_var = Trp_div_weight;
        case 'Trp_div_BMI'
            nutrition_var = Trp_div_BMI;
        case 'Trp_div_calories'
            nutrition_var = Trp_div_calories;
        case 'Trp_div_calories_per_BMI'
            nutrition_var = Trp_div_calories_per_BMI;
        case 'niacinEquiv'
            nutrition_var = niacin_equivalents;
        case 'niacinEquiv_div_weight'
            nutrition_var = niacinEquiv_div_weight;
        case 'niacinEquiv_div_BMI'
            nutrition_var = niacinEquiv_div_BMI;
        case 'niacinEquiv_div_calories'
            nutrition_var = niacinEquiv_div_calories;
        case 'niacinEquiv_div_calories_per_BMI'
            nutrition_var = niacinEquiv_div_calories_per_BMI;
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
        case 'niacin_div_calories_per_BMI'
            nutrition_var = niacin_div_calories_per_BMI;
        case 'Trp'
            nutrition_var = Trp;
        case 'Trp_div_weight'
            nutrition_var = Trp_div_weight;
        case 'Trp_div_BMI'
            nutrition_var = Trp_div_BMI;
        case 'Trp_div_calories'
            nutrition_var = Trp_div_calories;
        case 'Trp_div_calories_per_BMI'
            nutrition_var = Trp_div_calories_per_BMI;
        case 'niacinEquiv'
            nutrition_var = niacin_equivalents;
        case 'niacinEquiv_div_weight'
            nutrition_var = niacinEquiv_div_weight;
        case 'niacinEquiv_div_BMI'
            nutrition_var = niacinEquiv_div_BMI;
        case 'niacinEquiv_div_calories'
            nutrition_var = niacinEquiv_div_calories;
        case 'niacinEquiv_div_calories_per_BMI'
            nutrition_var = niacinEquiv_div_calories_per_BMI;
    end
    
    % one figure for all measures and one for integrated measures
    fig1 = fig; j_fig1 = 0;
    fig2 = fig; j_fig2 = 0;
    
    for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        % check good subjects for each measurement
        goodSubs = (~isnan(blood.(blood_m_nm)).*(~isnan(nutrition_var))) == 1;
        
        switch blood_m_nm
            case {'Nam','NMN','NR','NAD','NADH','NADP','NADPH','MeNam','MeXPY'}
                figure(fig1);
                j_fig1 = j_fig1 + 1;
                subplot(3,4,j_fig1);
            case {'NAD_div_NADH','NADP_div_NADPH',...
                    'total_NAD_precursors','total_NAD',...
                    'total_NAD_with_precursors',...
                    'total_NAD_with_byproducts','total_NAD_byproducts'}
                figure(fig2);
                j_fig2 = j_fig2 + 1;
                subplot(3,3,j_fig2);
        end
        hold on;
        hdl = scatter(nutrition_var, blood.(blood_m_nm));
        [nutrition_sorted, nutrition_sort_idx] = sort(nutrition_var(goodSubs));
        fit_blood_tmp = fit_blood.(blood_m_nm).(normType)(goodSubs);
        fit_blood_sorted = fit_blood_tmp(nutrition_sort_idx);
        fit_hdl = plot(nutrition_sorted, fit_blood_sorted);
        switch normType
            case 'niacin'
                xlabel('niacin (mg/week)');
            case 'niacin_div_weight'
                xlabel('niacin/weight');
            case 'niacin_div_BMI'
                xlabel('niacin/BMI');
            case 'niacin_div_calories'
                xlabel('niacin/calories');
            case 'niacin_div_calories_per_BMI'
                xlabel('niacin/calories/BMI');
            case 'Trp'
                xlabel('Trp (mg/week)');
            case 'Trp_div_weight'
                xlabel('Trp/weight');
            case 'Trp_div_BMI'
                xlabel('Trp/BMI');
            case 'Trp_div_calories'
                xlabel('Trp/calories');
            case 'Trp_div_calories_per_BMI'
                xlabel('Trp/calories/BMI');
            case 'niacinEquiv'
                xlabel('niacin eq. (mg/week)');
            case 'niacinEquiv_div_weight'
                xlabel('niacin eq./weight');
            case 'niacinEquiv_div_BMI'
                xlabel('niacin eq./BMI');
            case 'niacinEquiv_div_calories'
                xlabel('niacin eq./calories');
            case 'niacinEquiv_div_calories_per_BMI'
                xlabel('niacin eq./calories/BMI');
        end
        
        [blood_labelname] = blood_label(blood_m_nm);
        ylabel(blood_labelname);
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
        signif_correl.p005.p.(correl_nm) = corr_vals.p.(correl_nm);
        signif_correl.p005.b.(correl_nm) = corr_vals.b.(correl_nm);
    elseif corr_vals.p.(correl_nm) > 0.05 && corr_vals.p.(correl_nm) < 0.1
        signif_correl.p01.p.(correl_nm) = corr_vals.p.(correl_nm);
        signif_correl.p01.b.(correl_nm) = corr_vals.b.(correl_nm);
    end
end % loop through tests

%% check min/max subjects missing for each blood measure
n_maxNaN = 0;
for iB = 1:nNADHmetab
    blood_nm = blood_metabolites{iB};
    n_maxNaN = max(n_maxNaN, sum(isnan(blood.(blood_nm))));
end
disp(['maximum ',num2str(n_maxNaN),'/',num2str(NS),' subjects with NaN blood values.']);

n_minNaN = NS;
for iB = 1:nNADHmetab
    blood_nm = blood_metabolites{iB};
    n_minNaN = min(n_minNaN, sum(isnan(blood.(blood_nm))));
end
disp(['minimum ',num2str(NS-n_minNaN),'/',num2str(NS),' subjects with good blood values.']);
