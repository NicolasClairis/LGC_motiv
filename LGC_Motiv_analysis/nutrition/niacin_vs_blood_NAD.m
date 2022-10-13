% script to correlate niacin nutritional intake derived from the FFQ to
% blood measurements of the NAD metabolome.
%
%

%% working directory
pcRoot = LGCM_root_paths;
switch pcRoot
    case 'E:\'
        dataStoragePath = 'M:\human_data_private\analyzed_data\study1';
        gitPath = fullfile('C:','Users','clairis','Desktop','GitHub',...
            'LGC_motiv');
    case 'L:\human_data_private\raw_data_subject\'
        dataStoragePath = 'L:\human_data_private\analyzed_data\study1';
        gitPath = fullfile('C:','Users','Loco','Documents','GitHub',...
            'LGC_motiv');
end

%% select study and participants
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% extract all blood data
bloodFilePath = [dataStoragePath,filesep,'blood',filesep,...
    '20220915_preliminary_blood_NAD_report.xlsx'];
NAD_blood_table = readtable(bloodFilePath,...
        'Sheet','Final');
%% extract all niacin data
niacineFilePath = [gitPath, filesep, 'LGC_Motiv_analysis' filesep,...
    'nutrition',filesep,'niacine_scoring.xlsx'];
niacin_table = readtable(niacineFilePath,...
        'Sheet','Sheet1');
%% extract general infos of subjects
galFilePath = [gitPath, filesep, 'LGC_Motiv_results', filesep,...
    'summary_participants_infos.xlsx'];
sub_infos = readtable(galFilePath,...
        'Sheet','Sheet1');
%% loop through subjects to perform the correlation
blood_metabolites = {'Nam','MeXPY','NAD','NADP','NMN','MeNam','NR',...
    'NADH','NADPH','NAD_div_NADH','NADP_div_NADPH'};
nNADHmetab = length(blood_metabolites);
[niacin,...
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
    sub_infos_idx = find(strcmp(sub_infos.CID, sub_nm_bis));
    % extract niacin values
    if ~isempty(sub_niacin_idx) && size(sub_niacin_idx,1) == 1
        niacin(iS) = niacin_table.NiacineParSemaine_ug_(sub_niacin_idx);
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
normTypes = {'niacin',...
    'niacin_div_weight','niacin_div_BMI',...
    'niacin_div_calories'};
nNormTypes = length(normTypes);

%% test the correlations
for iNormType = 1:nNormTypes
    normType = normTypes{iNormType};
    switch normType
        case 'niacin'
            niacin_var = niacin;
        case 'niacin_div_weight'
            niacin_var = niacin_div_weight;
        case 'niacin_div_BMI'
            niacin_var = niacin_div_BMI;
        case 'niacin_div_calories'
            niacin_var = niacin_div_calories;
    end
    for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        % check good subjects for each measurement
        goodSubs = (~isnan(blood.(blood_m_nm)).*(~isnan(niacin_var))) == 1;
        niacin_f_blood_nm = [normType,'_f_',blood_m_nm];
        [b_tmp,~,stats_tmp] = glmfit(...
            niacin_var(goodSubs), blood.(blood_m_nm)(goodSubs),'normal');
        corr_vals.b.(niacin_f_blood_nm) = b_tmp(2);
        corr_vals.p.(niacin_f_blood_nm) = stats_tmp.p(2);
        fit_blood.(blood_m_nm).(normType) = NaN(1,NS);
        fit_blood.(blood_m_nm).(normType)(goodSubs) = glmval(b_tmp, niacin_var(goodSubs), 'identity');
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
            niacin_var = niacin;
        case 'niacin_div_weight'
            niacin_var = niacin_div_weight;
        case 'niacin_div_BMI'
            niacin_var = niacin_div_BMI;
        case 'niacin_div_calories'
            niacin_var = niacin_div_calories;
    end
    fig;
    for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        % check good subjects for each measurement
        goodSubs = (~isnan(blood.(blood_m_nm)).*(~isnan(niacin_var))) == 1;
        
        subplot(3,4,iBloodM);
        hold on;
        hdl = scatter(niacin_var, blood.(blood_m_nm));
        [niacin_sorted, niacin_sort_idx] = sort(niacin_var(goodSubs));
        fit_blood_tmp = fit_blood.(blood_m_nm).(normType)(goodSubs);
        fit_blood_sorted = fit_blood_tmp(niacin_sort_idx);
        fit_hdl = plot(niacin_sorted, fit_blood_sorted);
        switch normType
            case 'niacin'
                xlabel('niacin (μg/week)');
            case 'niacin_div_weight'
                xlabel('niacin/weight');
            case 'niacin_div_BMI'
                xlabel('niacin/BMI');
            case 'niacin_div_calories'
                xlabel('niacin/calories');
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