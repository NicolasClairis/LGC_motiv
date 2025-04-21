function[r_corr, pval, signif, NS_goodS] = brainMetabolites_f_nutrition(outlierF)
% [r_corr, pval, signif, NS_goodS] = brainMetabolites_f_nutrition(outlierF)
%
% INPUTS
% outlierF: filter outlier in the variables included in each correlation
% test (1) or not (0)
%
% OUTPUTS
% r_corr: structure with correlation coefficients
%
% pval: structure with p.value
%
% signif: structure with significant correlations
%
% NS_goodS: number of good subjects

%% outlier filtering
if ~exist('outlierF','var') || isempty(outlierF)
    outlierF_nm = questdlg('Outlier filtering?','Outlier filtering',...
        'No','Yes','Yes');
    switch outlierF_nm
        case 'Yes'
            outlierF = 1;
        case 'No'
            outlierF = 0;
    end
end % outlier filtering

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% working directory
root = LGCM_root_paths;
% root = ['E:',filesep];
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

%% load brain metabolites
[~, ~, brainMetabolites_bis] = metabolite_load(subject_id);
MRS_ROI_names = {'dmPFC','aIns'};
n_MRS_ROIs = length(MRS_ROI_names);
brainMetabolite_names = fieldnames(brainMetabolites_bis.dmPFC);
n_brainMb = length(brainMetabolite_names);

%% load nutrition data
FFQ_struct = load([nutritionPath,'FFQ_metabolites_extracted_results.mat'],...
    'FFQ_results',...
    'energy_kcalPerWeek','metabolite_names',...
    'subject_id');
energy_kcalPerWeek = FFQ_struct.energy_kcalPerWeek;
foodNutrient_names = FFQ_struct.metabolite_names;
nFoodNutrients = length(foodNutrient_names);
% filter relevant subjects from nutrition data + correct for total caloric
% intake
FFQ_subject_id = FFQ_struct.subject_id;
for iN = 1:nFoodNutrients
    foodNutrient_nm = foodNutrient_names{iN};
    % initialize variable of interest
    food_nutrients.(foodNutrient_nm) = NaN(1,NS);
    
    % extract relevant subjects
    for iS = 1:NS
        sub_nm = subject_id{iS};
        sub_nutrition_idx = strcmp(sub_nm, FFQ_subject_id);
        % extract nutrition info and normalize by total energy intake
        if ismember(foodNutrient_nm,...
                {'energy_kcalPerWeek','calories_kcalPerWeek',...
                'GJEnergie_kJ_portion__PerWeek',...
                'GCALZBEnergie_kcalAvecFibres_portion__PerWeek',...
                'GJZBEnergie_kJAvecFibres_portion__PerWeek'})
            food_nutrients.(foodNutrient_nm)(iS) = FFQ_struct.FFQ_results.(foodNutrient_nm)(sub_nutrition_idx);
        else
            food_nutrients.(foodNutrient_nm)(iS) = FFQ_struct.FFQ_results.(foodNutrient_nm)(sub_nutrition_idx)./energy_kcalPerWeek(sub_nutrition_idx);
        end
    end % subject loop
end % loop through nutrients

%% perform the correlations
n_total_brainMbs = n_MRS_ROIs*n_brainMb;
[r_corr_mtrx, pval_mtrx] = deal(NaN(n_total_brainMbs, nFoodNutrients));
brain_mb_names = cell(n_total_brainMbs,1);
[foodNutrient_short_names,...
    food_nutrient_short_names2] = deal(cell(nFoodNutrients,1));
% loop through brain metabolites
for iMRS_ROI = 1:n_MRS_ROIs
    MRS_ROI_nm = MRS_ROI_names{iMRS_ROI};
    for iBrainMb = 1:n_brainMb
        brain_mb_nm = brainMetabolite_names{iBrainMb};
        jBrain_mb = iBrainMb + n_brainMb.*(iMRS_ROI - 1); % brain metabolite index for the matrix
        brain_mb_names{jBrain_mb} = [MRS_ROI_nm,' - ',brain_mb_nm];
        
        % loop through food nutrients
        for iN = 1:nFoodNutrients
            foodNutrient_nm = foodNutrient_names{iN};
            end_nm_idx = strfind(foodNutrient_nm,'_');
            foodNutrient_short_names{iN} = foodNutrient_nm(1:(end_nm_idx(1)-1));
            food_nutrient_short_names2{iN} = nutrient_nm_converter(foodNutrient_short_names{iN});
%             if ismember(foodNutrient_nm,{'energy_kcalPerWeek','calories_kcalPerWeek',...
%                 'GJEnergie_kJ_portion__PerWeek',...
%                 'GCALZBEnergie_kcalAvecFibres_portion__PerWeek',...
%                 'GJZBEnergie_kJAvecFibres_portion__PerWeek'})
%                 food_nutrient_short_names2{iN} = 'energy';
%             else
%                 food_nutrient_short_names2{iN} = nutrient_nm_converter(foodNutrient_names{iN});
%             end
            corr_nm = [MRS_ROI_nm,'_',brain_mb_nm,'_f_',foodNutrient_nm];
            % remove (or not) the outliers in each measure
            switch outlierF
                case 0
                    food_nutrient_var_tmp = food_nutrients.(foodNutrient_nm);
                    brainMetabolism_var_tmp = brainMetabolites_bis.(MRS_ROI_nm).(brain_mb_nm);
                case 1 % remove outliers based on mean +/- 3*SD
                    [~, ~, food_nutrient_var_tmp] = rmv_outliers_3sd(food_nutrients.(foodNutrient_nm));
                    [~, ~, brainMetabolism_var_tmp] = rmv_outliers_3sd(brainMetabolites_bis.(MRS_ROI_nm).(brain_mb_nm));
            end
            goodS_tmp = ~isnan(food_nutrient_var_tmp.*brainMetabolism_var_tmp);
            NS_goodS.(corr_nm) = sum(goodS_tmp);
            [r_corr.(corr_nm), pval.(corr_nm)] = corr(food_nutrients.(foodNutrient_nm)(goodS_tmp)', brainMetabolites_bis.(MRS_ROI_nm).(brain_mb_nm)(goodS_tmp)');
            
            % store significant results
            if pval.(corr_nm) < 0.05
                signif.p005.(corr_nm).r_corr = r_corr.(corr_nm);
                signif.p005.(corr_nm).pval = pval.(corr_nm);
                if pval.(corr_nm) < 0.001
                    signif.p0001.(corr_nm).r_corr = r_corr.(corr_nm);
                    signif.p0001.(corr_nm).pval = pval.(corr_nm);
                end % very significant p.value
            end % significant p.value
            
            % store data in the big matrix
            r_corr_mtrx(jBrain_mb, iN) = r_corr.(corr_nm);
            pval_mtrx(jBrain_mb, iN) = pval.(corr_nm);
        end % food nutrient
    end % brain metabolite loop
end % brain area loop

%% display resulting correlation matrix
% general figure parameters
pSize = 10;

% correlation range
corr_range = [-1 1];

apply_pval_threshold = false; % display everything even not significant results
pval_threshold = []; % no pvalue threshold
disp_signif_stars = true; % display stars upon the significant correlations

% display general correlation matrix
corr_plot(r_corr_mtrx, pval_mtrx,...
    corr_range, food_nutrient_short_names2, brain_mb_names, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);
end % function