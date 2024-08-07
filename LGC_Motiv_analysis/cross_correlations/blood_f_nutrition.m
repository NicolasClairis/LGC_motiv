% function[r_corr, pval] = blood_f_nutrition(fig_disp)
% [r_corr, pval] = blood_f_nutrition(fig_disp)
% Try correlations between nutrition metabolites and blood metabolites (in
% particular Lactate which could correlate with some energy-related food
% nutrients).
%
% INPUTS
% fig_disp: display figure (1) or not (0)?

%% define default inputs
if ~exist('fig_disp','var') || isempty(fig_disp)
    fig_disp = 1;
end

%% define subjects
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
root = LGCM_root_paths;
% study
study_nm = 'study1';
switch root
    case ['E:',filesep]
        gitPath = fullfile('C:','Users','clairis','Desktop','GitHub');
    case [filesep,filesep,fullfile('sv-nas1.rcp.epfl.ch','Sandi-lab','human_data_private','raw_data_subject'),filesep]
        gitPath = fullfile('C:','Users','Nicolas Clairis','Documents','GitHub');
end
FFQ_resultsPath = [fullfile(gitPath, 'LGC_motiv','LGC_Motiv_results',...
    study_nm,'nutrition'),filesep];

%% load blood data (limit to Lactate on first instance)
[Lac_struct] = load_plasma_Lac;
% initialize variable for extraction
Lac = NaN(1,NS);

%% load nutritional data
FFQ_struct = load([FFQ_resultsPath,'FFQ_metabolites_extracted_results.mat'],...
    'FFQ_results',...
    'energy_kcalPerWeek','metabolite_names',...
    'subject_id');
energy_kcalPerWeek = FFQ_struct.energy_kcalPerWeek;
foodNutrient_names = FFQ_struct.metabolite_names;
nNutrients = length(foodNutrient_names);
FFQ_subject_id = FFQ_struct.subject_id;

% correct for caloric intake or not?
caloric_corr = questdlg('Correct food nutrients for total caloric intake?',...
    'Calories correction','yes','no','yes');
% extract data (and correct for total caloric intake if requested)
for iN = 1:nNutrients
    nutrient_nm = foodNutrient_names{iN};
    % prepare var of interest
    food_nutrients.(nutrient_nm) = NaN(1,NS);
    % normalize (or not) food data by energy levels
    switch caloric_corr
        case 'yes'
            food_nutrients_allData.(nutrient_nm) = FFQ_struct.FFQ_results.(nutrient_nm)./energy_kcalPerWeek;
        case 'no'
            food_nutrients_allData.(nutrient_nm) = FFQ_struct.FFQ_results.(nutrient_nm);
    end
end

%% loop through subjects to extract Lactate and FFQ info for each subject
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];

    %% load blood
    sub_blood_idx = strcmp(Lac_struct.CID, sub_nm_bis);
    Lac(iS) = Lac_struct.Lac(sub_blood_idx)./1000;
    %% load nutrition data
    sub_nutrition_idx = strcmp(sub_nm, FFQ_subject_id);
    % extract nutrition info
    for iN = 1:nNutrients
        nutrient_nm = foodNutrient_names{iN};
        food_nutrients.(nutrient_nm)(iS) = food_nutrients_allData.(nutrient_nm)(sub_nutrition_idx);
    end
end % subject loop

%% initialize correlation matrix data
[corr_mtrx, corr_mtrx_pval] = deal(NaN(1, nNutrients));
[foodNutrient_short_names,...
    foodNutrient_short_names2] = deal(cell(1, nNutrients));

%% test correlations
for iN = 1:nNutrients
    nutrient_nm = foodNutrient_names{iN};
    
    if strcmp(caloric_corr,'yes') && strcmp(nutrient_nm,'energy_kcalPerWeek')
    else
        % create short name by removing anything after '_'
        end_nm_idx = strfind(nutrient_nm,'_');
        foodNutrient_short_names{iN} = nutrient_nm(1:(end_nm_idx(1)-1));
        foodNutrient_short_names2{iN} = nutrient_nm_converter(foodNutrient_short_names{iN});
        corr_nm = ['Lac_f_',nutrient_nm];
        goodSubs.(corr_nm) = (~isnan(Lac)).*(~isnan(food_nutrients.(nutrient_nm))) == 1;
        [r_corr.(corr_nm), betas.(corr_nm), pval.(corr_nm),...
            ~, food_nutrients_fit_sorted.(corr_nm),...
            Lac_fit_foodNutrientSorted.(corr_nm),...
            pval_z.(corr_nm)] = glm_package(food_nutrients.(nutrient_nm)(goodSubs.(corr_nm)),...
            Lac(goodSubs.(corr_nm)),...
            'normal', 'on');
        pval_correl_tmp = pval.(corr_nm)(2);
        % store significant p.values for slope
        if pval_correl_tmp < 0.05
            signif.signif.pval.(corr_nm) = pval_correl_tmp;
            signif.signif.r_corr.(corr_nm) = r_corr.(corr_nm);
        elseif pval_correl_tmp > 0.05 && pval_correl_tmp < 0.1
            signif.almostSignif.pval.(corr_nm) = pval_correl_tmp;
            signif.almostSignif.r_corr.(corr_nm) = r_corr.(corr_nm);
        end
        
        % store all in a big correlation matrix (actually a vector)
        corr_mtrx(iN) = r_corr.(corr_nm);
        corr_mtrx_pval(iN) = pval_z.(corr_nm);
    end % filter caloric intake when variables corrected by caloric intake (just a series of 1)
end % nutrition loop

%% correlation and figure
% display global correlation matrix + significant correlations
if fig_disp == 1
    lWidth = 3;
    pSize = 25;
    black = 'k';
    orange = [254 75 3]./255;
    color_range_choices = redblue(45);
    
    % correlation range
    corr_range = [-1 1];
    apply_pval_threshold = 0;
    pval_threshold = [];
    disp_signif_stars = 2; % also display almost significant correlations
    
    %% correlation matrix
    corr_plot(corr_mtrx, corr_mtrx_pval,...
        corr_range, foodNutrient_short_names2, 'Lac','', '',...
        apply_pval_threshold, pval_threshold, disp_signif_stars);
    legend_size(10);
    
    %% show significant results
    for iN = 1:nNutrients
        nutrient_nm = foodNutrient_names{iN};
        short_nutrient_nm = foodNutrient_short_names2{iN};
        corr_nm = ['Lac_f_',nutrient_nm];
        if ~isempty(signif.signif.pval)
            signif_corrs = fieldnames(signif.signif.pval);
        else
            signif_corrs = {''};
        end
        if ~isempty(signif.almostSignif.pval)
            almostSignif_corrs = fieldnames(signif.almostSignif.pval);
        else
            almostSignif_corrs = {''};
        end
        
        if ismember(corr_nm,[almostSignif_corrs; signif_corrs])
            fig;
            scat_hdl = scatter(food_nutrients.(nutrient_nm), Lac);
            fit_hdl = plot(food_nutrients_fit_sorted.(corr_nm),...
                Lac_fit_foodNutrientSorted.(corr_nm));
            scat_hdl_upgrade(scat_hdl);
            fit_hdl_upgrade(fit_hdl);
            place_r_and_pval(corr_mtrx(iN), corr_mtrx_pval(iN));
            if (length(nutrient_nm) > 10 && strcmp(nutrient_nm(end-10:end),'kcalPerWeek')) ||...
                    (length(nutrient_nm) > 30 && strcmp(nutrient_nm(end-30:end),'kcalAvecFibres_portion__PerWeek'))
                xlabel([short_nutrient_nm,' (kCal/week)']);
            elseif (length(nutrient_nm) > 18 && strcmp(nutrient_nm(end-18:end),'kJ_portion__PerWeek')) ||...
                    (length(nutrient_nm) > 28 && strcmp(nutrient_nm(end-28:end),'kJAvecFibres_portion__PerWeek'))
                xlabel([short_nutrient_nm,' (kJ/week)']);
            elseif length(nutrient_nm) > 17 && (strcmp(nutrient_nm(end-17:end),'g_portion__PerWeek'))
                xlabel([short_nutrient_nm,' (g/week)']);
            elseif (length(nutrient_nm) > 18 && (strcmp(nutrient_nm(end-18:end),'mg_portion__PerWeek'))) ||...
                    (length(nutrient_nm) > 8 && (strcmp(nutrient_nm(end-8:end),'mgPerWeek')))
                xlabel([short_nutrient_nm,' (mg/week)']);
            elseif length(nutrient_nm) > 8 && (strcmp(nutrient_nm(end-8:end),'ugPerWeek'))
                xlabel([short_nutrient_nm,' (μg/week)']);
            else
                error(['case ',nutrient_nm,' to plan']);
            end
            ylabel('Plasma Lactate (mM)');
            legend_size(pSize);
        end % filter signif or almost signif correlations
    end % questionnaire loop
end % fig disp

% end % function