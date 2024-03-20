% function[] = FFQ_nutrition_f_sex()
% FFQ_nutrition_f_sex will compare nutritional components extracted with
% FFQ between the two sexes.
%
% See also
% C:\Users\clairis\Desktop\GitHub\LGC_motiv\LGC_Motiv_analysis\nutrition\BLS-Dateiaufbau.xlsx
% for the meaning of the acronyms based on the German food database

%% subject selection
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

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

%% load nutrition data
FFQ_struct = load([nutritionPath,'FFQ_metabolites_extracted_results.mat']);
energy_kcalPerWeek = FFQ_struct.energy_kcalPerWeek;
foodNutrient_names = FFQ_struct.metabolite_names;
nNutrients = length(foodNutrient_names);
FFQ_subject_id = FFQ_struct.subject_id;

%% extract all nutrients independently for each sex
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
    
    % initialize variable of interest
    for iF = 1:nNutrients
        fd_nm = foodNutrient_names{iF};
        [nutrition.(fd_nm).(sex_nm),...
            nutrition_f_calories.(fd_nm).(sex_nm)] = deal(NaN(1,NS));
    end % loop over food items
    
    for iS = 1:NS
        sub_nm = subject_id{iS};
        sub_idx = strcmp(FFQ_subject_id, sub_nm);
        % extract total calories for current subject
        total_caloric_intake_tmp = FFQ_struct.FFQ_results.calories_kcalPerWeek(sub_idx);
        % extract all nutrients (raw and normalized by total caloric
        % intake)
        for iF = 1:nNutrients
            fd_nm = foodNutrient_names{iF};
            nutrition.(fd_nm).(sex_nm)(iS) = FFQ_struct.FFQ_results.(fd_nm)(sub_idx);
            nutrition_f_calories.(fd_nm).(sex_nm)(iS) = nutrition.(fd_nm).(sex_nm)(iS)./total_caloric_intake_tmp;
        end % loop through nutrients
    end % subject loop
    
end % sex loop

%% compare male vs female nutritional intake
corr_method = 'bonferroni';
corr_method_nm = ['corr_',corr_method];
p_unc_raw_all = [];
p_unc_calNorm_all = [];
food_names_bis_unc_raw = {};
food_names_bis_unc_calNorm_all = {};
for iF = 1:nNutrients
    fd_nm = foodNutrient_names{iF};
    % uncorrected tests
    % raw nutrient intake
    [~,pval.uncorr.raw.(fd_nm)] = ttest2(nutrition.(fd_nm).male, nutrition.(fd_nm).female);
    
    if ~ismember(fd_nm,{'calories_kcalPerWeek',...
            'GJEnergie_kJ_portion__PerWeek',...
            'GCALZBEnergie_kcalAvecFibres_portion__PerWeek',...
            'GJZBEnergie_kJAvecFibres_portion__PerWeek'}) % pointless to include calories/calories
        % nutrient intake corrected by total caloric intake
        [~,pval.uncorr.norm_by_caloric_intake.(fd_nm)] = ttest2(nutrition_f_calories.(fd_nm).male, nutrition_f_calories.(fd_nm).female);
        
        % pool all p.value together for correction
        p_unc_raw_all = [p_unc_raw_all, pval.uncorr.raw.(fd_nm)];
        p_unc_calNorm_all = [p_unc_calNorm_all, pval.uncorr.norm_by_caloric_intake.(fd_nm)];
        % extract name of corresponding food item
        food_names_bis_unc_raw = [food_names_bis_unc_raw; fd_nm];
        food_names_bis_unc_calNorm_all = [food_names_bis_unc_calNorm_all; fd_nm];
    elseif strcmp(fd_nm,'calories_kcalPerWeek') % for raw intake include one measure of total caloric intake, but not all of them as it doesn't make sense to correct for the 4 measures
        p_unc_raw_all = [p_unc_raw_all, pval.uncorr.raw.(fd_nm)];
        food_names_bis_unc_raw = [food_names_bis_unc_raw; fd_nm];
    end
    
    % extract results that are significantly different
    if pval.uncorr.raw.(fd_nm) < 0.05
        signif.uncorr.raw.(fd_nm).pval = pval.uncorr.raw.(fd_nm);
        [signif.uncorr.raw.(fd_nm).male_m,...
            signif.uncorr.raw.(fd_nm).male_sem] = mean_sem_sd(nutrition.(fd_nm).male,2);
        [signif.uncorr.raw.(fd_nm).female_m,...
            signif.uncorr.raw.(fd_nm).female_sem] = mean_sem_sd(nutrition.(fd_nm).female,2);
    end
    % same for caloric intake normalized data
    if ~ismember(fd_nm,{'calories_kcalPerWeek',...
            'GJEnergie_kJ_portion__PerWeek',...
            'GCALZBEnergie_kcalAvecFibres_portion__PerWeek',...
            'GJZBEnergie_kJAvecFibres_portion__PerWeek'}) % pointless to include calories/calories
        if pval.uncorr.norm_by_caloric_intake.(fd_nm) < 0.05
            signif.uncorr.norm_by_caloric_intake.(fd_nm).pval = pval.uncorr.norm_by_caloric_intake.(fd_nm);
            [signif.uncorr.norm_by_caloric_intake.(fd_nm).male_m,...
                signif.uncorr.norm_by_caloric_intake.(fd_nm).male_sem] = mean_sem_sd(nutrition_f_calories.(fd_nm).male,2);
            [signif.uncorr.norm_by_caloric_intake.(fd_nm).female_m,...
                signif.uncorr.norm_by_caloric_intake.(fd_nm).female_sem] = mean_sem_sd(nutrition_f_calories.(fd_nm).female,2);
        end
    end
end % loop over nutrients

%% tests corrected for multiple comparisons
% correct pvalues
pval_raw_corr_all = pval_adjust(p_unc_raw_all,corr_method);
pval_calNorm_corr_all = pval_adjust(p_unc_calNorm_all,corr_method);
% extract pvalue for each nutrient
for iF = 1:length(food_names_bis_unc_raw)
    fd_nm = food_names_bis_unc_raw{iF};
    pval.(corr_method_nm).raw.(fd_nm) = pval_raw_corr_all(iF);
    
    % store data if significant
    if pval.(corr_method_nm).raw.(fd_nm) < 0.05
        signif.(corr_method_nm).raw.(fd_nm).pval = pval.(corr_method_nm).raw.(fd_nm);
        [signif.(corr_method_nm).raw.(fd_nm).male_m,...
            signif.(corr_method_nm).raw.(fd_nm).male_sem] = mean_sem_sd(nutrition.(fd_nm).male,2);
        [signif.(corr_method_nm).raw.(fd_nm).female_m,...
            signif.(corr_method_nm).raw.(fd_nm).female_sem] = mean_sem_sd(nutrition.(fd_nm).female,2);
    end
end

for iF = 1:length(food_names_bis_unc_calNorm_all)
    fd_nm = food_names_bis_unc_calNorm_all{iF};
    pval.(corr_method_nm).norm_by_caloric_intake.(fd_nm) = pval_calNorm_corr_all(iF);
    
    % store data if significant
    if pval.(corr_method_nm).norm_by_caloric_intake.(fd_nm) < 0.05
        signif.(corr_method_nm).norm_by_caloric_intake.(fd_nm).pval = pval.(corr_method_nm).norm_by_caloric_intake.(fd_nm);
        [signif.(corr_method_nm).norm_by_caloric_intake.(fd_nm).male_m,...
            signif.(corr_method_nm).norm_by_caloric_intake.(fd_nm).male_sem] = mean_sem_sd(nutrition.(fd_nm).male,2);
        [signif.(corr_method_nm).norm_by_caloric_intake.(fd_nm).female_m,...
            signif.(corr_method_nm).norm_by_caloric_intake.(fd_nm).female_sem] = mean_sem_sd(nutrition.(fd_nm).female,2);
    end
end

%% create one figure hilighting the data significantly different
% colour code
[~, ~, col] = general_fig_prm;
female_col = col.red;
male_col = col.blue_dark;

% uncorrected + raw data
fig;
jF = 0;
fd_label = {};
for iF = 1:nNutrients
    fd_nm = foodNutrient_names{iF};
    if isfield(signif.uncorr.raw, fd_nm) % if nutrient is significant
        jF = jF + 1;
        ok_males = ~isnan(nutrition.(fd_nm).male);
        Violin({nutrition.(fd_nm).male(ok_males)},jF,'ViolinColor',{male_col});
        
        jF = jF + 1;
        ok_females = ~isnan(nutrition.(fd_nm).female);
        Violin({nutrition.(fd_nm).female(ok_females)},jF,'ViolinColor',{female_col});
        
        %% add p.value
        [l_hdl, star_hdl] = add_pval_comparison(nutrition.(fd_nm).male,...
            nutrition.(fd_nm).female,...
            signif.uncorr.raw.(fd_nm).pval, jF-1, jF, '');
        
        % extract name for labels
        fd_nm_bis = fd_nm(1:4);
        fd_label = [fd_label, fd_nm_bis];
    end
end % nutrient loop

ylabel('Nutrient per week');
xticks(1.5:2:length(fd_label)*2);
xticklabels(fd_label);

% end % function