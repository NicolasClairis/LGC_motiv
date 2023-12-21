% function[r_corr, pval] = blood_f_nutrition()
% [r_corr, pval] = blood_f_nutrition()
% Try correlations between nutrition metabolites and blood metabolites (in
% particular Lactate which could correlate with some energy-related food
% nutrients).

%% define subjects
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
root = LGCM_root_paths;
% study
study_nm = 'study1';
switch root
    case 'E:'
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
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];

    %% load blood
    sub_blood_idx = strcmp(Lac_struct.CID, sub_nm_bis);
    Lac(iS) = Lac_struct.Lac(sub_blood_idx);
    %% load nutrition data
    sub_nutrition_idx = strcmp(sub_nm, FFQ_subject_id);
    % extract nutrition info
    for iN = 1:nNutrients
        nutrient_nm = foodNutrient_names{iN};
        food_nutrients.(nutrient_nm)(iS) = food_nutrients_allData.(nutrient_nm)(sub_nutrition_idx);
    end
end % subject loop

%% test correlations
for iN = 1:nNutrients
    nutrient_nm = foodNutrient_names{iN};
    corr_nm = ['Lac_f_',nutrient_nm];
    goodSubs.(corr_nm) = (~isnan(Lac)).*(~isnan(food_nutrients.(nutrient_nm))) == 1;
    [r_corr.(corr_nm), betas.(corr_nm), pval.(corr_nm),...
        ~, food_nutrients_fit_sorted.(corr_nm),...
        Lac_fit_foodNutrientSorted.(corr_nm)] = glm_package(food_nutrients.(nutrient_nm)(goodSubs.(corr_nm)),...
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
end % nutrition loop

%% correlation and figure
% TO BE DONE EVENTUALLY (although so many tests that not it really makes
% sense, unless maybe something Arthur-like with the correlation matrix)

% if figDisp == 1
%     lWidth = 3;
%     pSize = 25;
%     black = 'k';
%     orange = [254 75 3]./255;
%     %% show results
%     for iQuest = 1:nQuest
%         quest_nm = questToCheck{iQuest};
%         fig1 = fig; j_fig1 = 0;
%         fig2 = fig; j_fig2 = 0;
%         for iN = 1:n_BloodPrm
%             bloodMb_nm = bloodMb_names{iN};
%             corr_nm = [quest_nm,'_f_',bloodMb_nm];
%             % define y.label for questionnaires
%             switch quest_nm
%                 case 'JPI_RScore'
%                     quest_nm_bis = 'JPI-R';
%                 case 'MADRS_SCorrected'
%                     quest_nm_bis = 'MADRS-S';
%                 case 'MPSTEFSPhysicalTraitScore'
%                     quest_nm_bis = 'MPSTEFS physical';
%                 case 'MPSTEFSMentalTraitScore'
%                     quest_nm_bis = 'MPSTEFS mental';
%                 case 'PunishmentScore'
%                     quest_nm_bis = 'PANAS P';
%                 case 'RewardScore'
%                     quest_nm_bis = 'PANAS R';
%                 case 'IPAQ'
%                     quest_nm_bis = 'IPAQ activity';
%                 case 'IPAQInactivity'
%                     quest_nm_bis = 'IPAQ inactivity';
%             end
% 
%             switch bloodMb_nm
%                 case {'Nam','NMN','NR','NAD',...
%                         'NADH','NADP','NADPH','MeNam',...
%                         'MeXPY'}
%                     figure(fig1);
%                     j_fig1 = j_fig1 + 1;
%                     subplot(3,4,j_fig1);
%                 case {'NAD_div_NADH',...
%                         'NADP_div_NADPH',...
%                         'total_NAD_precursors',...
%                         'total_NAD',...
%                         'total_NAD_with_precursors',...
%                         'total_NAD_with_byproducts',...
%                         'total_NAD_byproducts'}
%                     figure(fig2);
%                     j_fig2 = j_fig2 + 1;
%                     subplot(3,3,j_fig2);
%             end
%             %% figure
%             hold on;
%             scat_hdl = scatter(bloodMb.(bloodMb_nm)(goodSubs.(corr_nm)),...
%                 quest_data.(quest_nm)(goodSubs.(corr_nm)));
%             plot_hdl = plot(bloodMb_sort.(corr_nm), quest_fit.(corr_nm));
%             scat_hdl.LineWidth = lWidth;
%             scat_hdl.MarkerEdgeColor = black;
%             plot_hdl.Color = orange;
%             plot_hdl.LineStyle = '--';
%             [blood_labelname] = blood_label(bloodMb_nm);
%             xlabel(blood_labelname);
%             ylabel(quest_nm_bis);
%             legend_size(pSize);
%         end % metabolite loop
%     end % questionnaire loop
% end % fig disp

% end % function