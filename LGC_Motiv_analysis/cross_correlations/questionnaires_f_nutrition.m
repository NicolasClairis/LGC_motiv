% function[r_corr, pval] = questionnaires_f_nutrition()
% [r_corr, pval] = questionnaires_f_nutrition()
% questionnaires_f_nutrition will tests the correlations between the
% nutritional metabolites and the behavioral questionnaires (with a focus
% on the motivation-related questionnaires like MADRS, JPI and MPSTEFS)

%% define subjects
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
root = LGCM_root_paths;
% study
study_nm = 'study1';
switch root
    case ['E:',filesep] % lab pc
        gitPath = fullfile('C:','Users','clairis','Desktop','GitHub');
    case [filesep,filesep,fullfile('sv-nas1.rcp.epfl.ch','Sandi-lab','human_data_private','raw_data_subject'),filesep] % LGC server
        gitPath = fullfile('C:','Users','Nicolas Clairis','Documents','GitHub');
    case ['F:',filesep] % hard drive
        gitPath = fullfile('C:','Users','Nicolas Clairis','Documents','GitHub');
end
FFQ_resultsPath = [fullfile(gitPath, 'LGC_motiv','LGC_Motiv_results',...
    study_nm,'nutrition'),filesep];

%% load behavioral questionnaires (limit to motivation questionnaires on first instance)
questionnairesToCheck = {'JPI_RScore', 'MADRS_SCorrected',...
    'MPSTEFSPhysicalTraitScore','MPSTEFSMentalTraitScore',...
    'PunishmentScore','RewardScore',...
    'IPAQ','IPAQInactivity'};
nQuest = length(questionnairesToCheck);
[excelReadQuestionnairesFile, sub_CID_list] = load_questionnaires_data;
for iQuest = 1:nQuest
    quest_nm = questionnairesToCheck{iQuest};
    quest_data.(quest_nm) = NaN(1,NS);
end

%% load nutritional data
FFQ_struct = load([FFQ_resultsPath,'FFQ_metabolites_extracted_results.mat'],...
    'FFQ_results',...
    'energy_kcalPerWeek','metabolite_names',...
    'subject_id');
FFQ_results = FFQ_struct.FFQ_results;
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
    
    %% load behavioral questionnaires
    sub_quest_idx = strcmp(sub_CID_list, sub_nm);
    for iQuest = 1:nQuest
        quest_nm = questionnairesToCheck{iQuest};
        quest_data.(quest_nm)(iS) = excelReadQuestionnairesFile.(quest_nm)(sub_quest_idx);
    end % questionnaire loop
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
    for iQuest = 1:nQuest
        quest_nm = questionnairesToCheck{iQuest};
        corr_nm = [quest_nm,'_f_',nutrient_nm];
        goodSubs.(corr_nm) = (~isnan(quest_data.(quest_nm))).*(~isnan(food_nutrients.(nutrient_nm))) == 1;
        [r_corr.(corr_nm), betas.(corr_nm), pval.(corr_nm),...
            ~, food_nutrients_fit_sorted.(corr_nm),...
            quest_fit_foodNutrientSorted.(corr_nm)] = glm_package(food_nutrients.(nutrient_nm)(goodSubs.(corr_nm)),...
            quest_data.(quest_nm)(goodSubs.(corr_nm)),...
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
    end % questionnaire
end % nutrition loop

%% correlation and figure
% TO BE DONE EVENTUALLY (although so many tests that not it really makes
% sense, unless maybe a big correlation matrix). For now, we will just
% display the correlations that are significant (p<0.05 uncorrected for
% multiple comparisons)
figDisp = 1; % display figure?
pval_thresh = 0.05;

if figDisp == 1
    lWidth = 3;
    pSize = 25;
    black = 'k';
    orange = [254 75 3]./255;
    %% show results
    for iQuest = 1:nQuest
        quest_nm = questionnairesToCheck{iQuest};
        for iN = 1:nNutrients
            nutrient_nm = foodNutrient_names{iN};
            [nutrient_nm2] = convert_nutrient_nm_for_label(nutrient_nm);
            corr_nm = [quest_nm,'_f_',nutrient_nm];
            % define y.label for questionnaires
            switch quest_nm
                case 'JPI_RScore'
                    quest_nm_bis = 'JPI-R';
                case 'MADRS_SCorrected'
                    quest_nm_bis = 'MADRS-S';
                case 'MPSTEFSPhysicalTraitScore'
                    quest_nm_bis = 'MPSTEFS physical';
                case 'MPSTEFSMentalTraitScore'
                    quest_nm_bis = 'MPSTEFS mental';
                case 'PunishmentScore'
                    quest_nm_bis = 'PANAS P';
                case 'RewardScore'
                    quest_nm_bis = 'PANAS R';
                case 'IPAQ'
                    quest_nm_bis = 'IPAQ activity';
                case 'IPAQInactivity'
                    quest_nm_bis = 'IPAQ inactivity';
            end

            %% figure
            if pval.(corr_nm)(2) < pval_thresh
                fig;
                scat_hdl = scatter(food_nutrients.(nutrient_nm)(goodSubs.(corr_nm)),...
                    quest_data.(quest_nm)(goodSubs.(corr_nm)));
                plot_hdl = plot(food_nutrients_fit_sorted.(corr_nm), quest_fit_foodNutrientSorted.(corr_nm));
                scat_hdl.LineWidth = lWidth;
                scat_hdl.MarkerEdgeColor = black;
                plot_hdl.Color = orange;
                plot_hdl.LineStyle = '--';
                xlabel(nutrient_nm2);
                ylabel(quest_nm_bis);
                legend_size(pSize);
                place_r_and_pval(r_corr.(corr_nm),pval.(corr_nm)(2));
            end % p.value threshold to decide which figure to display
        end % nutrient loop
    end % questionnaire loop
end % fig disp

% end % function