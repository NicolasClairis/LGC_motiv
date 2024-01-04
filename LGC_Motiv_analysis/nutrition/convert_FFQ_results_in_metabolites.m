% convert FFQ results in metabolites

%% working directory
% gitPath = fullfile('C:','Users','clairis','Desktop','GitHub','LGC_motiv');
gitPath = fullfile('C:','Users','Nicolas Clairis','Documents','GitHub','LGC_motiv');
resultsPath = [fullfile(gitPath,'LGC_Motiv_results','study1','nutrition'),filesep];
FFQ_conversion_path = [fullfile(gitPath,'LGC_Motiv_analysis','nutrition'),filesep];

%% load FFQ conversion table
food_to_metabolites_conversion_table = readtable([FFQ_conversion_path,'taux2.xlsx'],...
    'Sheet','AFTER CONVERSION');
food_questions = food_to_metabolites_conversion_table.Question;
nFoodQuestions = length(food_questions);

%% load frequency (per week) for each food item for each individual
% load results
FFQ_data = readtable([resultsPath,'FFQ.xlsx']);
NS = size(FFQ_data,1);
subject_id = FFQ_data.CID;
% convert frequency from strings to numerical values
for iQ1 = 1:nFoodQuestions
    foodQuestion_nm = food_questions{iQ1};
    result_freqPerWeek_perFoodItem.(foodQuestion_nm) = NaN(1,NS);
    for iS = 1:NS
        sub_nm = FFQ_data.CID{iS};
         result_foodItem_freq_nm = FFQ_data.(foodQuestion_nm){iS};
         % convert frequency label into number
         switch result_foodItem_freq_nm
             case 'Jamais ou <1/Mois'
                 result_foodItem_freq = 0;
             case '1-3x par Mois'
                 % 1-3x/month ~= 2times/Month
                 % 1 month ~= 4.3 weeks (based on this https://www.cuemath.com/questions/what-is-the-average-number-of-weeks-in-a-month/)
                 % => 2/4.3 times/week
                 result_foodItem_freq = 2/4.3;
             case '1x par Sem.'
                 result_foodItem_freq = 1;
             case '2-4x par Sem.'
                 % ~* 3 times/week
                 result_foodItem_freq = 3;
             case '5-6x par Sem.'
                 % 5-6 times per week ~= 5.5 times/week
                 result_foodItem_freq = 5.5;
             case '1x par Jour'
                 % 7 days/week
                 result_foodItem_freq = 7;
             case '2-3x par Jour'
                 % ~=2.5 times/day*7days/week
                 result_foodItem_freq = 2.5*7;
             case '4-5x par Jour'
                 % ~= 4.5 times/day*7days/week
                 result_foodItem_freq = 4.5*7;
             case '6x par Jour'
                 % 6 times*7days = 42 times/week
                 result_foodItem_freq = 6*7;
         end
         result_freqPerWeek_perFoodItem.(foodQuestion_nm)(iS) = result_foodItem_freq;
    end % subject loop
end % question loop

%% convert data to know the metabolites coming from nutrition and store it in a resulting table and structure

% first column of the new table = CID of the participants
FFQ_results = table(subject_id);

% extract name of all metabolites in the list
metabolite_names_v0 = food_to_metabolites_conversion_table.Properties.VariableNames(5:end);
n_mb = length(metabolite_names_v0);

% loop through each metabolite to extract and convert the FFQ into a corresponding
% concentration of metabolite
for i_mb = 1:n_mb
    mb_nm = metabolite_names_v0{i_mb};
    mb_perSub_perWeek = zeros(NS,1);
    
    % loop through subjects
    for iS = 1:NS
        
        % compute current metabolite dose for the corresponding subject
        for iQ2 = 1:nFoodQuestions
            foodItem_nm2 = food_questions{iQ2};
            % food item => metabolite conversion value
            foodItem_to_mb_conversionValue = food_to_metabolites_conversion_table.(mb_nm)(iQ2);
            % frequency food item is taken by the subject
            foodItem_freqPerWeek = result_freqPerWeek_perFoodItem.(foodItem_nm2)(iS);
            mb_perSub_perWeek(iS) = mb_perSub_perWeek(iS) + foodItem_freqPerWeek.*foodItem_to_mb_conversionValue;
        end
    end % subject loop
    
    %% add result in the table
    % for the most important metabolites, hilight the name of the
    % metabolite
    switch mb_nm
        case 'GCALEnergie_kcal_portion_' % calories per week
            FFQ_results.calories_kcalPerWeek = mb_perSub_perWeek;
            energy_kcalPerWeek = mb_perSub_perWeek;
        case 'ZF_graisses__mg_portion_' % total fat
            FFQ_results.fat_mgPerWeek = mb_perSub_perWeek;
            fat_mgPerWeek = mb_perSub_perWeek;
        case 'KMD_totalSucre__mg_portion_' % sugars
            FFQ_results.sugars_mgPerWeek = mb_perSub_perWeek;
            sugars_mgPerWeek = mb_perSub_perWeek;
        case 'ZE_prot_ines__mg_portion_' % proteins
            FFQ_results.proteins_mgPerWeek = mb_perSub_perWeek;
            proteins_mgPerWeek = mb_perSub_perWeek;
        case 'EGLY_Glycine__mg_portion_' % Glycine
            FFQ_results.Gly_mgPerWeek = mb_perSub_perWeek;
            Gly_mgPerWeek = mb_perSub_perWeek;
        case 'EGLU_acideGlutamatergique__mg_portion_' % Glutamate
            FFQ_results.Glu_mgPerWeek = mb_perSub_perWeek;
            Glu_mgPerWeek = mb_perSub_perWeek;
        case 'EASP_acideAspartique__mg_portion_' % Aspartate
            FFQ_results.Asp_mgPerWeek = mb_perSub_perWeek;
            Asp_mgPerWeek = mb_perSub_perWeek;
        case 'ETRP_Tryptophane__mg_portion_' % Tryptophan
            FFQ_results.Trp_mgPerWeek = mb_perSub_perWeek;
            Trp_mgPerweek = mb_perSub_perWeek;
        case 'ECYS_Cyst_ine__mg_portion_' % cystein
            FFQ_results.Cys_mgPerWeek = mb_perSub_perWeek;
            Cys_mgPerWeek = mb_perSub_perWeek;
        case 'VB3_vitamineB3_Niacine_acideNicotinique___g_portion_' % niacin
            FFQ_results.NA_ugPerWeek = mb_perSub_perWeek;
            NA_ugPerWeek = mb_perSub_perWeek;
        case 'VB3A_vitamineB3__quivalentsDeNiacine___g_portion_' % niacin equivalents
            FFQ_results.NE_ugPerWeek = mb_perSub_perWeek;
            NE_ugPerweek = mb_perSub_perWeek;
        case 'VB9G_vitamineB9_AcidesFoliques___g_portion_' % folic acid
            FFQ_results.FolicAcid_mgPerWeek = mb_perSub_perWeek;
            folicAcid_mgPerWeek = mb_perSub_perWeek;
        case 'FO3_omega3__mg_portion_' % omega 3
            FFQ_results.omega3_mgPerWeek = mb_perSub_perWeek;
            omega3_mgPerWeek = mb_perSub_perWeek;
        case 'ZK_glucidesAssimilables__mg_portion_' % glucides
            FFQ_results.glucids_mgPerWeek = mb_perSub_perWeek;
            Glucids_mgPerWeek = mb_perSub_perWeek;
        case 'MMG_Magn_sium__mg_portion_' % magnesium
            FFQ_results.Mg_mgPerWeek = mb_perSub_perWeek;
            Magnesium_mgPerWeek = mb_perSub_perWeek;
        case 'KMT_Glucose__mg_portion_' % glucose
            FFQ_results.Glc_mgPerWeek = mb_perSub_perWeek;
            Glc_mgPerWeek = mb_perSub_perWeek;
        case 'KDL_lactose__mg_portion_' % lactose
            FFQ_results.Lactose_mgPerWeek = mb_perSub_perWeek;
            Lactose_mgPerWeek = mb_perSub_perWeek;
        otherwise
           FFQ_results.([mb_nm,'_PerWeek']) = mb_perSub_perWeek;
    end
end % metabolite loop
% re-extract labels
metabolite_names = fieldnames(FFQ_results);
% remove subjects and non-nutrient labels
metabolite_names(ismember(metabolite_names,{'subject_id','Properties','Row','Variables'})) = [];
% save all relevant data in mat format
save([resultsPath,'FFQ_metabolites_extracted_results.mat'],'FFQ_results',...
    'energy_kcalPerWeek','fat_mgPerWeek','sugars_mgPerWeek',...
    'proteins_mgPerWeek','Gly_mgPerWeek','Glu_mgPerWeek','Asp_mgPerWeek',...
    'Trp_mgPerweek','Cys_mgPerWeek', 'NA_ugPerWeek', 'NE_ugPerweek',...
    'folicAcid_mgPerWeek', 'omega3_mgPerWeek', 'Glucids_mgPerWeek',...
    'Magnesium_mgPerWeek', 'Glc_mgPerWeek', 'Lactose_mgPerWeek',...
    'Lactose_mgPerWeek',...
    'FFQ_data','food_to_metabolites_conversion_table',...
    'subject_id','metabolite_names_v0','metabolite_names',...
    'food_questions','n_mb','nFoodQuestions','NS');

% save data as excel file
% careful git path may be blocked => you may need to write it elsewhere first and transfer it afterwards
excelfilename = [resultsPath,'FFQ_metabolites_extracted_results.xlsx']; 
writetable(FFQ_results, excelfilename,'Sheet',1)