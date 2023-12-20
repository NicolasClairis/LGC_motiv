% convert FFQ results in metabolites

%% working directory
resultsPath = [fullfile('C:','Users','clairis','Desktop','GitHub',...
    'LGC_motiv','LGC_Motiv_results','study1','nutrition'),filesep];
FFQ_conversion_path = [fullfile('C:','Users','clairis','Desktop','GitHub',...
    'LGC_motiv','LGC_Motiv_analysis','nutrition'),filesep];

%% load FFQ conversion table
food_to_metabolites_conversion_table = readtable([FFQ_conversion_path,'taux2.xlsx'],...
    'Sheet','AFTER CONVERSION');
food_questions = food_to_metabolites_conversion_table.Question;
nFoodQuestions = length(food_questions);

%% load frequency (per week) for each food item for each individual
% load results
FFQ_data = readtable([resultsPath,'FFQ.xlsx']);
NS = size(FFQ_data,1);
CID = FFQ_data.CID;
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
                 % ~* 2.5 times/week
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
FFQ_results = table(CID);

% extract name of all metabolites in the list
metabolite_names = food_to_metabolites_conversion_table.Properties.VariableNames(5:end);
n_mb = length(metabolite_names);

% loop through each metabolite to extract and convert the FFQ into a corresponding
% concentration of metabolite
for i_mb = 1:n_mb
    mb_nm = metabolite_names{i_mb};
    mb_perSub_perWeek = zeros(1,NS);
    
    % loop through subjects
    for iS = 1:NS
        
        % compute current metabolite dose for the corresponding subject
        for iQ2 = 1:nFoodQuestions
            foodItem_nm2 = food_questions{iQ2};
            % food item => metabolite conversion value
            foodItem_to_mb_conversionValue = food_to_metabolites_conversion_table.(mb_nm)(iQ2);
            % frequency food item is taken by the subject
            foodItem_freqPerWeek = result_freqPerWeek_perFoodItem.(foodItem_nm2)(iS);
            mb_perSub_perWeek(iS) = mb_perSub(iS) + foodItem_freqPerWeek.*foodItem_to_mb_conversionValue;
        end
    end % subject loop
    
    %% add result in the table
    % for the most important metabolites, hilight the name of the
    % metabolite
    switch mb_nm
        case 'GCALEnergie_kcal_portion_' % calories per week
            FFQ_results.calories_kcalPerWeek = mb_perSub_perWeek;
        case % total fat
            FFQ_results.fat_mgPerWeek = mb_perSub_perWeek;
        case % sugars
            FFQ_results.sugars_mgPerWeek = mb_perSub_perWeek;
        case % proteins
            FFQ_results.proteins_mgPerWeek = mb_perSub_perWeek;
        case % Glycine
            FFQ_results.Gly_mgPerWeek = mb_perSub_perWeek; 
        case % Glutamate
            FFQ_results.Glu_mgPerWeek = mb_perSub_perWeek; 
        case % niacin
            FFQ_results.NA_ugPerWeek = mb_perSub_perWeek; 
        case % niacin equivalents
            FFQ_results.NE_ugPerWeek = mb_perSub_perWeek; 
        case % Tryptophan
            FFQ_results.Trp_mgPerWeek = mb_perSub_perWeek; 
        case % omega 3
            FFQ_results.omega3_mgPerWeek = mb_perSub_perWeek; 
        case % cystein
            FFQ_results.Cys_mgPerWeek = mb_perSub_perWeek; 
        case % glucides
            FFQ_results.glucids_mgPerWeek = mb_perSub_perWeek; 
        case % magnesium
            FFQ_results.Mg_ugPerWeek = mb_perSub_perWeek; 
        case % glucose
            FFQ_results.Glc_mgPerWeek = mb_perSub_perWeek; 
        case 'EASP_acideAspartique__mg_portion_' % Aspartate
            FFQ_results.Asp_mgPerWeek = mb_perSub_perWeek; 
        otherwise
           FFQ_results.([mb_nm,'_PerWeek') = mb_perSub_perWeek; 
    end
end % metabolite loop

