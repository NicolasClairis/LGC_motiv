% function[] = blood_f_nutrition_metabolites()
% blood_f_nutrition_metabolites will check whether plasma levels of
% particular metabolites correspond to the evaluated nutritional intake
% based on the FFQ questionnaire.

%% subject selection
[study_nm, condition, gender, subject_id, NS] = sub_id;

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
plasmaPath = [fullfile('P:','boulot','postdoc_CarmenSandi','results',...
    'plasma'),filesep];

%% load nutrition data
FFQ_struct = load([nutritionPath,'FFQ_metabolites_extracted_results.mat'],...
    'FFQ_results',...
    'energy_kcalPerWeek','metabolite_names',...
    'subject_id');
energy_kcalPerWeek = FFQ_struct.energy_kcalPerWeek;
foodNutrient_names = FFQ_struct.metabolite_names;
nNutrients = length(foodNutrient_names);
FFQ_subject_id = FFQ_struct.subject_id;

%% load plasma levels
blood_AA_table = readtable([plasmaPath, 'Blood_AA.xlsx'],...
    'Sheet','filtered_data_clean');
%% initialize variables of interest
[plasma.L_Gly, plasma.L_Glu, plasma.L_Gln, plasma.Tau,...
    nutri.Gly, nutri.Cys, nutri.Glu, nutri.Tau,...
    nutri.calories,...
    nutri.Gly_div_totalCal,...
    nutri.Cys_div_totalCal,...
    nutri.Glu_div_totalCal,...
    nutri.Tau_div_totalCal] = deal(NaN(1,NS));
nutri_vars = fieldnames(nutri);
%% extract data for each
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    %% load nutrition scores
    sub_nutri_idx = find(strcmp(FFQ_subject_id, sub_nm));
    % extract nutrition intake values
    if ~isempty(sub_nutri_idx) && size(sub_nutri_idx,1) == 1
        nutri.Gly(iS) = FFQ_struct.FFQ_results.Gly_mgPerWeek(sub_nutri_idx);
        nutri.Cys(iS) = FFQ_struct.FFQ_results.Cys_mgPerWeek(sub_nutri_idx);
        nutri.Glu(iS) = FFQ_struct.FFQ_results.Glu_mgPerWeek(sub_nutri_idx);
        nutri.calories(iS) = FFQ_struct.FFQ_results.calories_kcalPerWeek(sub_nutri_idx);
    end
    
    % divide by calories
    if ~isnan(nutri.calories(iS))
        nutri.Gly_div_totalCal(iS) = nutri.Gly(iS)/nutri.calories(iS);
        nutri.Cys_div_totalCal(iS) = nutri.Cys(iS)/nutri.calories(iS);
        nutri.Glu_div_totalCal(iS) = nutri.Gly(iS)/nutri.calories(iS);
    end
    
    %% load plasma scores
    sub_plasma_idx = find(strcmp(blood_AA_table.CID,['CID',sub_nm]));
    plasma.L_Gly(iS) = blood_AA_table.L_Glycine(sub_plasma_idx);
    plasma.L_Glu(iS) = blood_AA_table.L_GlutamicA(sub_plasma_idx);
    plasma.L_Gln(iS) = blood_AA_table.L_Glutamine(sub_plasma_idx);
    plasma.Tau(iS)   = blood_AA_table.Taurine(sub_plasma_idx);
end % subject loop

%% correlate both

% uncorrected nutritional metabolites
[r_corr.Gly, betas.Gly, pval.Gly,...
    ~,...
    nutri_Gly_sorted, plasma_Gly_fit_nutriGlySorted] = glm_package(nutri.Gly', plasma.L_Gly', 'normal', 'on');
[r_corr.Glu, betas.Glu, pval.Glu,...
    ~,...
    nutri_Glu_sorted, plasma_Glu_fit_nutriGluSorted] = glm_package(nutri.Glu', plasma.L_Glu', 'normal', 'on');
[r_corr.Gln_f_nutriGlu, betas.Gln_f_nutriGlu, pval.Gln_f_nutriGlu,...
    ~,...
    nutri_Glu_sorted, plasma_Gln_fit_nutriGluSorted] = glm_package(nutri.Glu', plasma.L_Gln', 'normal', 'on');
% [r_corr.Tau, betas.Tau, pval.Tau,...
%     ~,...
%     nutri_Tau_sorted, plasma_Tau_fit_nutriTauSorted] = glm_package(nutri.Tau', plasma.Tau', 'normal', 'on');

% nutritional metabolites corrected by total caloric intake
[r_corr.Gly_div_totalCal, betas.Gly_div_totalCal, pval.Gly_div_totalCal,...
    ~,...
    nutri_Gly_div_totalCal_sorted, plasma_Gly_fit_nutriGly_div_totalCalSorted] = glm_package(nutri.Gly_div_totalCal', plasma.L_Gly', 'normal', 'on');
[r_corr.Glu_div_totalCal, betas.Glu_div_totalCal, pval.Glu_div_totalCal,...
    ~,...
    nutri_Glu_div_totalCal_sorted, plasma_Glu_fit_nutriGlu_div_totalCalSorted] = glm_package(nutri.Glu_div_totalCal', plasma.L_Glu', 'normal', 'on');
[r_corr.Gln_f_nutriGlu_div_totalCal, betas.Gln_f_nutriGlu_div_totalCal, pval.Gln_f_nutriGlu_div_totalCal,...
    ~,...
    nutri_Glu_div_totalCal_sorted, plasma_Gln_fit_nutriGlu_div_totalCalSorted] = glm_package(nutri.Glu_div_totalCal', plasma.L_Gln', 'normal', 'on');
% [r_corr.Tau_div_totalCal, betas.Tau_div_totalCal, pval.Tau_div_totalCal,...
%     ~,...
%     nutri_Tau_div_totalCal_sorted, plasma_Tau_fit_nutriTau_div_totalCalSorted] = glm_package(nutri.Tau_div_totalCal', plasma.Tau', 'normal', 'on');

%% display figures
[pSize, lWidth, col, mSize] = general_fig_prm;

%% Gly
fig;
Gly_scat_hdl = scatter(nutri.Gly./1000, plasma.L_Gly);
Gly_fit_hdl = plot(nutri_Gly_sorted./1000, plasma_Gly_fit_nutriGlySorted);
scat_hdl_upgrade(Gly_scat_hdl);
fit_hdl_upgrade(Gly_fit_hdl);
xlabel('Nutrition Glycine (g/week)');
ylabel('Plasma Glycine (μM)');
place_r_and_pval(r_corr.Gly, pval.Gly(2));
legend_size(pSize);

%% Glu
fig;
Glu_scat_hdl = scatter(nutri.Glu./1000, plasma.L_Glu);
Glu_fit_hdl = plot(nutri_Glu_sorted./1000, plasma_Glu_fit_nutriGluSorted);
scat_hdl_upgrade(Glu_scat_hdl);
fit_hdl_upgrade(Glu_fit_hdl);
xlabel('Nutrition Glutamate (g/week)');
ylabel('Plasma Glutamate (μM)');
place_r_and_pval(r_corr.Glu, pval.Glu(2));
legend_size(pSize);

%% plasmaGln = f(nutri Glu)
fig;
Gln_scat_hdl = scatter(nutri.Glu./1000, plasma.L_Gln);
Gln_fit_hdl = plot(nutri_Glu_sorted./1000, plasma_Gln_fit_nutriGluSorted);
scat_hdl_upgrade(Gln_scat_hdl);
fit_hdl_upgrade(Gln_fit_hdl);
xlabel('Nutrition Glutamate (g/week)');
ylabel('Plasma Glutamine (μM)');
place_r_and_pval(r_corr.Gln_f_nutriGlu, pval.Gln_f_nutriGlu(2));
legend_size(pSize);

% %% Taurine (nothing yet, as the BLS 3.02 doesn't contain any information
% regarding the nutritional intake of Taurine)

% fig;
% Tau_scat_hdl = scatter(nutri.Tau./1000, plasma.L_Tau);
% Tau_fit_hdl = plot(nutri_Tau_sorted./1000, plasma_Tau_fit_nutriTauSorted);
% scat_hdl_upgrade(Tau_scat_hdl);
% fit_hdl_upgrade(Tau_fit_hdl);
% xlabel('Nutrition Taurine (g/week)');
% ylabel('Plasma Taurine (μM)');
% place_r_and_pval(r_corr.Tau, pval.Tau(2));
% legend_size(pSize);

%% same but for nutrition corrected by total caloric intake
%% Gly
fig;
Gly_div_totalCal_scat_hdl = scatter(nutri.Gly_div_totalCal./1000, plasma.L_Gly);
Gly_div_totalCal_fit_hdl = plot(nutri_Gly_div_totalCal_sorted./1000, plasma_Gly_fit_nutriGly_div_totalCalSorted);
scat_hdl_upgrade(Gly_div_totalCal_scat_hdl);
fit_hdl_upgrade(Gly_div_totalCal_fit_hdl);
xlabel('Nutrition Glycine (g/week)');
ylabel('Plasma Glycine (μM)');
place_r_and_pval(r_corr.Gly_div_totalCal, pval.Gly_div_totalCal(2));
legend_size(pSize);

%% Glu
fig;
Glu_div_totalCal_scat_hdl = scatter(nutri.Glu_div_totalCal./1000, plasma.L_Glu);
Glu_div_totalCal_fit_hdl = plot(nutri_Glu_div_totalCal_sorted./1000, plasma_Glu_fit_nutriGlu_div_totalCalSorted);
scat_hdl_upgrade(Glu_div_totalCal_scat_hdl);
fit_hdl_upgrade(Glu_div_totalCal_fit_hdl);
xlabel('Nutrition Glutamate (g/week)');
ylabel('Plasma Glutamate (μM)');
place_r_and_pval(r_corr.Glu_div_totalCal, pval.Glu_div_totalCal(2));
legend_size(pSize);

%% plasmaGln = f(nutri Glu)
fig;
Gln_div_totalCal_scat_hdl = scatter(nutri.Glu_div_totalCal./1000, plasma.L_Gln);
Gln_div_totalCal_fit_hdl = plot(nutri_Glu_div_totalCal_sorted./1000, plasma_Gln_fit_nutriGlu_div_totalCalSorted);
scat_hdl_upgrade(Gln_div_totalCal_scat_hdl);
fit_hdl_upgrade(Gln_div_totalCal_fit_hdl);
xlabel('Nutrition Glutamate (g/week)');
ylabel('Plasma Glutamine (μM)');
place_r_and_pval(r_corr.Gln_f_nutriGlu_div_totalCal, pval.Gln_f_nutriGlu_div_totalCal(2));
legend_size(pSize);

% %% Taurine (nothing yet, as the BLS 3.02 doesn't contain any information
% regarding the nutritional intake of Taurine)

% fig;
% Tau_div_totalCal_scat_hdl = scatter(nutri.Tau./1000, plasma.L_Tau);
% Tau_div_totalCal_fit_hdl = plot(nutri_Tau_div_totalCal_sorted./1000, plasma_Tau_fit_nutriTau_div_totalCalSorted);
% scat_hdl_upgrade(Tau_div_totalCal_scat_hdl);
% fit_hdl_upgrade(Tau_div_totalCal_fit_hdl);
% xlabel('Nutrition Taurine (g/week)');
% ylabel('Plasma Taurine (μM)');
% place_r_and_pval(r_corr.Tau_div_totalCal, pval.Tau_div_totalCal(2));
% legend_size(pSize);

% end % function