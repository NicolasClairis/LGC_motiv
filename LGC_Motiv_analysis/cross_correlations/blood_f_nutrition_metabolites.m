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
    case {[fullfile('C:','Users','Loco','Downloads'),filesep],...
            [fullfile('L:','human_data_private','raw_data_subject'),filesep]}
        gitPath = fullfile('C:','Users','Loco','Downloads');
    otherwise
        error('case not ready yet');
end
nutritionPath = [fullfile(gitPath,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'nutrition'),filesep];
plasmaPath = [fullfile('P:','boulot','postdoc_CarmenSandi','results',...
    'plasma'),filesep];

%% load nutrition data
Gly_table = readtable([nutritionPath, 'Glycine_scoring.xlsx'],...
    'Sheet','Sheet1');
Cys_table = readtable([nutritionPath, 'cystein_scoring.xlsx'],...
    'Sheet','Sheet1');
Glu_table = readtable([nutritionPath, 'Glutamate_scoring.xlsx'],...
    'Sheet','Sheet1');
calories_table = readtable([nutritionPath, 'calories_scoring.xlsx'],...
    'Sheet','Sheet1');

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
    sub_Gly_idx = find(strcmp(Gly_table.CID, sub_nm));
    sub_Cys_idx = find(strcmp(Cys_table.CID, sub_nm));
    sub_Glu_idx = find(strcmp(Glu_table.CID, sub_nm));
    sub_calories_idx = find(strcmp(calories_table.CID, sub_nm));
    % extract nutrition intake values
    if ~isempty(sub_Gly_idx) && size(sub_Gly_idx,1) == 1
        nutri.Gly(iS) = Gly_table.GlyParSemaine_mg_(sub_Gly_idx);
    end
    if ~isempty(sub_Cys_idx) && size(sub_Cys_idx,1) == 1
        nutri.Cys(iS) = Cys_table.CysteineParSemaine_mg_(sub_Cys_idx);
    end
    if ~isempty(sub_Glu_idx) && size(sub_Glu_idx,1) == 1
        nutri.Glu(iS) = Glu_table.GlutamateParSemaine_mg_(sub_Glu_idx);
    end
    if ~isempty(sub_calories_idx) && size(sub_calories_idx,1) == 1
        nutri.calories(iS) = calories_table.CaloriesParSemaine_kcal_(sub_calories_idx);
    end
    
    % divide by calories
    if ~isnan(nutri.calories(iS))
        nutri.Gly_div_totalCal(iS) = nutri.Gly(iS)/nutri.calories(iS);
        nutri.Cys_div_totalCal(iS) = nutri.Cys(iS)/nutri.calories(iS);
        nutri.Glu_div_totalCal(iS) = nutri.Gly(iS)/nutri.calories(iS);
    end
    
    %% load plasma scores
    sub_plasma_idx = find(blood_AA_table.MasterSampleID == str2double(sub_nm));
    plasma.L_Gly(iS) = blood_AA_table.L_Glycine(sub_plasma_idx);
    plasma.L_Glu(iS) = blood_AA_table.L_GlutamicA(sub_plasma_idx);
    plasma.L_Gln(iS) = blood_AA_table.L_Glutamine(sub_plasma_idx);
    plasma.Tau(iS)   = blood_AA_table.Taurine(sub_plasma_idx);
end % subject loop

%% correlate both
[r_corr.Gly, betas.Gly, pval.Gly,...
    ~,...
    nutri_Gly_sorted, plasma_Gly_fit_nutriGlySorted] = glm_package(nutri.Gly', plasma.L_Gly', 'normal', 'on');
[r_corr.Glu, betas.Glu, pval.Glu,...
    ~,...
    nutri_Glu_sorted, plasma_Glu_fit_nutriGluSorted] = glm_package(nutri.Glu', plasma.L_Glu', 'normal', 'on');
[r_corr.Gln_f_nutriGlu, betas.Gln_f_nutriGlu, pval.Gln_f_nutriGlu,...
    ~,...
    nutri_Glu_sorted, plasma_Gln_fit_nutriGluSorted] = glm_package(nutri.Glu', plasma.L_Gln', 'normal', 'on');
[r_corr.Tau, betas.Tau, pval.Tau,...
    ~,...
    nutri_Tau_sorted, plasma_Tau_fit_nutriTauSorted] = glm_package(nutri.Tau', plasma.Tau', 'normal', 'on');

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

%% Taurine
fig;
Tau_scat_hdl = scatter(nutri.Tau./1000, plasma.L_Tau);
Tau_fit_hdl = plot(nutri_Tau_sorted./1000, plasma_Tau_fit_nutriTauSorted);
scat_hdl_upgrade(Tau_scat_hdl);
fit_hdl_upgrade(Tau_fit_hdl);
xlabel('Nutrition Taurine (g/week)');
ylabel('Plasma Taurine (μM)');
place_r_and_pval(r_corr.Tau, pval.Tau(2));
legend_size(pSize);

% end % function