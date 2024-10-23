% function[GLM_b_F_vs_mb, pval_b_F_vs_mb, signif, N_goodS] = fatigue_f_metabolism_and_other_ctrl_vars(outlierF)
% [corr_F_vs_mb, pval_F_vs_mb, signif, N_goodS] = fatigue_f_metabolism_and_other_ctrl_vars(outlierF)
% fatigue_f_metabolism will look at the correlation matrix between the
% different fatigue metrics from questionnaires + behavior and the
% metabolic measures performed in the plasma, whole-blood and brain.
%
% INPUTS
% outlierF: filter outlier in the variables included in each correlation
% test (1) or not (0)
%
% OUTPUTS
% corr_F_vs_mb: n_fatigue_variables*n_metabolism_variables correlation matrix
% with the correlation coefficients corresponding to each test
%
% pval_F_vs_mb: n_fatigue_variables*n_metabolism_variables matrix with the
% corresponding p.values for corr_F_vs_mb
%
% signif: structure with the significant tests (p < 0.05 uncorrected for 
% multiple comparisons)
%
% N_goodS: number of good subjects for each correlation test
%
% See also fatigue_measurements_crosscorrel.m & fatigue_nm_rename.m

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

%% load general variables (BMI, age, sex)
[BMI, age, sex] = deal(NaN(NS,1));
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
sub_gal = excelReadGeneralFile.CID;
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = strcmp(['CID',sub_nm],sub_gal);
    BMI(iS) = excelReadGeneralFile.BMI(sub_idx); % BMI
    age(iS) = excelReadGeneralFile.Age_yearsOld_(sub_idx);% age in years
    sex(iS) = strcmp(excelReadGeneralFile.Sexe_femaleF_maleM_(sub_idx),'F') +...
        -strcmp(excelReadGeneralFile.Sexe_femaleF_maleM_(sub_idx),'M'); % -1=male/+1=female
end

%% load fatigue metrics
[fatigue_measures] = fatigue_pool(study_nm, condition, subject_id, NS, genderFilter);
fatigue_var_names = fieldnames(fatigue_measures);
fatigue_var_names(strcmp(fatigue_var_names,'sub_selection')) = [];
n_F_vars = length(fatigue_var_names);
fatigue_vars_bis = cell(n_F_vars,1);
% rename fatigue variables for the figure
for iF = 1:n_F_vars
    [fatigue_vars_bis{iF}] = fatigue_nm_rename(fatigue_var_names{iF});
end % loop over fatigue variables

%% load whole-blood metabolites
wholeBlood_mb = load_blood_NAD(study_nm, subject_id);
wholeB_mb_names = fieldnames(wholeBlood_mb);
n_wholeB_mb = length(wholeB_mb_names);
for iWholeB_Mb = 1:n_wholeB_mb
    wholeB_mb_nm = wholeB_mb_names{iWholeB_Mb};
    metabolism.(['wholeB_',wholeB_mb_nm]) = wholeBlood_mb.(wholeB_mb_nm);
end

%% pool all metabolic measures together
% load plasma metabolites
[plasmaM, plasma_mb_names, n_plasma_mb] = load_plasma_metabolites(subject_id);
for iPlasma_Mb = 1:n_plasma_mb
    plasma_mb_nm = plasma_mb_names{iPlasma_Mb};
    metabolism.(['plasma_',plasma_mb_nm]) = plasmaM.(plasma_mb_nm);
end

%% load brain metabolites
[~, ~, brain_metabolites_bis] = metabolite_load(subject_id);
brain_mb_names = fieldnames(brain_metabolites_bis.dmPFC);
n_brain_mb = length(brain_mb_names);
brain_mb_names_bis = cell(n_brain_mb*2,1);
% add dmPFC/dACC
for iBrain_Mb = 1:n_brain_mb
    brain_mb_nm = brain_mb_names{iBrain_Mb};
    metabolism.(['dmPFC_',brain_mb_nm]) = brain_metabolites_bis.dmPFC.(brain_mb_nm);
end
% add aIns
for iBrain_Mb = 1:n_brain_mb
    brain_mb_nm = brain_mb_names{iBrain_Mb};
    metabolism.(['aIns_',brain_mb_nm]) = brain_metabolites_bis.aIns.(brain_mb_nm);
end

%% general correlation matrix
metabolite_names = fieldnames(metabolism);
n_mb_vars = length(metabolite_names);
N_goodS = NaN(n_F_vars, n_mb_vars); % number of good subjects for each correlation

% general matrix
[GLM_b_F_vs_sex_mb, GLM_b_F_vs_age_mb, GLM_b_F_vs_BMI_mb, GLM_b_F_vs_mb,...
    pval_b_F_vs_sex_mb, pval_b_F_vs_age_mb, pval_b_F_vs_BMI_mb, pval_b_F_vs_mb] = deal(NaN(n_F_vars, n_mb_vars));
% matrix for fatigue & whole-blood NADomics
[GLM_b_F_vs_sex_wbNADomics, GLM_b_F_vs_age_wbNADomics, GLM_b_F_vs_BMI_wbNADomics, GLM_b_F_vs_wbNADomics,...
    pval_b_F_vs_sex_wbNADomics, pval_b_F_vs_age_wbNADomics, pval_b_F_vs_BMI_wbNADomics, pval_b_F_vs_wbNADomics] = deal(NaN(n_F_vars, n_wholeB_mb));
% matrix for fatigue & plasma metabolites
[GLM_b_F_vs_sex_plasmaMb, GLM_b_F_vs_age_plasmaMb, GLM_b_F_vs_BMI_plasmaMb, GLM_b_F_vs_plasmaMb,...
    pval_b_F_vs_sex_plasmaMb, pval_b_F_vs_age_plasmaMb, pval_b_F_vs_BMI_plasmaMb, pval_b_F_vs_plasmaMb] = deal(NaN(n_F_vars, n_plasma_mb));
% matrix for fatigue & brain metabolites
[GLM_b_F_vs_sex_brainMb, GLM_b_F_vs_age_brainMb, GLM_b_F_vs_BMI_brainMb, GLM_b_F_vs_brainMb,...
    pval_b_F_vs_sex_brainMb, pval_b_F_vs_age_brainMb, pval_b_F_vs_BMI_brainMb, pval_b_F_vs_brainMb] = deal(NaN(n_F_vars, n_brain_mb*2));

for iF = 1:n_F_vars
    F_nm = fatigue_var_names{iF};
    
    %% loop across all metabolites
    for iMb = 1:n_mb_vars
        mb_nm = metabolite_names{iMb};
        % remove (or not) the outliers in each measure
        switch outlierF
            case 0
                fatigue_var_tmp = fatigue_measures.(F_nm);
                metabolism_var_tmp = metabolism.(mb_nm);
            case 1 % remove outliers based on mean +/- 3*SD
                [~, ~, fatigue_var_tmp] = rmv_outliers_3sd(fatigue_measures.(F_nm));
                [~, ~, metabolism_var_tmp] = rmv_outliers_3sd(metabolism.(mb_nm));
        end
        goodS_tmp = ~isnan(BMI.*age.*sex.*metabolism_var_tmp'.*fatigue_var_tmp');
        N_goodS(iF, iMb) = sum(goodS_tmp);
        
        % perform the GLM
        X_all_mb = [sex(goodS_tmp),...
            zscore(age(goodS_tmp)),...
            zscore(BMI(goodS_tmp)),...
            zscore(metabolism.(mb_nm)(goodS_tmp))'];
        [b_tmp, ~, stats_tmp] = glmfit(X_all_mb, zscore(fatigue_measures.(F_nm)(goodS_tmp))','normal');
        GLM_b_F_vs_sex_mb(iF, iMb) = b_tmp(2);
        GLM_b_F_vs_age_mb(iF, iMb) = b_tmp(3);
        GLM_b_F_vs_BMI_mb(iF, iMb) = b_tmp(4);
        GLM_b_F_vs_mb(iF, iMb) = b_tmp(5);
        pval_b_F_vs_sex_mb(iF, iMb) = stats_tmp.p(2);
        pval_b_F_vs_age_mb(iF, iMb) = stats_tmp.p(3);
        pval_b_F_vs_BMI_mb(iF, iMb) = stats_tmp.p(4);
        pval_b_F_vs_mb(iF, iMb) = stats_tmp.p(5);
        % store corresponding fit
        corr_nm = [F_nm,'_f_',mb_nm];
        [z_metabolism_fit.(corr_nm), idx_sort_tmp] = sort(zscore(metabolism.(mb_nm)(goodS_tmp)));
        fatigue_measures_fit.(corr_nm) = glmval(b_tmp, X_all_mb(idx_sort_tmp,:), 'identity', 'Constant','on');
        
        % store the significant correlations
        if pval_b_F_vs_mb(iF,iMb) < 0.05
            signif.([F_nm,'_f_',mb_nm]).p005.b = GLM_b_F_vs_mb(iF, iMb);
            signif.([F_nm,'_f_',mb_nm]).p005.p = pval_b_F_vs_mb(iF, iMb);
        end
        if pval_b_F_vs_mb(iF,iMb) < 0.01
            signif.([F_nm,'_f_',mb_nm]).p001.b = GLM_b_F_vs_mb(iF, iMb);
            signif.([F_nm,'_f_',mb_nm]).p001.p = pval_b_F_vs_mb(iF, iMb);
        end
        if pval_b_F_vs_mb(iF,iMb) < 0.005
            signif.([F_nm,'_f_',mb_nm]).p0005.b = GLM_b_F_vs_mb(iF, iMb);
            signif.([F_nm,'_f_',mb_nm]).p0005.p = pval_b_F_vs_mb(iF, iMb);
        end
        if pval_b_F_vs_mb(iF,iMb) < 0.001
            signif.([F_nm,'_f_',mb_nm]).p0001.b = GLM_b_F_vs_mb(iF, iMb);
            signif.([F_nm,'_f_',mb_nm]).p0001.p = pval_b_F_vs_mb(iF, iMb);
        end
    end % metabolite loop
    
    %% loop across whole-blood metabolites
    for iWholeB_Mb = 1:n_wholeB_mb
        wholeB_mb_nm = wholeB_mb_names{iWholeB_Mb};
        % remove (or not) the outliers in each measure
        switch outlierF
            case 0
                fatigue_var_tmp = fatigue_measures.(F_nm);
                wholeBloodMb_var_tmp = wholeBlood_mb.(wholeB_mb_nm);
            case 1 % remove outliers based on mean +/- 3*SD
                [~, ~, fatigue_var_tmp] = rmv_outliers_3sd(fatigue_measures.(F_nm));
                [~, ~, wholeBloodMb_var_tmp] = rmv_outliers_3sd(wholeBlood_mb.(wholeB_mb_nm));
        end
        wB_goodS_tmp = ~isnan(BMI.*age.*sex.*wholeBloodMb_var_tmp'.*fatigue_var_tmp');
        X_wholeB = [sex(wB_goodS_tmp),...
            zscore(age(wB_goodS_tmp)),...
            zscore(BMI(wB_goodS_tmp)),...
            zscore(wholeBlood_mb.(wholeB_mb_nm)(wB_goodS_tmp))'];
        [b_wB_tmp, ~, stats_wB] = glmfit(X_wholeB, zscore(fatigue_measures.(F_nm)(wB_goodS_tmp))','normal');
        GLM_b_F_vs_sex_wbNADomics(iF, iWholeB_Mb)  = b_wB_tmp(2);
        GLM_b_F_vs_age_wbNADomics(iF, iWholeB_Mb)  = b_wB_tmp(3);
        GLM_b_F_vs_BMI_wbNADomics(iF, iWholeB_Mb)  = b_wB_tmp(4);
        GLM_b_F_vs_wbNADomics(iF, iWholeB_Mb)      = b_wB_tmp(5);
        pval_b_F_vs_sex_wbNADomics(iF, iWholeB_Mb) = stats_wB.p(2);
        pval_b_F_vs_age_wbNADomics(iF, iWholeB_Mb) = stats_wB.p(3);
        pval_b_F_vs_BMI_wbNADomics(iF, iWholeB_Mb) = stats_wB.p(4);
        pval_b_F_vs_wbNADomics(iF, iWholeB_Mb)     = stats_wB.p(5);
    end % whole-blood NADomics
        
    %% loop across plasma metabolites
    for iPlasma_Mb = 1:n_plasma_mb
        plasma_mb_nm = plasma_mb_names{iPlasma_Mb};
        % remove (or not) the outliers in each measure
        switch outlierF
            case 0
                fatigue_var_tmp = fatigue_measures.(F_nm);
                plasmaMb_var_tmp = plasmaM.(plasma_mb_nm);
            case 1 % remove outliers based on mean +/- 3*SD
                [~, ~, fatigue_var_tmp] = rmv_outliers_3sd(fatigue_measures.(F_nm));
                [~, ~, plasmaMb_var_tmp] = rmv_outliers_3sd(plasmaM.(plasma_mb_nm));
        end
        plasma_goodS_tmp = ~isnan(BMI.*age.*sex.*fatigue_var_tmp'.*plasmaMb_var_tmp');
        X_fatigue_plasmaMb = [sex(plasma_goodS_tmp),...
            zscore(age(plasma_goodS_tmp)),...
            zscore(BMI(plasma_goodS_tmp)),...
            zscore(plasmaM.(plasma_mb_nm)(plasma_goodS_tmp)')];
        [b_pMb_tmp, ~, stats_pMb] = glmfit(X_fatigue_plasmaMb, zscore(fatigue_measures.(F_nm)(plasma_goodS_tmp))','normal');
        GLM_b_F_vs_sex_plasmaMb(iF, iPlasma_Mb)  = b_pMb_tmp(2);
        GLM_b_F_vs_age_plasmaMb(iF, iPlasma_Mb)  = b_pMb_tmp(3);
        GLM_b_F_vs_BMI_plasmaMb(iF, iPlasma_Mb)  = b_pMb_tmp(4);
        GLM_b_F_vs_plasmaMb(iF, iPlasma_Mb)      = b_pMb_tmp(5);
        pval_b_F_vs_sex_plasmaMb(iF, iPlasma_Mb) = stats_pMb.p(2);
        pval_b_F_vs_age_plasmaMb(iF, iPlasma_Mb) = stats_pMb.p(3);
        pval_b_F_vs_BMI_plasmaMb(iF, iPlasma_Mb) = stats_pMb.p(4);
        pval_b_F_vs_plasmaMb(iF, iPlasma_Mb)     = stats_pMb.p(5);
    end % plasma metabolites
    
    %% loop across brain metabolites
    for iBrain_Mb = 1:n_brain_mb
        brain_mb_nm = brain_mb_names{iBrain_Mb};
        % dmPFC/dACC
        j_dmPFC_mb = 1 + 2*(iBrain_Mb - 1);
        brain_mb_names_bis{j_dmPFC_mb} = ['dmPFC ',brain_mb_nm];
        % remove (or not) the outliers in each measure
        switch outlierF
            case 0
                fatigue_var_tmp = fatigue_measures.(F_nm);
                dmPFC_mb_var_tmp = brain_metabolites_bis.dmPFC.(brain_mb_nm);
                aIns_mb_var_tmp = brain_metabolites_bis.aIns.(brain_mb_nm);
            case 1 % remove outliers based on mean +/- 3*SD
                [~, ~, fatigue_var_tmp] = rmv_outliers_3sd(fatigue_measures.(F_nm));
                [~, ~, dmPFC_mb_var_tmp] = rmv_outliers_3sd(brain_metabolites_bis.dmPFC.(brain_mb_nm));
                [~, ~, aIns_mb_var_tmp] = rmv_outliers_3sd(brain_metabolites_bis.aIns.(brain_mb_nm));
        end
        dmPFC_goodS_tmp = ~isnan(BMI.*age.*sex.*fatigue_var_tmp'.*dmPFC_mb_var_tmp');
        X_fatigue_dmPFCMb = [sex(dmPFC_goodS_tmp),...
            zscore(age(dmPFC_goodS_tmp)),...
            zscore(BMI(dmPFC_goodS_tmp)),...
            zscore(brain_metabolites_bis.dmPFC.(brain_mb_nm)(dmPFC_goodS_tmp)')];
        [b_dmPFCmb_tmp, ~, stats_dmPFC_mb] = glmfit(X_fatigue_dmPFCMb, zscore(fatigue_measures.(F_nm)(dmPFC_goodS_tmp))','normal');
        GLM_b_F_vs_sex_brainMb(iF, j_dmPFC_mb)  = b_dmPFCmb_tmp(2);
        GLM_b_F_vs_age_brainMb(iF, j_dmPFC_mb)  = b_dmPFCmb_tmp(3);
        GLM_b_F_vs_BMI_brainMb(iF, j_dmPFC_mb)  = b_dmPFCmb_tmp(4);
        GLM_b_F_vs_brainMb(iF, j_dmPFC_mb)      = b_dmPFCmb_tmp(5);
        pval_b_F_vs_sex_brainMb(iF, j_dmPFC_mb) = stats_dmPFC_mb.p(2);
        pval_b_F_vs_age_brainMb(iF, j_dmPFC_mb) = stats_dmPFC_mb.p(3);
        pval_b_F_vs_BMI_brainMb(iF, j_dmPFC_mb) = stats_dmPFC_mb.p(4);
        pval_b_F_vs_brainMb(iF, j_dmPFC_mb)     = stats_dmPFC_mb.p(5);
        
        % anterior insula
        j_aIns_mb = 2 + 2*(iBrain_Mb - 1);
        brain_mb_names_bis{j_aIns_mb} = ['aIns ',brain_mb_nm];
        aIns_goodS_tmp = ~isnan(BMI.*age.*sex.*fatigue_var_tmp'.*aIns_mb_var_tmp');
        X_fatigue_aInsMb = [sex(aIns_goodS_tmp),...
            zscore(age(aIns_goodS_tmp)),...
            zscore(BMI(aIns_goodS_tmp)),...
            zscore(brain_metabolites_bis.aIns.(brain_mb_nm)(aIns_goodS_tmp)')];
        [b_aInsMb_tmp, ~, stats_aIns_mb] = glmfit(X_fatigue_aInsMb,...
            zscore(fatigue_measures.(F_nm)(aIns_goodS_tmp))','normal');
        GLM_b_F_vs_sex_brainMb(iF, j_aIns_mb)   = b_aInsMb_tmp(2);
        GLM_b_F_vs_age_brainMb(iF, j_aIns_mb)   = b_aInsMb_tmp(3);
        GLM_b_F_vs_BMI_brainMb(iF, j_aIns_mb)   = b_aInsMb_tmp(4);
        GLM_b_F_vs_brainMb(iF, j_aIns_mb)       = b_aInsMb_tmp(5);
        pval_b_F_vs_sex_brainMb(iF, j_aIns_mb)  = stats_aIns_mb.p(2);
        pval_b_F_vs_age_brainMb(iF, j_aIns_mb)  = stats_aIns_mb.p(3);
        pval_b_F_vs_BMI_brainMb(iF, j_aIns_mb)  = stats_aIns_mb.p(4);
        pval_b_F_vs_brainMb(iF, j_aIns_mb)      = stats_aIns_mb.p(5);
    end % brain metabolites
end % fatigue loop

%% display resulting correlation matrix
% general figure parameters
pSize = 15;

% correlation range
corr_range = [-1 1];

apply_pval_threshold = false; % display everything even not significant results
pval_threshold = []; % no pvalue threshold
disp_signif_stars = true; % display stars upon the significant correlations

%% general correlation matrix
corr_plot(GLM_b_F_vs_mb, pval_b_F_vs_mb,...
    corr_range, [], fatigue_vars_bis, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

%% correlation matrix whole-blood NADomics
corr_plot(GLM_b_F_vs_wbNADomics, pval_b_F_vs_wbNADomics,...
    corr_range, wholeB_mb_names, fatigue_vars_bis, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

%% correlation matrix plasma metabolites
corr_plot(GLM_b_F_vs_plasmaMb, pval_b_F_vs_plasmaMb,...
    corr_range, plasma_mb_names, fatigue_vars_bis, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

%% correlation matrix brain metabolites
corr_plot(GLM_b_F_vs_brainMb, pval_b_F_vs_brainMb,...
    corr_range, brain_mb_names_bis, fatigue_vars_bis, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

%% additional correlation plots
pSize = 40;
detail_plots = questdlg('Do you want to zoom in some of these correlations?','correlation plots?','No','Yes','Yes');
switch detail_plots
    case 'Yes'
        n_plots_str = inputdlg('How many plots do you want to see?','N.plots');
        n_plots = str2double(n_plots_str);
        for iPlot = 1:n_plots
            % select variable
            jMb = listdlg('PromptString','Select metabolite',...
                'SelectionMode','single','ListString',metabolite_names);
            mb_var_nm = metabolite_names{jMb};
            mb_var = metabolism.(mb_var_nm);
            jF = listdlg('PromptString','Select fatigue variable',...
                'SelectionMode','single','ListString',fatigue_var_names);
            fatigue_var_nm = fatigue_var_names{jF};
            fatigue_var = fatigue_measures.(fatigue_var_nm);
            goodS_tmp = ~isnan(mb_var.*fatigue_var);
            
            % display corresponding figure + correlation coefficient and
            % p.value
            fig;
            scat_hdl = scatter(zscore(mb_var(goodS_tmp)), zscore(fatigue_var(goodS_tmp)));
            scat_hdl_upgrade(scat_hdl);
            corr_nm = [fatigue_var_nm,'_f_',mb_var_nm];
            fit_hdl = plot(z_metabolism_fit.(corr_nm), fatigue_measures_fit.(corr_nm));
            fit_hdl_upgrade(fit_hdl);
            xlabel(['z(',strrep(mb_var_nm,'_',' '),')']);
            ylabel(['z(',strrep(fatigue_var_nm,'_',' '),')']);
            place_r_and_pval(GLM_b_F_vs_mb(jF, jMb), pval_b_F_vs_mb(jF, jMb));
            legend_size(pSize);
        end % plot loop
    case 'No'
end % show correlation plots in more details

% end % function