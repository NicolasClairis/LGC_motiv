function[corr_F_vs_mb, pval_F_vs_mb, signif, N_goodS] = fatigue_f_metabolism(outlierF)
% [corr_F_vs_mb, pval_F_vs_mb, signif, N_goodS] = fatigue_f_metabolism(outlierF)
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

%% load fatigue metrics
[fatigue_measures] = fatigue_pool(study_nm, condition, subject_id, NS, genderFilter);
fatigue_vars = fieldnames(fatigue_measures);
fatigue_vars(strcmp(fatigue_vars,'sub_selection')) = [];
n_F_vars = length(fatigue_vars);
fatigue_vars_bis = cell(n_F_vars,1);
% rename fatigue variables for the figure
for iF = 1:n_F_vars
    [fatigue_vars_bis{iF}] = fatigue_nm_rename(fatigue_vars{iF});
end % loop over fatigue variables

%% load whole-blood metabolites
[wholeBlood_mb, sub_List] = load_blood_NAD(study_nm, subject_id);
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
[corr_F_vs_mb, pval_F_vs_mb] = deal(NaN(n_F_vars, n_mb_vars));
% matrix for fatigue & whole-blood NADomics
[corr_F_vs_wbNADomics, pval_F_vs_wbNADomics] = deal(NaN(n_F_vars, n_wholeB_mb));
% matrix for fatigue & plasma metabolites
[corr_F_vs_plasmaMb, pval_F_vs_plasmaMb] = deal(NaN(n_F_vars, n_plasma_mb));
% matrix for fatigue & brain metabolites
[corr_F_vs_brainMb, pval_F_vs_brainMb] = deal(NaN(n_F_vars, n_brain_mb*2));
for iF = 1:n_F_vars
    F_nm = fatigue_vars{iF};
    
    % loop across all metabolites
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
        goodS_tmp = ~isnan(fatigue_var_tmp.*metabolism_var_tmp);
        N_goodS(iF, iMb) = sum(goodS_tmp);
        % perform the correlation
        [corr_F_vs_mb(iF, iMb), pval_F_vs_mb(iF, iMb)] = corr(fatigue_measures.(F_nm)(goodS_tmp)', metabolism.(mb_nm)(goodS_tmp)');
        % store the significant correlations
        if pval_F_vs_mb(iF,iMb) < 0.05
            signif.([F_nm,'_f_',mb_nm]).p005.r = corr_F_vs_mb(iF, iMb);
            signif.([F_nm,'_f_',mb_nm]).p005.p = pval_F_vs_mb(iF, iMb);
        end
        if pval_F_vs_mb(iF,iMb) < 0.01
            signif.([F_nm,'_f_',mb_nm]).p001.r = corr_F_vs_mb(iF, iMb);
            signif.([F_nm,'_f_',mb_nm]).p001.p = pval_F_vs_mb(iF, iMb);
        end
        if pval_F_vs_mb(iF,iMb) < 0.005
            signif.([F_nm,'_f_',mb_nm]).p0005.r = corr_F_vs_mb(iF, iMb);
            signif.([F_nm,'_f_',mb_nm]).p0005.p = pval_F_vs_mb(iF, iMb);
        end
        if pval_F_vs_mb(iF,iMb) < 0.001
            signif.([F_nm,'_f_',mb_nm]).p0001.r = corr_F_vs_mb(iF, iMb);
            signif.([F_nm,'_f_',mb_nm]).p0001.p = pval_F_vs_mb(iF, iMb);
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
        goodS_tmp = ~isnan(fatigue_var_tmp.*wholeBloodMb_var_tmp);
        [corr_F_vs_wbNADomics(iF, iWholeB_Mb), pval_F_vs_wbNADomics(iF, iWholeB_Mb)] = corr(fatigue_measures.(F_nm)(goodS_tmp)', wholeBlood_mb.(wholeB_mb_nm)(goodS_tmp)');
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
        goodS_tmp = ~isnan(fatigue_var_tmp.*plasmaMb_var_tmp);
        [corr_F_vs_plasmaMb(iF, iPlasma_Mb), pval_F_vs_plasmaMb(iF, iPlasma_Mb)] = corr(fatigue_measures.(F_nm)(goodS_tmp)', plasmaM.(plasma_mb_nm)(goodS_tmp)');
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
        dmPFC_goodS_tmp = ~isnan(fatigue_var_tmp.*dmPFC_mb_var_tmp);
        [corr_F_vs_brainMb(iF, j_dmPFC_mb), pval_F_vs_brainMb(iF, j_dmPFC_mb)] = corr(fatigue_measures.(F_nm)(dmPFC_goodS_tmp)', brain_metabolites_bis.dmPFC.(brain_mb_nm)(dmPFC_goodS_tmp)');
        % anterior insula
        j_aIns_mb = 2 + 2*(iBrain_Mb - 1);
        brain_mb_names_bis{j_aIns_mb} = ['aIns ',brain_mb_nm];
        aIns_goodS_tmp = ~isnan(fatigue_var_tmp.*aIns_mb_var_tmp);
        [corr_F_vs_brainMb(iF, j_aIns_mb), pval_F_vs_brainMb(iF, j_aIns_mb)] = corr(fatigue_measures.(F_nm)(aIns_goodS_tmp)', brain_metabolites_bis.aIns.(brain_mb_nm)(aIns_goodS_tmp)');
    end % brain metabolites
end % fatigue loop

%% display resulting correlation matrix
% general figure parameters
[~, ~, col] = general_fig_prm;
pSize = 15;
color_range_choices = redblue(45);

% correlation range
corr_range = [-1 1];

apply_pval_threshold = false; % display everything even not significant results
pval_threshold = []; % no pvalue threshold
disp_signif_stars = true; % display stars upon the significant correlations

%% general correlation matrix
corr_plot(corr_F_vs_mb, pval_F_vs_mb,...
    corr_range, [], fatigue_vars_bis, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

%% correlation matrix whole-blood NADomics
corr_plot(corr_F_vs_wbNADomics, pval_F_vs_wbNADomics,...
    corr_range, wholeB_mb_names, fatigue_vars_bis, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

%% correlation matrix plasma metabolites
corr_plot(corr_F_vs_plasmaMb, pval_F_vs_plasmaMb,...
    corr_range, plasma_mb_names, fatigue_vars_bis, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

%% correlation matrix brain metabolites
corr_plot(corr_F_vs_brainMb, pval_F_vs_brainMb,...
    corr_range, brain_mb_names_bis, fatigue_vars_bis, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

end % function