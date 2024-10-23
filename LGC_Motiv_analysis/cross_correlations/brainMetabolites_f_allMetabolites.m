% function[r_corr, pval, signif, NS_goodS] = brainMetabolites_f_allMetabolites(outlierF)
% [r_corr, pval, signif, NS_goodS] = brainMetabolites_f_allMetabolites(outlierF)
% brainMetabolites_f_allMetabolites test the correlation between one
% selected brain metabolite and all peripheral metabolites
%
% INPUTS
% outlierF: filter outlier in the variables included in each correlation
% test (1) or not (0)
%
% OUTPUTS
% r_corr: 1*n_metabolites correlation vector with the correlation 
% coefficients corresponding to each test
%
% pval: 1*n_metabolites vector with the corresponding p.values for r_corr
%
% signif: structure with the significant tests (p < 0.05 uncorrected for 
% multiple comparisons)
%
% NS_goodS: number of good subjects for each correlation test

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

%% load brain metabolite of interest
[brain_mb_allSubs, MRS_ROI_nm, brain_mb_nm] = metabolite_extraction(study_nm, subject_id);

%% load all peripheral variables
% load whole-blood NADomics
[wholeBlood_mb] = load_blood_NAD(study_nm, subject_id);
wholeB_mb_names = fieldnames(wholeBlood_mb);
n_wholeB_mb = length(wholeB_mb_names);
% load plasma
[plasmaM, plasma_mb_names, n_plasma_mb] = load_plasma_metabolites(subject_id);

% pool all metabolic measures together
for iWholeB_Mb = 1:n_wholeB_mb
    wholeB_mb_nm = wholeB_mb_names{iWholeB_Mb};
    periph_metabolism.(['wholeB_',wholeB_mb_nm]) = wholeBlood_mb.(wholeB_mb_nm);
end
for iPlasma_Mb = 1:n_plasma_mb
    plasma_mb_nm = plasma_mb_names{iPlasma_Mb};
    periph_metabolism.(['plasma_',plasma_mb_nm]) = plasmaM.(plasma_mb_nm);
end
periph_mb_names = fieldnames(periph_metabolism);
n_periph_mb = length(periph_mb_names);

%% perform correlations
[r_corr_vec, pval_vec] = deal(NaN(1, n_periph_mb));
signif = struct;
for iPeriphMb = 1:n_periph_mb
    periph_mb_nm = periph_mb_names{iPeriphMb};
    corr_nm = [brain_mb_nm,'_f_',periph_mb_nm];
    periph_mb_var = periph_metabolism.(periph_mb_nm);
    
    % filter outliers
    switch outlierF
        case 0
            brain_mb_var_tmp = brain_mb_allSubs;
            periph_mb_var_tmp = periph_mb_var;
        case 1
            [~,~,brain_mb_var_tmp] = rmv_outliers_3sd(brain_mb_allSubs);
            [~,~,periph_mb_var_tmp] = rmv_outliers_3sd(periph_mb_var);
    end
    
    goodS_tmp = ~isnan(periph_mb_var_tmp.*brain_mb_var_tmp);
    NS_goodS.(corr_nm) = sum(goodS_tmp);
    % perform the correlation
    [r_corr.(corr_nm), pval.(corr_nm)] = corr(periph_mb_var(goodS_tmp)', brain_mb_allSubs(goodS_tmp)');
    r_corr_vec(iPeriphMb) = r_corr.(corr_nm);
    pval_vec(iPeriphMb) = pval.(corr_nm);
    
    % store significant results
    if pval.(corr_nm) >= 0.05 && pval.(corr_nm) < 0.1
        signif.almostSignif.(corr_nm).r_corr = r_corr.(corr_nm);
        signif.almostSignif.(corr_nm).pval = pval.(corr_nm);
    elseif pval.(corr_nm) < 0.05
        signif.p005.(corr_nm).r_corr = r_corr.(corr_nm);
        signif.p005.(corr_nm).pval = pval.(corr_nm);
        % store also case when p.value really hyper significant
        if pval.(corr_nm) < 0.001
            signif.p0001.(corr_nm).r_corr = r_corr.(corr_nm);
            signif.p0001.(corr_nm).pval = pval.(corr_nm);
        end % p.value
    end % p.value
end % peripheral metabolites loop

%% show general graph
% general figure parameters
pSize = 15;
apply_pval_threshold = false; % display everything even not significant results
pval_threshold = []; % no pvalue threshold
disp_signif_stars = true; % display stars upon the significant correlations
% correlation range
corr_range = [-1 1];

corr_plot(r_corr_vec, pval_vec,...
    corr_range, periph_mb_names, [MRS_ROI_nm,' - ', brain_mb_nm], [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

% end % function