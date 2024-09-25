function[corr_F, pval_F, signif, N_goodS] = fatigue_measurements_crosscorrel(outlierF)
% [corr_F, pval_F, signif, N_goodS] = fatigue_measurements_crosscorrel(outlierF)
% fatigue_measurements_crosscorrel will look at the cross-correlation 
% between the different fatigue metrics from questionnaires + behavior.
%
% INPUTS
% outlierF: filter outlier in the variables included in each correlation
% test (1) or not (0)
%
% OUTPUTS
% corr_F: n_fatigue_variables*n_fatigue_variables cross-correlation matrix
% with the correlation coefficients corresponding to each test
%
% pval_F: n_fatigue_variables*n_fatigue_variables matrix with the
% corresponding p.values for corr_F
%
% signif: structure with the significant tests (p < 0.05 uncorrected for 
% multiple comparisons)
%
% N_goodS: number of good subjects for each correlation test

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
% improve labels
for iF1 = 1:n_F_vars
    F_nm1 = fatigue_vars{iF1};
    switch F_nm1
        case 'MPSTEFS_physical'
            fatigue_vars_bis{iF1} = 'MPSTEFS_p';
        case 'MPSTEFS_mental'
            fatigue_vars_bis{iF1} = 'MPSTEFS_m';
        case 'n_covid'
            fatigue_vars_bis{iF1} = 'n.covid';
        case 'F4_F1'
            fatigue_vars_bis{iF1} = 'F4-F1';
        case 'F4_F3'
            fatigue_vars_bis{iF1} = 'F4-F3';
        case 'MVC_R1a'
            fatigue_vars_bis{iF1} = 'MVC r1a';
        case 'MVC_R1b'
            fatigue_vars_bis{iF1} = 'MVC r1b';
        case 'MVC_R2a'
            fatigue_vars_bis{iF1} = 'MVC r2a';
        case 'MVC_R2b'
            fatigue_vars_bis{iF1} = 'MVC r2b';
        case 'MVC_dR1'
            fatigue_vars_bis{iF1} = 'MVC Δr1';
        case 'MVC_dR2'
            fatigue_vars_bis{iF1} = 'MVC Δr2';
        case 'MVC_dTask'
            fatigue_vars_bis{iF1} = 'MVC ΔTask';
        case 'perf_decrease_slope_E0'
            fatigue_vars_bis{iF1} = 'perf dE0';
        case 'perf_decrease_slope_E1'
            fatigue_vars_bis{iF1} = 'perf dE1';
        case 'perf_decrease_slope_E2'
            fatigue_vars_bis{iF1} = 'perf dE2';
        case 'perf_decrease_slope_E3'
            fatigue_vars_bis{iF1} = 'perf dE3';
        case 'choices_kFp'
            fatigue_vars_bis{iF1} = 'kFp';
        otherwise
            fatigue_vars_bis{iF1} = fatigue_vars{iF1};
    end % fatigue variable
end % loop over fatigue variables

%% cross-correlation matrix
N_goodS = NaN(n_F_vars, n_F_vars); % number of good subjects for each correlation
% general matrix
[corr_F, pval_F] = deal(NaN(n_F_vars, n_F_vars));
for iF1 = 1:n_F_vars
    F_nm1 = fatigue_vars{iF1};
    
    % loop across all metabolites
    for iF2 = 1:n_F_vars
        F_nm2 = fatigue_vars{iF2};
        % remove (or not) the outliers in each measure
        switch outlierF
            case 0
                fatigue_var1_tmp = fatigue_measures.(F_nm1);
                fatigue_var2_tmp = fatigue_measures.(F_nm2);
            case 1 % remove outliers based on mean +/- 3*SD
                [~, ~, fatigue_var1_tmp] = rmv_outliers_3sd(fatigue_measures.(F_nm1));
                [~, ~, fatigue_var2_tmp] = rmv_outliers_3sd(fatigue_measures.(F_nm2));
        end
        goodS_tmp = ~isnan(fatigue_var1_tmp.*fatigue_var2_tmp);
        N_goodS(iF1, iF2) = sum(goodS_tmp);
        % perform the correlation
        [corr_F(iF1, iF2), pval_F(iF1, iF2)] = corr(fatigue_measures.(F_nm1)(goodS_tmp)', fatigue_measures.(F_nm2)(goodS_tmp)');
        % store the significant correlations
        if pval_F(iF1, iF2) < 0.05
            signif.([F_nm1,'_f_',F_nm2]).p005.r = corr_F(iF1, iF2);
            signif.([F_nm1,'_f_',F_nm2]).p005.p = pval_F(iF1, iF2);
        end
        if pval_F(iF1, iF2) < 0.01
            signif.([F_nm1,'_f_',F_nm2]).p001.r = corr_F(iF1, iF2);
            signif.([F_nm1,'_f_',F_nm2]).p001.p = pval_F(iF1, iF2);
        end
        if pval_F(iF1, iF2) < 0.005
            signif.([F_nm1,'_f_',F_nm2]).p0005.r = corr_F(iF1, iF2);
            signif.([F_nm1,'_f_',F_nm2]).p0005.p = pval_F(iF1, iF2);
        end
        if pval_F(iF1, iF2) < 0.001
            signif.([F_nm1,'_f_',F_nm2]).p0001.r = corr_F(iF1, iF2);
            signif.([F_nm1,'_f_',F_nm2]).p0001.p = pval_F(iF1, iF2);
        end
    end % fatigue loop
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
corr_plot(corr_F, pval_F,...
    corr_range, fatigue_vars_bis, fatigue_vars_bis, [], [],...
    apply_pval_threshold, pval_threshold, disp_signif_stars);
legend_size(pSize);

end % function