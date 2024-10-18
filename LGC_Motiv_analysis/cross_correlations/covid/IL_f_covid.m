function[r_corr, pval, signif, N_goodS] = IL_f_covid()
% [r_corr, pval, signif, N_goodS] = IL_f_covid()
% IL_f_covid will display the correlation between the
% number of covid infections, the distance to the last infection and
% salivary interleukin measurements.
%
% INPUTS
% 
%
% OUTPUTS
% r_corr: structure with correlation coefficient for each test
%
% pval: structure with p.value for each correlation
%
% signif: significant correlations (based on p < 0.05)
%
% N_goodS: structure with number of good subjects for each test (will vary
% depending on if outlier filtering or not

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

%% load number of covid infections and date of infection
[questionnaires, categ_quests, n_categ, subject_id] = extract_questionnaires(study_nm, subject_id, NS);
n_covid = questionnaires.general.n_covid;
covid_dist = questionnaires.general.covid_months_since_last_infection; % distance with last covid in months

%% load interleukins
[IL_data] = load_IL(study_nm, subject_id);
IL_names = {'IL1b','IL6','IL18'};
n_IL = length(IL_names);

%% loop through interleukins
for iIL = 1:n_IL
    IL_nm = IL_names{iIL};
    
    %% remove (or not) the outliers in each measure
    switch outlierF
        case 0
            n_covid_tmp = n_covid;
            covid_dist_tmp = covid_dist;
            IL_var_tmp = IL_data.(IL_nm);
        case 1 % remove outliers based on mean +/- 3*SD
            [~, ~, n_covid_tmp] = rmv_outliers_3sd(n_covid);
            [~, ~, covid_dist_tmp] = rmv_outliers_3sd(covid_dist);
            [~, ~, IL_var_tmp] = rmv_outliers_3sd(IL_data.(IL_nm));
    end
    
    %% perform correlations
    fig;
    
    % n.covid
    n_covid_f_IL_nm = ['nCovid_f_',IL_nm];
    goodS_nCovid_tmp = ~isnan(n_covid_tmp.*IL_var_tmp);
    N_goodS.(n_covid_f_IL_nm) = sum(goodS_nCovid_tmp);
    [r_corr.(n_covid_f_IL_nm)(iIL), pval.(n_covid_f_IL_nm)(iIL)] = corr(n_covid(goodS_nCovid_tmp)', IL_data.(IL_nm)(goodS_nCovid_tmp));
    % store the significant correlations
    if pval.(n_covid_f_IL_nm)(iIL) < 0.05
        signif.(n_covid_f_IL_nm).p005.r = r_corr.(n_covid_f_IL_nm)(iIL);
        signif.(n_covid_f_IL_nm).p005.p = pval.(n_covid_f_IL_nm)(iIL);
    end
    % figure
    subplot(1,2,1);
    nCovid_hdl = scatter(n_covid(goodS_nCovid_tmp)', IL_data.(IL_nm)(goodS_nCovid_tmp));
    scat_hdl_upgrade(nCovid_hdl);
    xlabel('Number of covid infections');
    ylabel([IL_nm,' (pg/mL)']);
    place_r_and_pval(r_corr.(n_covid_f_IL_nm)(iIL), pval.(n_covid_f_IL_nm)(iIL));
    
    % covid distance
    covidDist_f_IL_nm = ['nCovid_f_',IL_nm];
    goodS_covidDist_tmp = ~isnan(covid_dist_tmp.*IL_var_tmp);
    N_goodS.(covidDist_f_IL_nm)(iIL) = sum(goodS_covidDist_tmp);
    [r_corr.(covidDist_f_IL_nm)(iIL), pval.(covidDist_f_IL_nm)(iIL)] = corr(covid_dist(goodS_covidDist_tmp)', IL_data.(IL_nm)(goodS_covidDist_tmp));
    % store the significant correlations
    if pval.(covidDist_f_IL_nm)(iIL) < 0.05
        signif.(covidDist_f_IL_nm).p005.r = r_corr.(covidDist_f_IL_nm)(iIL);
        signif.(covidDist_f_IL_nm).p005.p = pval.(covidDist_f_IL_nm)(iIL);
    end
    % figure
    subplot(1,2,2);
    covidDist_hdl = scatter(covid_dist(goodS_nCovid_tmp)', IL_data.(IL_nm)(goodS_nCovid_tmp));
    scat_hdl_upgrade(covidDist_hdl);
    xlabel('Covid distance (months)');
    ylabel([IL_nm,' (mM)']);
    place_r_and_pval(r_corr.(covidDist_f_IL_nm)(iIL), pval.(covidDist_f_IL_nm)(iIL));
end % interleukin loop

end % function