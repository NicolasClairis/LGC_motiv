function[r_corr, pval, signif, N_goodS] = covid_f_wholeBlood_NADomics()
% [r_corr, pval, signif, N_goodS] = covid_f_wholeBlood_NADomics()
% covid_f_wholeBlood_NADomics will display the correlation between the
% number of covid infections, the distance to the last infection and each 
% metabolic measure extracted with whole-blood NAD-omics.
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

%% load whole-blood metabolites
[wholeBlood_mb, sub_List] = load_blood_NAD(study_nm, subject_id);
wholeB_mb_names = fieldnames(wholeBlood_mb);
n_wholeB_mb = length(wholeB_mb_names);

%% loop through metabolites
[N_goodS.nCovid, N_goodS.covidDist,...
    r_corr.nCovid_vs_wbNADomics, pval.nCovid_vs_wbNADomics,...
    r_corr.covidDist_vs_wbNADomics, pval.covidDist_vs_wbNADomics] = deal(NaN(1, n_wholeB_mb));
for iWB_mb = 1:n_wholeB_mb
    wholeB_mb_nm = wholeB_mb_names{iWB_mb};
    
    %% remove (or not) the outliers in each measure
    switch outlierF
        case 0
            n_covid_tmp = n_covid;
            covid_dist_tmp = covid_dist;
            wholeBlood_mb_var_tmp = wholeBlood_mb.(wholeB_mb_nm);
        case 1 % remove outliers based on mean +/- 3*SD
            [~, ~, n_covid_tmp] = rmv_outliers_3sd(n_covid);
            [~, ~, covid_dist_tmp] = rmv_outliers_3sd(covid_dist);
            [~, ~, wholeBlood_mb_var_tmp] = rmv_outliers_3sd(wholeBlood_mb.(wholeB_mb_nm));
    end
    
    %% perform correlations
    fig;
    
    % n.covid
    goodS_nCovid_tmp = ~isnan(n_covid_tmp.*wholeBlood_mb_var_tmp);
    N_goodS.nCovid(iWB_mb) = sum(goodS_nCovid_tmp);
    [r_corr.nCovid_vs_wbNADomics(iWB_mb), pval.nCovid_vs_wbNADomics(iWB_mb)] = corr(n_covid(goodS_nCovid_tmp)', wholeBlood_mb.(wholeB_mb_nm)(goodS_nCovid_tmp)');
    % store the significant correlations
    if pval.nCovid_vs_wbNADomics(iWB_mb) < 0.05
        signif.(['nCovid_f_',wholeB_mb_nm]).p005.r = r_corr.nCovid_vs_wbNADomics(iWB_mb);
        signif.(['nCovid_f_',wholeB_mb_nm]).p005.p = pval.nCovid_vs_wbNADomics(iWB_mb);
    end
    % figure
    subplot(1,2,1);
    nCovid_hdl = scatter(n_covid(goodS_nCovid_tmp)', wholeBlood_mb.(wholeB_mb_nm)(goodS_nCovid_tmp)');
    scat_hdl_upgrade(nCovid_hdl);
    xlabel('Number of covid infections');
    ylabel([wholeB_mb_nm,' (mM)']);
    place_r_and_pval(r_corr.nCovid_vs_wbNADomics(iWB_mb), pval.nCovid_vs_wbNADomics(iWB_mb));
    
    % covid distance
    goodS_covidDist_tmp = ~isnan(covid_dist_tmp.*wholeBlood_mb_var_tmp);
    N_goodS.covidDist(iWB_mb) = sum(goodS_covidDist_tmp);
    [r_corr.covidDist_vs_wbNADomics(iWB_mb), pval.covidDist_vs_wbNADomics(iWB_mb)] = corr(covid_dist(goodS_covidDist_tmp)', wholeBlood_mb.(wholeB_mb_nm)(goodS_covidDist_tmp)');
    % store the significant correlations
    if pval.covidDist_vs_wbNADomics(iWB_mb) < 0.05
        signif.(['covidDist_f_',wholeB_mb_nm]).p005.r = r_corr.covidDist_vs_wbNADomics(iWB_mb);
        signif.(['covidDist_f_',wholeB_mb_nm]).p005.p = pval.covidDist_vs_wbNADomics(iWB_mb);
    end
    % figure
    subplot(1,2,2);
    covidDist_hdl = scatter(covid_dist(goodS_nCovid_tmp)', wholeBlood_mb.(wholeB_mb_nm)(goodS_nCovid_tmp)');
    scat_hdl_upgrade(covidDist_hdl);
    xlabel('Covid distance (months)');
    ylabel([wholeB_mb_nm,' (mM)']);
    place_r_and_pval(r_corr.covidDist_vs_wbNADomics(iWB_mb), pval.covidDist_vs_wbNADomics(iWB_mb));
end % metabolite loop

end % function