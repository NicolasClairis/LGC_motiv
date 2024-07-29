% script to prepare correlation comparison in R script
% compare_correlations.R in order to identify specificity of the
% correlations identified in the paper.

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% initialize variables of interest (to force correct dimensions)
[dmPFC_Lac.allSubs,...
    aIns_Lac.allSubs,...
    plasma_Lac.allSubs,...
    HE.allSubs,...
    HPE.allSubs,...
    HME.allSubs,...
    kEp.allSubs,...
    kEm.allSubs,...
    kR.allSubs,...
    kP.allSubs,...
    kFp.allSubs,...
    kLm.allSubs,...
    kBias.allSubs] = deal(NaN(NS,1));

%% load dmPFC/dACC, aIns and plasma lactate
% load brain lactate
[metabolites, CRLB, metabolites_bis] = metabolite_load(subject_id);
dmPFC_Lac.allSubs(:) = metabolites.dmPFC.Lac;
aIns_Lac.allSubs(:) = metabolites.aIns.Lac;
% load plasma lactate
[plasmaStruct] = load_plasma_Lac(subject_id);
plasma_Lac.allSubs(:) = plasmaStruct.Lac;

%% load % high-effort choices for each task (HPE and HME)
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, 0);
HE.allSubs(:) = choice_hE.EpEm;
HPE.allSubs(:) = choice_hE.Ep;
HME.allSubs(:) = choice_hE.Em;

%% load behavioral parameters
[prm] = prm_extraction(study_nm, subject_id, 'bayesian', '5');
kEp.allSubs(:) = prm.kEp;
kEm.allSubs(:) = prm.kEm;
kR.allSubs(:) = prm.kR;
kP.allSubs(:) = prm.kP;
kFp.allSubs(:) = prm.kFp;
kLm.allSubs(:) = prm.kLm;
kBias.allSubs(:) = prm.kBias;

%% remove outliers within each variable and extract vector for each to know which subjects to remove
[~,~,dmPFC_Lac.noOutliers, idx_goodSubs.dmpfc_Lac] = rmv_outliers_3sd(dmPFC_Lac.allSubs');
[~,~,aIns_Lac.noOutliers, idx_goodSubs.ains_Lac] = rmv_outliers_3sd(aIns_Lac.allSubs');
[~,~,plasma_Lac.noOutliers, idx_goodSubs.plasma_Lac] = rmv_outliers_3sd(plasma_Lac.allSubs');
[~,~,HE.noOutliers, idx_goodSubs.HE] = rmv_outliers_3sd(HE.allSubs');
[~,~,HPE.noOutliers, idx_goodSubs.HPE] = rmv_outliers_3sd(HPE.allSubs');
[~,~,HME.noOutliers, idx_goodSubs.HME] = rmv_outliers_3sd(HME.allSubs');
[~,~,kEp.noOutliers, idx_goodSubs.kEp] = rmv_outliers_3sd(kEp.allSubs');
[~,~,kEm.noOutliers, idx_goodSubs.kEm] = rmv_outliers_3sd(kEm.allSubs');
[~,~,kR.noOutliers, idx_goodSubs.kR] = rmv_outliers_3sd(kR.allSubs');
[~,~,kP.noOutliers, idx_goodSubs.kP] = rmv_outliers_3sd(kP.allSubs');
[~,~,kFp.noOutliers, idx_goodSubs.kFp] = rmv_outliers_3sd(kFp.allSubs');
[~,~,kLm.noOutliers, idx_goodSubs.kLm] = rmv_outliers_3sd(kLm.allSubs');
[~,~,kBias.noOutliers, idx_goodSubs.kBias] = rmv_outliers_3sd(kBias.allSubs');

%% extract meaningful correlation coefficients with all data or after filtering outliers
lac_vars = {'dmpfc','ains','plasma'};
n_lac_vars = length(lac_vars);
bhv_vars = {'HE','HPE','HME',...
    'kEp','kEm','kR','kP','kFp','kLm','kBias'};
n_bhv_vars = length(bhv_vars);
%% extract correlations
for iLac = 1:n_lac_vars
    lac_var_nm = lac_vars{iLac};
    switch lac_var_nm
        case 'dmpfc'
            lac_var = dmPFC_Lac.allSubs;
        case 'ains'
            lac_var = aIns_Lac.allSubs;
        case 'plasma'
            lac_var = plasma_Lac.allSubs;
    end
    for iBhv = 1:n_bhv_vars
        bhv_var_nm = bhv_vars{iBhv};
        switch bhv_var_nm
            case 'HE'
                bhv_var = HE.allSubs;
            case 'HPE'
                bhv_var = HPE.allSubs;
            case 'HME'
                bhv_var = HME.allSubs;
            case 'kEp'
                bhv_var = kEp.allSubs;
            case 'kEm'
                bhv_var = kEm.allSubs;
            case 'kR'
                bhv_var = kR.allSubs;
            case 'kP'
                bhv_var = kP.allSubs;
            case 'kFp'
                bhv_var = kFp.allSubs;
            case 'kLm'
                bhv_var = kLm.allSubs;
            case 'kBias'
                bhv_var = kBias.allSubs;
        end
        corr_nm = [lac_var_nm,'Lac_vs_',bhv_var_nm];
        
        % raw data
        idx_goodSubs.raw.(corr_nm) = ~isnan(lac_var.*bhv_var);
        goodS_tmp1 = idx_goodSubs.raw.(corr_nm);
        r_corr.raw.(corr_nm) = corr(lac_var(goodS_tmp1), bhv_var(goodS_tmp1)); % extract correlation coefficient
        NS_goodS.raw.(corr_nm) = sum(goodS_tmp1); % extract number of subjects
        
        % correct for outliers in each variable included in the
        % correlation
        idx_goodSubs.noOutliers.(corr_nm) = (idx_goodSubs.([lac_var_nm,'_Lac']).*idx_goodSubs.(bhv_var_nm)) == 1;
        goodS_tmp2 = idx_goodSubs.noOutliers.(corr_nm); % create short variable with index of good subjects for practical reasons
        r_corr.noOutliers.(corr_nm) = corr(lac_var(goodS_tmp2), bhv_var(goodS_tmp2)); % extract correlation
        NS_goodS.noOutliers.(corr_nm) = sum(goodS_tmp2); % extract number of good subjects included in the current correlation
        
        % correct for outliers present in dmPFC/dACC AND in aIns (to
        % test specificity of dmPFC/dACC vs aIns)
        if ismember(lac_var_nm,{'dmpfc','ains'})
            outlier_condition_dmPFC_vs_aIns = 'noOutliers_dmPFC_vs_aIns';
            idx_goodSubs.(outlier_condition_dmPFC_vs_aIns).(corr_nm) = (idx_goodSubs.dmpfc_Lac.*idx_goodSubs.ains_Lac.*idx_goodSubs.(bhv_var_nm)) == 1;
            goodS_tmp2 = idx_goodSubs.(outlier_condition_dmPFC_vs_aIns).(corr_nm); % create short variable with index of good subjects for practical reasons
            r_corr.(outlier_condition_dmPFC_vs_aIns).(corr_nm) = corr(lac_var(goodS_tmp2), bhv_var(goodS_tmp2)); % extract correlation
            NS_goodS.(outlier_condition_dmPFC_vs_aIns).(corr_nm) = sum(goodS_tmp2); % extract number of good subjects included in the current correlation
            % also correlate directly dmPFC/dACC with aIns lactate to inform cocor
            r_corr.(outlier_condition_dmPFC_vs_aIns).dmpfcLac_vs_ainsLac = corr(dmPFC_Lac.allSubs(goodS_tmp2), aIns_Lac.allSubs(goodS_tmp2));
        end % filter dmPFC/dACC and aIns only
        
        % correct for outliers present in HPE AND in HME (to
        % test specificity of HPE vs HME)
        if ismember(bhv_var_nm,{'HPE','HME'})
            outlier_condition_HPE_vs_HME = (['noOutliers_',lac_var_nm,'_HPE_vs_HME']);
            idx_goodSubs.(outlier_condition_HPE_vs_HME).(corr_nm) = (idx_goodSubs.([lac_var_nm,'_Lac']).*idx_goodSubs.HPE.*idx_goodSubs.HME) == 1;
            goodS_tmp2 = idx_goodSubs.(outlier_condition_HPE_vs_HME).(corr_nm); % create short variable with index of good subjects for practical reasons
            r_corr.(outlier_condition_HPE_vs_HME).(corr_nm) = corr(lac_var(goodS_tmp2), bhv_var(goodS_tmp2)); % extract correlation
            NS_goodS.(outlier_condition_HPE_vs_HME).(corr_nm) = sum(goodS_tmp2); % extract number of good subjects included in the current correlation
            % also correlate directly HPE with HME to inform cocor
            r_corr.(outlier_condition_HPE_vs_HME).HPE_vs_HME = corr(kEp.allSubs(goodS_tmp2), kEm.allSubs(goodS_tmp2));
        end % filter HPE and HME
        
        % correct for outliers present in kEp AND in other behavioral prm (to
        % test specificity of kEp vs other behavioral prm)
        if ismember(bhv_var_nm,{'kEm','kR','kP','kFp','kLm','kBias'})
            outlier_condition = (['noOutliers_',lac_var_nm,'_kEp_vs_',bhv_var_nm]);
            idx_goodSubs.(outlier_condition).(corr_nm) = (idx_goodSubs.([lac_var_nm,'_Lac']).*idx_goodSubs.kEp.*idx_goodSubs.(bhv_var_nm)) == 1;
            goodS_tmp2 = idx_goodSubs.(outlier_condition).(corr_nm); % create short variable with index of good subjects for practical reasons
            NS_goodS.(outlier_condition).(corr_nm) = sum(goodS_tmp2); % extract number of good subjects included in the current correlation
            
            % extract correlation for kEp (with subjects filtered based on
            % both conditions)
            r_corr.(outlier_condition).([lac_var_nm,'Lac_vs_kEp']) = corr(lac_var(goodS_tmp2), kEp.allSubs(goodS_tmp2)); % extract correlation
            % extract correlation for behavioral parameter 2
            r_corr.(outlier_condition).(corr_nm) = corr(lac_var(goodS_tmp2), bhv_var(goodS_tmp2)); % extract correlation
            
            % also correlate directly kEp with bhv parameter to inform cocor
            r_corr.(outlier_condition).(['kEp_vs_',bhv_var_nm]) = corr(kEp.allSubs(goodS_tmp2), bhv_var(goodS_tmp2));
        end % filter kEp and other behavioral parameter
    end % loop over behavioral parameters
end % loop over lactate measures

%% save for R processing


%% extract all elements to perform all 2-by-2 comparisons in cocor toolbox (http://comparingcorrelations.org/)
prm2_names = {'kEm','kR','kP','kFp','kLm','kBias'};
nPrm_bis = length(prm2_names);

% dmPFC/dACC-lactate/HE vs aIns-lactate/HE
disp(' ');
correl_nm = 'noOutliers_dmPFC_vs_aIns';
disp(['r(dmPFC-Lac - HE) = ',num2str(r_corr.(correl_nm).dmpfcLac_vs_HE)]); % correlation A<=>B
disp(['r(aIns-Lac - HE) = ',num2str(r_corr.(correl_nm).ainsLac_vs_HE)]); % correlation A<=>C
disp(['r(dmPFC-Lac - aIns-Lac) = ',num2str(r_corr.(correl_nm).dmpfcLac_vs_ainsLac)]); % correlation B<=>C
disp(['NS(dmPFC-Lac/HE*aIns-Lac/HE) = ',num2str(NS_goodS.(correl_nm).dmpfcLac_vs_HE)]); % number of subjects included

% dmPFC/dACC-lactate HPE vs HME
disp(' ');
lac_var_nm = 'dmpfc';
correl_nm = ['noOutliers_',lac_var_nm,'_HPE_vs_HME'];
disp(['r(dmPFC-Lac - HPE) = ',num2str(r_corr.(correl_nm).dmpfcLac_vs_HPE)]); % correlation A<=>B
disp(['r(dmPFC-Lac - HME) = ',num2str(r_corr.(correl_nm).dmpfcLac_vs_HME)]); % correlation A<=>C
disp(['r(HPE - HME) = ',num2str(r_corr.(correl_nm).HPE_vs_HME)]); % correlation B<=>C
disp(['NS(dmPFC-Lac-HPE*HME) = ',num2str(NS_goodS.(correl_nm).dmpfcLac_vs_HPE)]); % number of subjects included

% dmPFC/dACC-lactate/HPE vs aIns-lactate/HPE
disp(' ');
correl_nm = 'noOutliers_dmPFC_vs_aIns';
disp(['r(dmPFC-Lac - HPE) = ',num2str(r_corr.(correl_nm).dmpfcLac_vs_HPE)]); % correlation A<=>B
disp(['r(aIns-Lac - HPE) = ',num2str(r_corr.(correl_nm).ainsLac_vs_HPE)]); % correlation A<=>C
disp(['r(dmPFC-Lac - aIns-Lac) = ',num2str(r_corr.(correl_nm).dmpfcLac_vs_ainsLac)]); % correlation B<=>C
disp(['NS(dmPFC-Lac/HPE*aIns-Lac/HPE) = ',num2str(NS_goodS.(correl_nm).dmpfcLac_vs_HPE)]); % number of subjects included

% dmPFC/dACC-lactate kEp vs other parameters
lac_var_nm = 'dmpfc';
for iPrm = 1:nPrm_bis
    prm2_nm = prm2_names{iPrm};
    disp(' ');
    correl_nm = ['noOutliers_',lac_var_nm,'_kEp_vs_',prm2_nm];
    disp(['r(dmPFC-Lac-kEp) = ',num2str(r_corr.(correl_nm).dmpfcLac_vs_kEp)]); % correlation A<=>B
    disp(['r(dmPFC-Lac-',prm2_nm,') = ',num2str(r_corr.(correl_nm).(['dmpfcLac_vs_',prm2_nm]))]); % correlation A<=>C
    disp(['r(kEp-',prm2_nm,') = ',num2str(r_corr.(correl_nm).(['kEp_vs_',prm2_nm]))]); % correlation B<=>C
    disp(['NS(dmPFC-Lac-kEp*',prm2_nm,') = ',num2str(NS_goodS.(correl_nm).(['dmpfcLac_vs_',prm2_nm]))]); % number of subjects included
end % loop on parameters

% dmPFC/dACC-lactate/kEp vs aIns-lactate/kEp
disp(' ');
correl_nm = 'noOutliers_dmPFC_vs_aIns';
disp(['r(dmPFC-Lac - kEp) = ',num2str(r_corr.(correl_nm).dmpfcLac_vs_kEp)]); % correlation A<=>B
disp(['r(aIns-Lac - kEp) = ',num2str(r_corr.(correl_nm).ainsLac_vs_kEp)]); % correlation A<=>C
disp(['r(dmPFC-Lac - aIns-Lac) = ',num2str(r_corr.(correl_nm).dmpfcLac_vs_ainsLac)]); % correlation B<=>C
disp(['NS(dmPFC-Lac/kEp*aIns-Lac/kEp) = ',num2str(NS_goodS.(correl_nm).dmpfcLac_vs_kEp)]); % number of subjects included
