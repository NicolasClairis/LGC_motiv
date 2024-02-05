function[betas, pval, r_corr] = subj_fatigue_f_lactate()
% [betas, pval, rho] = subj_fatigue_f_lactate()
% subj_fatigue_f_lactate will compare plasma and brain levels of lactate to
% subjective ratings of fatigue across the experiment.
%
% OUTPUTS
% betas: regression estimates for correlations
%
% pval: structure with p.value for correlations
%
% rho: correlation coefficient for correlations

%% subject selection
[study_nm, ~, ~, subject_id, NS] = sub_id;

%% initialize variable of interest
[plasma_Lac, dmPFC_Lac, aIns_Lac] = deal(NaN(1,NS));
fatigue_rating_names = {'preMRS','postMRS','prefMRI','postfMRI',...
    'postfMRI_min_preMRS','postfMRI_min_prefMRI'};
n_fatigue_rtgs = length(fatigue_rating_names);
subj_fatigue = NaN(n_fatigue_rtgs, NS); % 1 measure pre-MRS, 1 post-MRS, 1 pre-fMRI and 1 post-fMRI + delta pre-MRS/post-fMRI + delta pre-fMRI/post-fMRI

%% load subjective fatigue
[~,...
    subj_fatigue(1,:), subj_fatigue(2,:),...
    subj_fatigue(3,:), subj_fatigue(4,:),...
    subj_fatigue(5,:), subj_fatigue(6,:)] = extract_subjective_fatigue_ratings(study_nm, subject_id, NS);

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac(:) = metabolites.dmPFC.Lac;
aIns_Lac(:) = metabolites.aIns.Lac;

%% load plasma lactate
[plasma_Lac_struct] = load_plasma_Lac(subject_id);
plasma_Lac(:) = plasma_Lac_struct.Lac;

%% test correlations
for iF = 1:n_fatigue_rtgs
    fatigue_nm = fatigue_rating_names{iF};
    
    fatigue_plasma_nm = [fatigue_nm,'_fatigue_f_plasma_Lac'];
    [r_corr.(fatigue_plasma_nm), betas.(fatigue_plasma_nm),...
        pval.(fatigue_plasma_nm),...
        ~, plasma_Lac_sorted,...
        subj_fatigue_fit.(fatigue_plasma_nm)] = glm_package(plasma_Lac, subj_fatigue(iF,:), 'normal');
    fatigue_dmPFC_nm = [fatigue_nm,'_fatigue_f_dmPFC_Lac'];
    [r_corr.(fatigue_dmPFC_nm), betas.(fatigue_dmPFC_nm),...
        pval.(fatigue_dmPFC_nm),...
        ~, dmPFC_Lac_sorted,...
        subj_fatigue_fit.(fatigue_dmPFC_nm)] = glm_package(dmPFC_Lac, subj_fatigue(iF,:), 'normal');
    fatigue_aIns_nm = [fatigue_nm,'_fatigue_f_aIns_Lac'];
    [r_corr.(fatigue_aIns_nm), betas.(fatigue_aIns_nm),...
        pval.(fatigue_aIns_nm),...
        ~, aIns_Lac_sorted,...
        subj_fatigue_fit.(fatigue_aIns_nm)] = glm_package(aIns_Lac, subj_fatigue(iF,:), 'normal');
end % loop over fatigue periods
end % function