function[betas, pval, r_corr] = grip_MVC_f_lactate()
% [betas, pval, rho] = grip_MVC_f_lactate()
% grip_MVC_f_lactate will compare plasma and brain levels of lactate to
% maximum voluntary contraction (MVC) force in the physical effort task.
%
% OUTPUTS
% betas: regression estimates for correlations
%
% pval: structure with p.value for correlations
%
% rho: correlation coefficient for correlations

%% subject selection
[~, ~, ~, subject_id, NS] = sub_id;

%% initialize variable of interest
[plasma_Lac, dmPFC_Lac, aIns_Lac,...
    MVC, theoretical_MVC] = deal(NaN(1,NS));

%% load MVC
[~, ~, MVC_subject_id, ~, ~,...
    MVC_all, PCSA_all] = grip_MVC_vs_PCSA(0);
% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    MVC_sub_idx = find(strcmp(sub_nm,MVC_subject_id));
    MVC(iS) = MVC_all(MVC_sub_idx);
    theoretical_MVC(iS) = PCSA_all(MVC_sub_idx);
end % subject loop

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac(:) = metabolites.dmPFC.Lac;
aIns_Lac(:) = metabolites.aIns.Lac;

%% load plasma lactate
[plasma_Lac_struct] = load_plasma_Lac(subject_id);
plasma_Lac(:) = plasma_Lac_struct.Lac;

%% test correlations
% tests with MVC
[r_corr.MVC_f_plasma_Lac, betas.MVC_f_plasma_Lac,...
    pval.MVC_f_plasma_Lac,...
    ~, plasma_Lac_sorted,...
    MVC_fit_f_plasma_Lac] = glm_package(plasma_Lac, MVC, 'normal');
[r_corr.MVC_f_dmPFC_Lac, betas.MVC_f_dmPFC_Lac,...
    pval.MVC_f_dmPFC_Lac,...
    ~, ~,...
    MVC_fit_f_dmPFC_Lac] = glm_package(dmPFC_Lac, MVC, 'normal');
[r_corr.MVC_f_aIns_Lac, betas.MVC_f_aIns_Lac,...
    pval.MVC_f_aIns_Lac,...
    ~, ~,...
    MVC_fit_f_aIns_Lac] = glm_package(aIns_Lac, MVC, 'normal');
% same but with theoretical MVC
[r_corr.theoretical_MVC_f_plasma_Lac, betas.theoretical_MVC_f_plasma_Lac,...
    pval.theoretical_MVC_f_plasma_Lac,...
    ~, plasma_Lac_sorted,...
    theoretical_MVC_fit_f_plasma_Lac] = glm_package(plasma_Lac, theoretical_MVC, 'normal');
[r_corr.theoretical_MVC_f_dmPFC_Lac, betas.theoretical_MVC_f_dmPFC_Lac,...
    pval.theoretical_MVC_f_dmPFC_Lac,...
    ~, ~,...
    theoretical_MVC_fit_f_dmPFC_Lac] = glm_package(dmPFC_Lac, theoretical_MVC, 'normal');
[r_corr.theoretical_MVC_f_aIns_Lac, betas.theoretical_MVC_f_aIns_Lac,...
    pval.theoretical_MVC_f_aIns_Lac,...
    ~, ~,...
    theoretical_MVC_fit_f_aIns_Lac] = glm_package(aIns_Lac, theoretical_MVC, 'normal');

end % function