% function[] = fatigue_f_plasmaMb()
% [] = fatigue_f_plasmaMb()
% fatigue_f_plasmaMb will look at the correlation matrix between the
% different fatigue metrics from questionnaires + behavior and the
% metabolic measures performed in the plasma and in the brain
%
% INPUTS
%
% OUTPUTS
%

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% load fatigue metrics
[fatigue_measures] = fatigue_pool(study_nm, condition, subject_id, NS, genderFilter);
fatigue_vars = fieldnames(fatigue_measures);
fatigue_vars(strcmp(fatigue_vars,'sub_selection')) = [];
n_F_vars = length(fatigue_vars);

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
[metabolites, CRLB, metabolites_bis] = metabolite_load(subject_id);
brain_mb_names = fieldnames(metabolites_bis.dmPFC);
n_brain_mb = length(brain_mb_names);
% add dmPFC/dACC
for iBrain_Mb = 1:n_brain_mb
    brain_mb_nm = brain_mb_names{iBrain_Mb};
    metabolism.(['dmPFC_',brain_mb_nm]) = metabolites_bis.dmPFC.(brain_mb_nm);
end
% add aIns
for iBrain_Mb = 1:n_brain_mb
    brain_mb_nm = brain_mb_names{iBrain_Mb};
    metabolism.(['aIns_',brain_mb_nm]) = metabolites_bis.aIns.(brain_mb_nm);
end

%% correlation matrix
metabolite_names = fieldnames(metabolism);
n_mb_vars = length(metabolite_names);
[corr_F_vs_mb, pval_F_vs_mb] = deal(NaN(n_F_vars, n_mb_vars));
for iF = 1:n_F_vars
    F_nm = fatigue_vars{iF};
    for iMb = 1:n_mb_vars
        mb_nm = metabolite_names{iMb};
        [corr_F_vs_mb(iF, iMb), pval_F_vs_mb(iF, iMb)] = corr(fatigue_measures.(F_nm)', metabolism.(mb_nm)');
        if pval_F_vs_mb(iF,iMb) < 0.05
            signif.([F_nm,'_f_',mb_nm]).r = corr_F_vs_mb(iF, iMb);
            signif.([F_nm,'_f_',mb_nm]).p = pval_F_vs_mb(iF, iMb);
        end
    end % metabolite loop
end % fatigue loop

%% display resulting correlation matrix


% % end % function