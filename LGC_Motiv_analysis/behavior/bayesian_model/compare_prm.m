function[r_corr, pval, NS_goodS, prm, mdlType, mdlN] = compare_prm()
% [r_corr, pval, NS_goodS, prm, mdlType, mdlN] = compare_prm()
% compare_prm will test the correlation between parameters
%
% OUTPUTS
% r_corr
%
% pval:
%
% NS_goodS: structure indicating number of subjects included in each
% correlation as well as for each parameter before and after removing
% outliers
%
% prm: structure with parameter values for each subject
%
% mdlType: bayesian or non-bayesian model?
%
% mdlN: model number

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% load paramaters
[prm, mdlType, mdlN] = prm_extraction(study_nm, subject_id);
prm = rmfield(prm,'CID');
prm_names = fieldnames(prm);
nPrm = length(prm_names);

%% remove outliers
for iPrm = 1:nPrm
    prm_nm = prm_names{iPrm};
    NS_goodS.prm.raw.(prm_nm) = sum(~isnan(prm.(prm_nm)));
    [~, ~, prm.noOutliers.(prm_nm), idx_goodSubs.prm.(prm_nm)] = rmv_outliers_3sd(prm.(prm_nm));
    NS_goodS.prm.noOutliers.(prm_nm) = sum(idx_goodSubs.prm.(prm_nm));
end % parameter loop

%% perform correlations
for iPrm1 = 1:nPrm
    prm_nm1 = prm_names{iPrm1};
    for iPrm2 = 1:nPrm
        prm_nm2 = prm_names{iPrm2};
        corr_nm = [prm_nm1,'_vs_',prm_nm2];
        % correlation for raw data
        [r_corr.raw.(corr_nm), pval.raw.(corr_nm)] = corr(prm.(prm_nm1)', prm.(prm_nm2)');
        NS_goodS.raw.(corr_nm) = NS;
        
        % correlation for data after removing outliers
        % 1) identify subjects ok for both variables
        idx_goodS.(corr_nm) = (idx_goodSubs.prm.(prm_nm1).*idx_goodSubs.prm.(prm_nm2)) == 1;
        NS_goodS.noOutliers.(corr_nm) = sum(idx_goodS.(corr_nm));
        % 2) perform the correlation
        [r_corr.noOutliers.(corr_nm),...
            pval.noOutliers.(corr_nm)] = corr(prm.(prm_nm1)(idx_goodS.(corr_nm))', prm.(prm_nm2)(idx_goodS.(corr_nm))');
    end % parameter loop 2
end % parameter loop 1

%% correlation heatmap
warning('TO ADD');

end % function