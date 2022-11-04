%% script to correlate brain metabolites with behavioral scores
%
% designed by N.Clairis - 2022

%% define all subjects
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% define metabolite and ROI you want to focus on
[metabolite_allSubs, ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);

%% extract behavioral parameters
prm = prm_extraction(subject_id);
parameters = fieldnames(prm);

%% perform GLM and correlation between behavioral parameters and spectroscopy
n_prm = length(parameters);
% remove NaN subjects
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    if ~strcmp(prm_nm,'CID') % no sense in this case
        goodSubs = (~isnan(metabolite_allSubs).*~isnan(prm.(prm_nm)))==true;
        var_nm = [prm_nm,'_f_',metabolite_nm];
        [betas.(var_nm), ~,stats_tmp] = glmfit(metabolite_allSubs(goodSubs), prm.(prm_nm)(goodSubs),'normal');
        pval.(var_nm) = stats_tmp.p;
        prm_f_metabolite_fit.(var_nm) = glmval(betas.(var_nm),metabolite_allSubs(goodSubs),'identity');
    end
end % parameters loop

%% display results
pSize = 50;
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    if ~strcmp(prm_nm,'CID') % no sense in this case
        var_nm = [prm_nm,'_f_',metabolite_nm];
        fig;
        hdl=scatter(metabolite_allSubs, prm.(prm_nm),...
            'MarkerEdgeColor','k','SizeData',100,'LineWidth',3);
        hold on;
        [metabolite_allSubs_sorted, idx_sorted] = sort(metabolite_allSubs(goodSubs));
        plot(metabolite_allSubs_sorted, prm_f_metabolite_fit.(var_nm)(idx_sorted),...
            'LineStyle','--','LineWidth',3);
        xlabel([ROI_nm,' ',metabolite_nm]);
        ylabel(prm_nm);
        legend_size(pSize);
    end
end % parameter loop