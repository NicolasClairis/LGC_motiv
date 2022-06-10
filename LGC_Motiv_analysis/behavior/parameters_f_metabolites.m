%% script to correlate brain metabolites with behavioral scores
%
% designed by N.Clairis - 2022

%% define metabolite and ROI you want to focus on
% ROI
ROIs = {'dmPFC','aIns'};
nROIs = length(ROIs);
ROI_idx = spm_input('Metabolites in which brain area?',1,'m',...
    ROIs,1:nROIs,0);
ROI_nm = ROIs{ROI_idx};
% select metabolite of interest
metabolites = {'Mac','Ala','Asp','PCho','Cr','PCr','GABA',...
    'Gln','Glu','GSH','Gly','Ins','Lac','NAA','Scyllo','Tau',...
    'Asc','Glc','NAAG','GPC','PE','Ser',...
    'NAA_NAAG','Glu_Gln','GPC_PCho','Cr_PCr','Gly_Ins','Gln_div_Glu'};
n_met = length(metabolites);
metabolite_idx = spm_input('Which metabolite to focus on?',1,'m',...
    metabolites,1:n_met,0);
metabolite_nm = metabolites{metabolite_idx};

%% define all subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection('study1', condition);

%% extract all metabolites
[metabolites] = metabolite_load(subject_id);
% focus on metabolite and brain area selected
metabolite_allSubs = metabolites.(ROI_nm).(metabolite_nm);

%% extract behavioral parameters
%% perform behavioral model
figDispGroup = 0;
dispMoneyOrLevels = 'levels';
[betas_fullList, pvalues_fullList] = logitfit_choices_group(figDispGroup, dispMoneyOrLevels);

[prm.kRp, prm.kPp, prm.kEp, prm.kFp,...
    prm.kRm, prm.kPm, prm.kEm, prm.kFm] = deal(NaN(1,NS));
mdl_nm = 'mdl_4';
for iS = 1:NS
    sub_nm = subject_id{iS};
    % extract physical task parameters
    sub_idx = strcmp(betas_fullList.subList,sub_nm);
    prm.kRp(iS) = betas_fullList.Ep.(mdl_nm).kR(sub_idx);
    prm.kPp(iS) = betas_fullList.Ep.(mdl_nm).kP(sub_idx);
    prm.kEp(iS) = betas_fullList.Ep.(mdl_nm).kEffort(sub_idx);
    prm.kFp(iS) = betas_fullList.Ep.(mdl_nm).kFatigue(sub_idx);
    % extract mental task parameters
    prm.kRm(iS) = betas_fullList.Em.(mdl_nm).kR(sub_idx);
    prm.kPm(iS) = betas_fullList.Em.(mdl_nm).kP(sub_idx);
    prm.kEm(iS) = betas_fullList.Em.(mdl_nm).kEffort(sub_idx);
    prm.kFm(iS) = betas_fullList.Em.(mdl_nm).kFatigue(sub_idx);
end % subject list

%% perform GLM and correlation between behavioral parameters and spectroscopy
parameters = {'kRp','kPp','kEp','kFp',...
    'kRm','kPm','kEm','kFm'};
n_prm = length(parameters);
% remove NaN subjects
goodSubs = ~isnan(metabolite_allSubs);
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    var_nm = [prm_nm,'_f_',metabolite_nm];
    [betas.(var_nm), ~,stats_tmp] = glmfit(metabolite_allSubs(goodSubs), prm.(prm_nm)(goodSubs),'normal');
    pval.(var_nm) = stats_tmp.p;
    prm_f_metabolite_fit.(var_nm) = glmval(betas.(var_nm),metabolite_allSubs(goodSubs),'identity');
end % parameters loop

%% display results
pSize = 50;
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
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

end % parameter loop