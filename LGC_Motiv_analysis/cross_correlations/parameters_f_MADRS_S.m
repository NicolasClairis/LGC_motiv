%% script to correlate basic parameters from the behavioral model
% with MADRS-S scores
%
% designed by N.Clairis - 2022

%% perform behavioral model
computerRoot = LGCM_root_paths;
study_nm = 'study1';
figDispGroup = 0;
figDispIndiv = 0;
dispMoneyOrLevels = 'levels';
n_NV_bins = 6;
n_trialN_bins = 6;
[betas_fullList, pvalues_fullList] = logitfit_choices_group(computerRoot, study_nm,...
    figDispGroup, figDispIndiv, dispMoneyOrLevels, n_NV_bins, n_trialN_bins);

%% load MADRS-S scores
MADRS_S_fullList = getfield(load([computerRoot,filesep,...
    study_nm,filesep,...
    'MADRS_S_temporary.mat'],'MADRS_S'),'MADRS_S');

%% extract relevant subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

[MADRS_S_score,...
    prm.kRp, prm.kPp, prm.kEp, prm.kFp,...
    prm.kRm, prm.kPm, prm.kEm, prm.kFm] = deal(NaN(1,NS));
mdl_nm = 'mdl_4';
for iS = 1:NS
    sub_nm = subject_id{iS};
    MADRS_S_score(iS) = MADRS_S_fullList.MADRS_S_score(strcmp(MADRS_S_fullList.subject,sub_nm));
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
%% perform GLM and correlation between behavioral parameters and MADRS-S
parameters = {'kRp','kPp','kEp','kFp',...
    'kRm','kPm','kEm','kFm'};
n_prm = length(parameters);
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    var_nm = [prm_nm,'_f_MADRS'];
    [betas.(var_nm), ~,stats_tmp] = glmfit(MADRS_S_score, prm.(prm_nm),'normal');
    pval.(var_nm) = stats_tmp.p;
    prm_f_MADRS_S_fit.(var_nm) = glmval(betas.(var_nm),MADRS_S_score,'identity');
end % parameters loop

%% display results
pSize = 50;
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    fig;
    hdl=scatter(MADRS_S_score, prm.(prm_nm),...
        'MarkerEdgeColor','k','SizeData',100,'LineWidth',3);
    hold on;
    [MADRS_S_score_sorted, idx_sorted] = sort(MADRS_S_score);
    plot(MADRS_S_score_sorted, prm_f_MADRS_S_fit.(var_nm)(idx_sorted),...
        'LineStyle','--','LineWidth',3);
%     ylim_vals = ylim();
%     line([4 4],ylim_vals,...
%         'LineStyle','-','LineWidth',3,'Color','k');
    xlabel('MADRS-S score');
    ylabel(prm_nm);
    legend_size(pSize);

end % parameter loop