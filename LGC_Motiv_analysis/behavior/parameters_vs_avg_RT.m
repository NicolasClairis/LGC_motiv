function[betas, pval] = parameters_vs_avg_RT()
% [betas, pval] = parameters_vs_avg_RT()
% parameters_vs_avg_RT will look at the correlation between model
% parameters and average reaction times
%
% INPUTS
%
% OUTPUTS
% betas: structure with information about betas from correlations
%
% pval: structure with p.values

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% working directories
computerRoot = LGCM_root_paths;
studyFolder = [computerRoot, study_nm, filesep];

%% extract RT
[avg_RT] = avg_RT_group(study_nm, subject_id, condition);
%% extract parameters
mdlType = 'bayesian';
mdlN = '3';
[prm] = prm_extraction(study_nm, subject_id, mdlType, mdlN);
prm = rmfield(prm,'CID');
prm_names = fieldnames(prm);
nPrm = length(prm_names);

%% correlation
pSize = 30;
for iP = 1:nPrm
    prm_nm = prm_names{iP};
    
    % Ep
    Ep_field_nm = ['RT_Ep_f_',prm_nm];
    Ep_okSubs = ~isnan(prm.(prm_nm).*avg_RT.Ep);
    [betas.(Ep_field_nm),~,stats_Ep_tmp] = glmfit(prm.(prm_nm)(Ep_okSubs),...
        avg_RT.Ep(Ep_okSubs),'normal');
    pval.(Ep_field_nm) = stats_Ep_tmp.p;
    prm_sort_for_fit = sort(prm.(prm_nm));
    avg_RT_fit.Ep = glmval(betas.(Ep_field_nm), prm_sort_for_fit,'identity');
    
    % Em
    Em_field_nm = ['RT_Em_f_',prm_nm];
    Em_okSubs = ~isnan(prm.(prm_nm).*avg_RT.Em);
    [betas.(Em_field_nm),~,stats_Em_tmp] = glmfit(prm.(prm_nm)(Em_okSubs),...
        avg_RT.Em(Em_okSubs),'normal');
    pval.(Em_field_nm) = stats_Em_tmp.p;
    avg_RT_fit.Em = glmval(betas.(Em_field_nm), prm_sort_for_fit,'identity');
    
    % global
    EpEm_field_nm = ['RT_EpEm_f_',prm_nm];
    EpEm_okSubs = ~isnan(prm.(prm_nm).*avg_RT.EpEm);
    [betas.(EpEm_field_nm),~,stats_EpEm_tmp] = glmfit(prm.(prm_nm)(EpEm_okSubs),...
        avg_RT.EpEm(EpEm_okSubs),'normal');
    pval.(EpEm_field_nm) = stats_EpEm_tmp.p;
    avg_RT_fit.EpEm = glmval(betas.(EpEm_field_nm), prm_sort_for_fit,'identity');
    
    % display corresponding figure
    fig;
    % Ep
    subplot(1,3,1); hold on;
    Ep_hdl = scatter(prm.(prm_nm), avg_RT.Ep);
    Ep_hdl.MarkerEdgeColor = 'k';
    fit_Ep_hdl = plot(prm_sort_for_fit, avg_RT_fit.Ep);
    fit_hdl_upgrade(fit_Ep_hdl);
    legend_size(pSize);
    xlabel(prm_nm);
    ylabel('RT Ep');
    
    
    % Em
    subplot(1,3,2); hold on;
    Em_hdl = scatter(prm.(prm_nm), avg_RT.Em);
    Em_hdl.MarkerEdgeColor = 'k';
    fit_Em_hdl = plot(prm_sort_for_fit, avg_RT_fit.Em);
    fit_hdl_upgrade(fit_Em_hdl);
    legend_size(pSize);
    xlabel(prm_nm);
    ylabel('RT Em');
    
    
    % Ep+Em
    subplot(1,3,3); hold on;
    EpEm_hdl = scatter(prm.(prm_nm), avg_RT.EpEm);
    EpEm_hdl.MarkerEdgeColor = 'k';
    fit_EpEm_hdl = plot(prm_sort_for_fit, avg_RT_fit.EpEm);
    fit_hdl_upgrade(fit_EpEm_hdl);
    legend_size(pSize);
    xlabel(prm_nm);
    ylabel('RT Ep+Em');
end % parameter loop
end % function