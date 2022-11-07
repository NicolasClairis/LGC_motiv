% function[] = behavPrm_vs_inflammation()
% will check correlation between behavioral parameters and inflammatory
% markers (interleukins)

%% define subject list
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load interleukins
[IL_data] = load_IL;

%% load behavioral parameter
[prm, mdlType, mdlN] = prm_extraction(subject_id, [], []);
% remove CID from parameters
prm_CID = prm.CID;
prm = rmfield(prm, 'CID');
% extract full list of parameters of the model selected
parameter_names = fieldnames(prm);
nPrm = length(parameter_names);
%% correlate both
[IL1b, IL6, IL18, IL_sum] = deal(NaN(1,NS));
for iPrm = 1:nPrm
    prm_nm = parameter_names{iPrm};
    parameter.(prm_nm) = deal(NaN(1,NS));
end % parameter loop
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    behavPrm_CID = strcmp(sub_nm, prm_CID);
    IL_CID = strcmp(['CID',sub_nm], IL_data.CID);
    
    % extract the data
    for iPrm = 1:nPrm
        prm_nm = parameter_names{iPrm};
        parameter.(prm_nm)(iS) = prm.(prm_nm)(behavPrm_CID);
    end % parameter loop
    IL1b(iS) = IL_data.IL1b(IL_CID);
    IL6(iS) = IL_data.IL6(IL_CID);
    IL18(iS) = IL_data.IL18(IL_CID);
    IL_sum(iS) = IL1b(iS) + IL6(iS) + IL18(iS);
end % subject loop

%% correlations and figures
pSize = 50;
lWidth = 3;
black = [0 0 0];
grey = [143 143 143]./255;

for iPrm = 1:nPrm
    prm_nm = parameter_names{iPrm};
    
    %% perform correlations
    IL1b_goodSubs = ~isnan(IL1b).*~isnan(parameter.(prm_nm)) == 1;
    [b_IL1b.(prm_nm),~,stats_tmp] = glmfit(IL1b(IL1b_goodSubs),...
        parameter.(prm_nm)(IL1b_goodSubs), 'normal');
    pval.IL1b.(prm_nm) = stats_tmp.p;
    IL1b_fit = sort(IL1b(IL1b_goodSubs));
    IL1b_prm_fit = glmval(b_IL1b.(prm_nm), IL1b_fit, 'identity');
    
    IL6_goodSubs = ~isnan(IL6).*~isnan(parameter.(prm_nm)) == 1;
    [b_IL6.(prm_nm),~,stats_tmp] = glmfit(IL1b(IL6_goodSubs),...
        parameter.(prm_nm)(IL6_goodSubs), 'normal');
    pval.IL6.(prm_nm) = stats_tmp.p;
    IL6_fit = sort(IL6(IL6_goodSubs));
    IL6_prm_fit = glmval(b_IL6.(prm_nm), IL6_fit, 'identity');
    
    IL18_goodSubs = ~isnan(IL18).*~isnan(parameter.(prm_nm)) == 1;
    [b_IL18.(prm_nm),~,stats_tmp] = glmfit(IL18(IL18_goodSubs),...
        parameter.(prm_nm)(IL18_goodSubs), 'normal');
    pval.IL18.(prm_nm) = stats_tmp.p;
    IL18_fit = sort(IL18(IL18_goodSubs));
    IL18_prm_fit = glmval(b_IL18.(prm_nm), IL18_fit, 'identity');
    
    ILsum_goodSubs = ~isnan(IL_sum).*~isnan(parameter.(prm_nm)) == 1;
    [b_ILsum.(prm_nm),~,stats_tmp] = glmfit(IL_sum(ILsum_goodSubs),...
        parameter.(prm_nm)(ILsum_goodSubs), 'normal');
    pval.IL_sum.(prm_nm) = stats_tmp.p;
    ILsum_fit = sort(IL_sum(ILsum_goodSubs));
    ILsum_prm_fit = glmval(b_ILsum.(prm_nm), ILsum_fit, 'identity');
    
    %% figure
    fig;
    
    % IL 1b
    subplot(2,2,1);
    hold on;
    scat_hdl = scatter(IL1b(IL1b_goodSubs),...
        parameter.(prm_nm)(IL1b_goodSubs));
    fit_hdl = plot(IL1b_fit, IL1b_prm_fit);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl.LineStyle = '--';
    fit_hdl.Color = grey;
    ylabel(prm_nm);
    xlabel('IL-1b (pg/mL)');
    legend_size(pSize);
    
    % IL 6
    subplot(2,2,2);
    hold on;
    scat_hdl = scatter(IL6(IL6_goodSubs),...
        parameter.(prm_nm)(IL6_goodSubs));
    fit_hdl = plot(IL6_fit, IL6_prm_fit);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl.LineStyle = '--';
    fit_hdl.Color = grey;
    ylabel(prm_nm);
    xlabel('IL-6 (pg/mL)');
    legend_size(pSize);
    
    % IL 18
    subplot(2,2,3);
    hold on;
    scat_hdl = scatter(IL18(IL18_goodSubs),...
        parameter.(prm_nm)(IL18_goodSubs));
    fit_hdl = plot(IL18_fit, IL18_prm_fit);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl.LineStyle = '--';
    fit_hdl.Color = grey;
    ylabel(prm_nm);
    xlabel('IL-18 (pg/mL)');
    legend_size(pSize);
    
    % sum IL
    subplot(2,2,4);
    hold on;
    scat_hdl = scatter(IL_sum(ILsum_goodSubs),...
        parameter.(prm_nm)(ILsum_goodSubs));
    fit_hdl = plot(ILsum_fit, ILsum_prm_fit);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl.LineStyle = '--';
    fit_hdl.Color = grey;
    ylabel(prm_nm);
    xlabel('IL-1b + IL6 + IL18 (pg/mL)');
    legend_size(pSize);
    
end % parameter loop
% end % function