% function[] = allMetabolites_f_inflammation()

%% define subjects to include
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, 'all');

%% load metabolites for all individuals and all brain areas
[metabolites] = metabolite_load(subject_id);
switch study_nm
    case 'study1'
        MRS_ROIs = {'dmPFC','aIns'};
    case 'study2'
        error('not ready yet');
end
nROIs = length(MRS_ROIs);
for iROI = 1:nROIs
    metabolite_names.(MRS_ROIs{iROI}) = fieldnames(metabolites.(MRS_ROIs{iROI}));
    n_metabolites.(MRS_ROIs{iROI}) = length(metabolite_names.(MRS_ROIs{iROI}));
end % roi loop

%% load inflammatory markers
[IL_data] = load_IL(study_nm);

%% extract correlation data
[IL1b, IL6, IL18, IL_sum] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_line = strcmp(['CID',sub_nm], IL_data.CID);
    IL1b(iS) = IL_data.IL1b(sub_line);
    IL6(iS) = IL_data.IL6(sub_line);
    IL18(iS) = IL_data.IL18(sub_line);
    IL_sum(iS) = IL1b(iS) + IL6(iS) + IL18(iS);
end % subject loop
%% correlation and figure
pSize = 50;
lWidth = 3;
black = [0 0 0];
grey = [143 143 143]./255;

for iROI = 1:nROIs
    MRS_ROI_nm = MRS_ROIs{iROI};
    for iMb = 1:n_metabolites.(MRS_ROI_nm)
        metab_nm = metabolite_names.(MRS_ROI_nm){iMb};
        
        %% correlations
        IL1b_goodSubs = ~isnan(IL1b).*~isnan(metabolites.(MRS_ROI_nm).(metab_nm)) == 1;
        [b_IL1b.(MRS_ROI_nm).(metab_nm),~,stats_tmp] = glmfit(IL1b(IL1b_goodSubs),...
            metabolites.(MRS_ROI_nm).(metab_nm)(IL1b_goodSubs), 'normal');
        pval.IL1b.(MRS_ROI_nm).(metab_nm) = stats_tmp.p;
        IL1b_fit = sort(IL1b(IL1b_goodSubs));
        IL1b_metab_fit = glmval(b_IL1b.(MRS_ROI_nm).(metab_nm), IL1b_fit, 'identity');
        if stats_tmp.p(2) < 0.05
            pval.signif.IL1b.(MRS_ROI_nm).(metab_nm) = stats_tmp.p(2);
        end
        
        IL6_goodSubs = ~isnan(IL6).*~isnan(metabolites.(MRS_ROI_nm).(metab_nm)) == 1;
        [b_IL6.(MRS_ROI_nm).(metab_nm),~,stats_tmp] = glmfit(IL6(IL6_goodSubs),...
            metabolites.(MRS_ROI_nm).(metab_nm)(IL6_goodSubs), 'normal');
        pval.IL6.(MRS_ROI_nm).(metab_nm) = stats_tmp.p;
        IL6_fit = sort(IL6(IL6_goodSubs));
        IL6_metab_fit = glmval(b_IL6.(MRS_ROI_nm).(metab_nm), IL6_fit, 'identity');
        if stats_tmp.p(2) < 0.05
            pval.signif.IL6.(MRS_ROI_nm).(metab_nm) = stats_tmp.p(2);
        end
        
        IL18_goodSubs = ~isnan(IL18).*~isnan(metabolites.(MRS_ROI_nm).(metab_nm)) == 1;
        [b_IL18.(MRS_ROI_nm).(metab_nm),~,stats_tmp] = glmfit(IL18(IL18_goodSubs),...
            metabolites.(MRS_ROI_nm).(metab_nm)(IL18_goodSubs), 'normal');
        pval.IL18.(MRS_ROI_nm).(metab_nm) = stats_tmp.p;
        IL18_fit = sort(IL18(IL18_goodSubs));
        IL18_metab_fit = glmval(b_IL18.(MRS_ROI_nm).(metab_nm), IL18_fit, 'identity');
        if stats_tmp.p(2) < 0.05
            pval.signif.IL18.(MRS_ROI_nm).(metab_nm) = stats_tmp.p(2);
        end
        
        ILsum_goodSubs = ~isnan(IL_sum).*~isnan(metabolites.(MRS_ROI_nm).(metab_nm)) == 1;
        [b_ILsum.(MRS_ROI_nm).(metab_nm),~,stats_tmp] = glmfit(IL_sum(ILsum_goodSubs),...
            metabolites.(MRS_ROI_nm).(metab_nm)(ILsum_goodSubs), 'normal');
        pval.ILsum.(MRS_ROI_nm).(metab_nm) = stats_tmp.p;
        ILsum_fit = sort(IL_sum(ILsum_goodSubs));
        ILsum_metab_fit = glmval(b_ILsum.(MRS_ROI_nm).(metab_nm), ILsum_fit, 'identity');
        if stats_tmp.p(2) < 0.05
            pval.signif.ILsum.(MRS_ROI_nm).(metab_nm) = stats_tmp.p(2);
        end
    
        %% figure
        fig;
        
        subplot(2,2,1);
        hold on;
        IL1b_hdl = scatter(IL1b, metabolites.(MRS_ROI_nm).(metab_nm));
        fit_hdl = plot(IL1b_fit, IL1b_metab_fit);
        IL1b_hdl.LineWidth = lWidth;
        IL1b_hdl.MarkerEdgeColor = black;
        fit_hdl.LineStyle = '--';
        fit_hdl.Color = grey;
        xlabel('IL-1b (pg/mL)');
        ylabel([MRS_ROI_nm,' - ',metab_nm]);
        legend_size(pSize);
        
        subplot(2,2,2);
        hold on;
        IL6_hdl = scatter(IL6, metabolites.(MRS_ROI_nm).(metab_nm));
        fit_hdl = plot(IL6_fit, IL6_metab_fit);
        IL6_hdl.LineWidth = lWidth;
        IL6_hdl.MarkerEdgeColor = black;
        fit_hdl.LineStyle = '--';
        fit_hdl.Color = grey;
        xlabel('IL-6 (pg/mL)');
        ylabel([MRS_ROI_nm,' - ',metab_nm]);
        legend_size(pSize);
        
        subplot(2,2,3);
        hold on;
        IL18_hdl = scatter(IL18, metabolites.(MRS_ROI_nm).(metab_nm));
        fit_hdl = plot(IL18_fit, IL18_metab_fit);
        IL18_hdl.LineWidth = lWidth;
        IL18_hdl.MarkerEdgeColor = black;
        fit_hdl.LineStyle = '--';
        fit_hdl.Color = grey;
        xlabel('IL-18 (pg/mL)');
        ylabel([MRS_ROI_nm,' - ',metab_nm]);
        legend_size(pSize);
        
        subplot(2,2,4);
        hold on;
        ILsum_hdl = scatter(IL_sum, metabolites.(MRS_ROI_nm).(metab_nm));
        fit_hdl = plot(ILsum_fit, ILsum_metab_fit);
        ILsum_hdl.LineWidth = lWidth;
        ILsum_hdl.MarkerEdgeColor = black;
        fit_hdl.LineStyle = '--';
        fit_hdl.Color = grey;
        xlabel('IL-1b + IL-6 + IL-18 (pg/mL)');
        ylabel([MRS_ROI_nm,' - ',metab_nm]);
        legend_size(pSize);
        
    end % metabolite loop
end % ROI loop
% end % function