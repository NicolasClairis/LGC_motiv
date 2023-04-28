% script to look at correlation between metabolites and blood NAD/NADH

%% define subjects
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% extract all metabolites
[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);
metabolite_nm_bis = metab_div_rnm(metabolite_nm);
full_mb_nm = [MRS_ROI_nm,'_',metabolite_nm];
full_mb_nm_bis = [MRS_ROI_nm,' ',metabolite_nm_bis];
%% extract all blood data
[bloodTable, blood_NAD_sub_List] = load_blood_NAD(study_nm);
bloodMb_names = fieldnames(bloodTable);
n_BloodPrm = length(bloodMb_names);
for iBlood = 1:n_BloodPrm
    bloodMb.(bloodMb_names{iBlood}) = NaN(1,NS);
end % metabolite loop
% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    
    %% load blood
    sub_blood_idx = strcmp(sub_nm_bis, blood_NAD_sub_List);
    % extract blood info
    for iBlood = 1:n_BloodPrm
        bloodMb_nm = bloodMb_names{iBlood};
        bloodMb.(bloodMb_nm)(iS) = bloodTable.(bloodMb_nm)(sub_blood_idx);
    end
end % subject loop

%% test correlations
for iBlood = 1:n_BloodPrm
    bloodMb_nm = bloodMb_names{iBlood};
    corr_nm = [full_mb_nm,'_f_',bloodMb_nm];
    goodSubs.(corr_nm) = (~isnan(metabolite_allSubs)).*(~isnan(bloodMb.(bloodMb_nm))) == 1;
    [beta.(corr_nm),~,stats_tmp] =...
        glmfit(bloodMb.(bloodMb_nm)(goodSubs.(corr_nm)),...
        metabolite_allSubs(goodSubs.(corr_nm)), 'normal');
    pval.(corr_nm) = stats_tmp.p;
    % store significant p.values for slope
    if stats_tmp.p(2) < 0.05
        pval.signif.(corr_nm) = stats_tmp.p(2);
    elseif stats_tmp.p(2) > 0.05 && stats_tmp.p(2) < 0.1
        pval.almostSignif.(corr_nm) = stats_tmp.p(2);
    end
    bloodMb_sort.(corr_nm) = sort(bloodMb.(bloodMb_nm)(goodSubs.(corr_nm)));
    metabolite_fit.(corr_nm) = glmval(beta.(corr_nm), bloodMb_sort.(corr_nm), 'identity');
end % blood loop

%% correlation and figure
% if figDisp == 1
lWidth = 3;
pSize = 25;
black = 'k';
orange = [254 75 3]./255;
%% show results
fig1 = fig; j_fig1 = 0;
fig2 = fig; j_fig2 = 0;
for iBlood = 1:n_BloodPrm
    bloodMb_nm = bloodMb_names{iBlood};
    corr_nm = [full_mb_nm,'_f_',bloodMb_nm];
    
    switch bloodMb_nm
        case {'Nam','NMN','NR','NAD',...
                'NADH','NADP','NADPH','MeNam',...
                'MeXPY'}
            figure(fig1);
            j_fig1 = j_fig1 + 1;
            subplot(3,4,j_fig1);
        case {'NAD_div_NADH',...
                'NADP_div_NADPH',...
                'total_NAD_precursors',...
                'total_NAD',...
                'total_NAD_with_precursors',...
                'total_NAD_with_byproducts',...
                'total_NAD_byproducts'}
            figure(fig2);
            j_fig2 = j_fig2 + 1;
            subplot(3,3,j_fig2);
    end
    %% figure
    hold on;
    scat_hdl = scatter(bloodMb.(bloodMb_nm)(goodSubs.(corr_nm)),...
        metabolite_allSubs(goodSubs.(corr_nm)));
    plot_hdl = plot(bloodMb_sort.(corr_nm), metabolite_fit.(corr_nm));
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    plot_hdl.Color = orange;
    plot_hdl.LineStyle = '--';
    [blood_labelname] = blood_label(bloodMb_nm);
    xlabel(blood_labelname);
    ylabel(full_mb_nm_bis);
    legend_size(pSize);
end % metabolite loop
% end % fig disp