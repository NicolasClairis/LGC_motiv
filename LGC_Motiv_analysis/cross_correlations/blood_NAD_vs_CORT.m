%% compare blood NAD and salivary interleukins

%% define subjects of interest
condition = subject_condition;
study_nm = 'study1';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load blood NAD
[blood, blood_sub_List] = load_blood_NAD(study_nm);
bloodMb_names = fieldnames(blood);
n_Mb = length(bloodMb_names);
%% load IL
[CORT_data] = load_CORT(study_nm);

%% extract correspondancy between both data
CORT.AUCg = NaN(1,NS);
for iMb = 1:n_Mb
    bloodMb.(bloodMb_names{iMb}) = NaN(1,NS);
end % metabolite loop

for iS = 1:NS
    % subject info
    sub_nm = ['CID',subject_id{iS}];
    sub_blood_idx = strcmp(sub_nm, blood_sub_List);
    sub_CORT_idx = strcmp(sub_nm, CORT_data.CID);
    % extract blood info
    for iMb = 1:n_Mb
        bloodMb_nm = bloodMb_names{iMb};
        bloodMb.(bloodMb_nm)(iS) = blood.(bloodMb_nm)(sub_blood_idx);
    end
    % extract CORT data
    CORT.AUCg(iS) = CORT_data.AUCg(sub_CORT_idx);
end % subject loop

%% correlation and figure
lWidth = 3;
pSize = 25;
black = 'k';
orange = [254 75 3]./255;
figCORT1.AUCg = fig; j_fig1 = 0;
figCORT2.AUCg = fig; j_fig2 = 0;
CORT_nm = 'AUCg';
for iMb = 1:n_Mb
    bloodMb_nm = bloodMb_names{iMb};
    
    switch bloodMb_nm
            case {'Nam','NMN','NR','NAD',...
                    'NADH','NADP','NADPH','MeNam',...
                    'MeXPY'}
                figure(figCORT1.AUCg);
                j_fig1 = j_fig1 + 1;
                subplot(3,4,j_fig1);
            case {'NAD_div_NADH',...
                    'NADP_div_NADPH',...
                    'total_NAD_precursors',...
                    'total_NAD',...
                    'total_NAD_with_precursors',...
                    'total_NAD_with_byproducts',...
                    'total_NAD_byproducts'}
                figure(figCORT2.AUCg);
                j_fig2 = j_fig2 + 1;
                subplot(3,3,j_fig2);
    end
        
    %% test correlations
    goodSubs.(CORT_nm) = (~isnan(CORT.(CORT_nm))).*(~isnan(bloodMb.(bloodMb_nm))) == 1;
    corr_nm.(CORT_nm) = ['CORT_',CORT_nm,'_f_',bloodMb_nm];
    [beta.(corr_nm.(CORT_nm)),~,stats_tmp] =...
        glmfit(bloodMb.(bloodMb_nm)(goodSubs.(CORT_nm)),...
        CORT.(CORT_nm)(goodSubs.(CORT_nm)), 'normal');
    pval.(corr_nm.(CORT_nm)) = stats_tmp.p;
    % store significant p.values for slope
    if stats_tmp.p(2) < 0.05
        pval.signif.(corr_nm.(CORT_nm)) = stats_tmp.p(2);
    elseif stats_tmp.p(2) > 0.05 && stats_tmp.p(2) < 0.1
        pval.almostSignif.(corr_nm.(CORT_nm)) = stats_tmp.p(2);
    end
    bloodMb_sort = sort(bloodMb.(bloodMb_nm)(goodSubs.(CORT_nm)));
    CORT_fit = glmval(beta.(corr_nm.(CORT_nm)), bloodMb_sort, 'identity');
    %% figure
    hold on;
    scat_hdl = scatter(bloodMb.(bloodMb_nm)(goodSubs.(CORT_nm)),...
        CORT.(CORT_nm)(goodSubs.(CORT_nm)));
    plot_hdl = plot(bloodMb_sort, CORT_fit);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    plot_hdl.Color = orange;
    plot_hdl.LineStyle = '--';
    [blood_labelname] = blood_label(bloodMb_nm);
    xlabel(blood_labelname);
    ylabel(['CORT ',(CORT_nm),' (μg/dL)']);
    legend_size(pSize);
end % metabolite loop