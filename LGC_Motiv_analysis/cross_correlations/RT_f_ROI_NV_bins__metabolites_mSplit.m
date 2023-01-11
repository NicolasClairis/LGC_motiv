%% extract RT in function of ROI and input parameter
RT_f_ROI_NV_bins;

%% load group indexes depending on the level of metabolites
[met_subs.low, met_subs.high, metabolite_nm, ROI_nm] = medSplit_metabolites(study_nm, subject_id);
met_groups = {'low','high'};
n_met_grps = length(met_groups);
%% average data according to which group (low/high metabolite)
% basic data
for iGrp = 1:n_met_grps
    grp_nm = met_groups{iGrp};
    grp_nm_bis = [grp_nm,'_',metabolite_nm];
    % extract data
    [m_input_f_input.(grp_nm_bis),...
        sem_input_f_input.(grp_nm_bis)] = mean_sem_sd(input_f_input_bin(:,met_subs.(grp_nm)), 2);
    [m_ROI_f_input.(grp_nm_bis),...
        sem_ROI_f_input.(grp_nm_bis)] = mean_sem_sd(ROI_f_input_bin(:,met_subs.(grp_nm)), 2);
    [m_RT_f_input.(grp_nm_bis),...
        sem_RT_f_input.(grp_nm_bis)] = mean_sem_sd(RT_f_input_bin(:,met_subs.(grp_nm)), 2);
    % low ROI median split
    [m_input_f_input_low_ROI.(grp_nm_bis),...
        sem_input_f_input_low_ROI.(grp_nm_bis)] = mean_sem_sd(input_f_input_low_ROI(:,met_subs.(grp_nm)), 2);
    [m_ROI_f_input_low_ROI.(grp_nm_bis),...
        sem_ROI_f_input_low_ROI.(grp_nm_bis)] = mean_sem_sd(ROI_f_input_low_ROI(:,met_subs.(grp_nm)), 2);
    [m_RT_f_input_low_ROI.(grp_nm_bis),...
        sem_RT_f_input_low_ROI.(grp_nm_bis)] = mean_sem_sd(RT_f_input_low_ROI(:,met_subs.(grp_nm)), 2);
    % high ROI median split
    [m_input_f_input_high_ROI.(grp_nm_bis),...
        sem_input_f_input_high_ROI.(grp_nm_bis)] = mean_sem_sd(input_f_input_high_ROI(:,met_subs.(grp_nm)), 2);
    [m_ROI_f_input_high_ROI.(grp_nm_bis),...
        sem_ROI_f_input_high_ROI.(grp_nm_bis)] = mean_sem_sd(ROI_f_input_high_ROI(:,met_subs.(grp_nm)), 2);
    [m_RT_f_input_high_ROI.(grp_nm_bis),...
        sem_RT_f_input_high_ROI.(grp_nm_bis)] = mean_sem_sd(RT_f_input_high_ROI(:,met_subs.(grp_nm)), 2);
end % loop median split

%% figures
if dispFig == true
    lWidth = 3;
    black = [0 0 0];
    grey = [143 143 143]./255;
    purple = [153 142 195]./255;
    orange = [241 163 64]./255;
    pSize = 50;
    input_prm_nm = strrep(input_prm_nm,'_',' ');
    low_met_nm = ['low_',metabolite_nm];
    high_met_nm = ['high_',metabolite_nm];
    low_met_nm_bis = ['low ',metabolite_nm];
    high_met_nm_bis = ['high ',metabolite_nm];
    
    % look at the general figure (RT = f(inputs), ROI=f(inputs)
    fig;
    % RT = f(inputs)
    subplot(1,2,1);
    hold on;
    gal_data_hdl_low_met = errorbar(m_input_f_input.(low_met_nm),...
        m_RT_f_input.(low_met_nm),...
        sem_RT_f_input.(low_met_nm));
    gal_data_hdl_high_met = errorbar(m_input_f_input.(high_met_nm),...
        m_RT_f_input.(high_met_nm),...
        sem_RT_f_input.(high_met_nm));
    gal_data_hdl_low_met.LineStyle = '--';
    gal_data_hdl_low_met.LineWidth = lWidth;
    gal_data_hdl_low_met.MarkerEdgeColor = grey;
    gal_data_hdl_high_met.LineStyle = '-';
    gal_data_hdl_high_met.LineWidth = lWidth;
    gal_data_hdl_high_met.MarkerEdgeColor = black;
    legend([gal_data_hdl_high_met, gal_data_hdl_low_met],...
        {high_met_nm_bis, low_met_nm_bis});
    legend('boxoff');
    xlabel({input_prm_nm; full_bhv_taskName});
    ylabel({'RT (s)'; full_bhv_taskName});
    legend_size(pSize);

    subplot(1,2,2);
    % ROI = f(inputs)
    hold on;
    ROI_f_input_hdl_low_met = errorbar(m_input_f_input.(low_met_nm),...
        m_ROI_f_input.(low_met_nm),...
        sem_ROI_f_input.(low_met_nm));
    ROI_f_input_hdl_high_met = errorbar(m_input_f_input.(high_met_nm),...
        m_ROI_f_input.(high_met_nm),...
        sem_ROI_f_input.(high_met_nm));
    ROI_f_input_hdl_low_met.LineStyle = '--';
    ROI_f_input_hdl_low_met.LineWidth = lWidth;
    ROI_f_input_hdl_low_met.MarkerEdgeColor = grey;
    ROI_f_input_hdl_high_met.LineStyle = '-';
    ROI_f_input_hdl_high_met.LineWidth = lWidth;
    ROI_f_input_hdl_high_met.MarkerEdgeColor = black;
    legend([ROI_f_input_hdl_high_met, ROI_f_input_hdl_low_met],...
        {high_met_nm_bis, low_met_nm_bis});
    legend('boxoff');
    xlabel({input_prm_nm; full_bhv_taskName});
    ylabel({[ROI_short_nm,' BOLD'];...
        full_ROI_taskName});
    legend_size(pSize);

    %% sanity check: ROI = f(inputs) + ROI median split (did it work?)
    fig;
    hold on;
    ROI_f_input_low_ROI_low_met_hdl = errorbar(m_input_f_input_low_ROI.(low_met_nm),...
        m_ROI_f_input_low_ROI.(low_met_nm),...
        sem_ROI_f_input_low_ROI.(low_met_nm));
    ROI_f_input_high_ROI_low_met_hdl = errorbar(m_input_f_input_high_ROI.(low_met_nm),...
        m_ROI_f_input_high_ROI.(low_met_nm),...
        sem_ROI_f_input_high_ROI.(low_met_nm));
    ROI_f_input_low_ROI_high_met_hdl = errorbar(m_input_f_input_low_ROI.(high_met_nm),...
        m_ROI_f_input_low_ROI.(high_met_nm),...
        sem_ROI_f_input_low_ROI.(high_met_nm));
    ROI_f_input_high_ROI_high_met_hdl = errorbar(m_input_f_input_high_ROI.(high_met_nm),...
        m_ROI_f_input_high_ROI.(high_met_nm),...
        sem_ROI_f_input_high_ROI.(high_met_nm));
    ROI_f_input_low_ROI_low_met_hdl.LineStyle = '--';
    ROI_f_input_low_ROI_low_met_hdl.LineWidth = lWidth;
    ROI_f_input_low_ROI_low_met_hdl.MarkerEdgeColor = purple;
    ROI_f_input_high_ROI_low_met_hdl.LineStyle = '--';
    ROI_f_input_high_ROI_low_met_hdl.LineWidth = lWidth;
    ROI_f_input_high_ROI_low_met_hdl.MarkerEdgeColor = orange;
    ROI_f_input_low_ROI_high_met_hdl.LineStyle = '-';
    ROI_f_input_low_ROI_high_met_hdl.LineWidth = lWidth;
    ROI_f_input_low_ROI_high_met_hdl.MarkerEdgeColor = purple;
    ROI_f_input_high_ROI_high_met_hdl.LineStyle = '-';
    ROI_f_input_high_ROI_high_met_hdl.LineWidth = lWidth;
    ROI_f_input_high_ROI_high_met_hdl.MarkerEdgeColor = orange;
    xlabel({input_prm_nm; full_bhv_taskName});
    ylabel({[ROI_short_nm,' BOLD'];...
        full_ROI_taskName});
    legend([ROI_f_input_low_ROI_low_met_hdl, ROI_f_input_high_ROI_low_met_hdl,...
        ROI_f_input_low_ROI_high_met_hdl, ROI_f_input_high_ROI_high_met_hdl],...
        {['low ',ROI_short_nm,' - low ',metabolite_nm],...
        ['high ',ROI_short_nm,' - low ',metabolite_nm],...
        ['low ',ROI_short_nm,' - high ',metabolite_nm],...
        ['high ',ROI_short_nm,' - high ',metabolite_nm]});
    legend('boxoff');
    legend_size(pSize);

    %% main figure: RT = f(inputs) + ROI median split
    fig;
    hold on;
    RT_f_inputs_low_ROI_data_low_met_hdl = errorbar(m_input_f_input_low_ROI.(low_met_nm),...
        m_RT_f_input_low_ROI.(low_met_nm),...
        sem_RT_f_input_low_ROI.(low_met_nm));
    RT_f_inputs_high_ROI_data_low_met_hdl = errorbar(m_input_f_input_high_ROI.(low_met_nm),...
        m_RT_f_input_high_ROI.(low_met_nm),...
        sem_RT_f_input_high_ROI.(low_met_nm));
    RT_f_inputs_low_ROI_data_high_met_hdl = errorbar(m_input_f_input_low_ROI.(high_met_nm),...
        m_RT_f_input_low_ROI.(high_met_nm),...
        sem_RT_f_input_low_ROI.(high_met_nm));
    RT_f_inputs_high_ROI_data_high_met_hdl = errorbar(m_input_f_input_high_ROI.(high_met_nm),...
        m_RT_f_input_high_ROI.(high_met_nm),...
        sem_RT_f_input_high_ROI.(high_met_nm));
    RT_f_inputs_low_ROI_data_low_met_hdl.LineStyle = '--';
    RT_f_inputs_low_ROI_data_low_met_hdl.LineWidth = lWidth;
    RT_f_inputs_low_ROI_data_low_met_hdl.Color = purple;
    RT_f_inputs_high_ROI_data_low_met_hdl.LineStyle = '--';
    RT_f_inputs_high_ROI_data_low_met_hdl.LineWidth = lWidth;
    RT_f_inputs_high_ROI_data_low_met_hdl.Color = orange;
    RT_f_inputs_low_ROI_data_high_met_hdl.LineStyle = '-';
    RT_f_inputs_low_ROI_data_high_met_hdl.LineWidth = lWidth;
    RT_f_inputs_low_ROI_data_high_met_hdl.Color = purple;
    RT_f_inputs_high_ROI_data_high_met_hdl.LineStyle = '-';
    RT_f_inputs_high_ROI_data_high_met_hdl.LineWidth = lWidth;
    RT_f_inputs_high_ROI_data_high_met_hdl.Color = orange;
    legend([RT_f_inputs_low_ROI_data_low_met_hdl,...
        RT_f_inputs_high_ROI_data_low_met_hdl,...
        RT_f_inputs_low_ROI_data_high_met_hdl,...
        RT_f_inputs_high_ROI_data_high_met_hdl],...
        {['low ',ROI_short_nm,' - low ',metabolite_nm],...
        ['high ',ROI_short_nm,' - low ',metabolite_nm],...
        ['low ',ROI_short_nm,' - high ',metabolite_nm],...
        ['high ',ROI_short_nm,' - high ',metabolite_nm]});
    legend('boxoff');
    xlabel({input_prm_nm; full_bhv_taskName});
    ylabel({'RT (s)'; full_bhv_taskName});
    legend_size(pSize);
end

%% compare slopes: should do an ANOVA or GLM with dmPFC and GSH as factors and look at the interaction(s)

% [slope_low, slope_high,...
%     intercept_low, intercept_high] = deal(NaN(1,NS));
% for iS = 1:NS
%     b1 = glmfit(1:nBins, RT_f_input_low_ROI(:,iS),'normal');
%     intercept_low(iS) = b1(1);
%     slope_low(iS) = b1(2);
%     b2 = glmfit(1:nBins, RT_f_input_high_ROI(:,iS),'normal');
%     intercept_high(iS) = b2(1);
%     slope_high(iS) = b2(2);
% end
% [~,pval.intercept]=ttest(intercept_low,intercept_high);
% [~,pval.slope]=ttest(slope_low,slope_high);