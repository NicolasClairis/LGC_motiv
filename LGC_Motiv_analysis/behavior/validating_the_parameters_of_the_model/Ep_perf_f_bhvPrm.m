% performance in physical task (latency, overshoot, peak force, etc.)
% depending on behavioral parameters (kEp, kFp) for the selected subjects.

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
computerRoot = LGCM_root_paths;
dataRoot = [computerRoot, filesep, study_nm, filesep];

%% load behavioral parameters
[mdlType, mdlN] = behavioral_model_selection;
prm = prm_extraction(study_nm, subject_id, mdlType, mdlN);
% extract parameter of interest
prmToTest = {'kFp','kEp'};
nPrmToTest = length(prmToTest);

%% load force variables
[latency.avg.allTrials,...
    AUC.avg.allTrials, forcePeak.avg.allTrials, AUC_overshoot.avg.allTrials,...
    AUC_N.avg.allTrials, forcePeak_N.avg.allTrials, AUC_overshoot_N.avg.allTrials] = deal(NaN(1,NS));
n_hE_levels = 3;
nBins_Ech = 2;
[latency.avg.choice_hE, latency.avg.choice_lE,...
    AUC.avg.choice_hE, AUC.avg.choice_lE,...
    forcePeak.avg.choice_hE, forcePeak.avg.choice_lE,...
    AUC_overshoot.avg.choice_hE, AUC_overshoot.avg.choice_lE,...
    AUC_N.avg.choice_hE, AUC_N.avg.choice_lE,...
    forcePeak_N.avg.choice_hE, forcePeak_N.avg.choice_lE,...
    AUC_overshoot_N.avg.choice_hE, AUC_overshoot_N.avg.choice_lE] = deal(NaN(n_hE_levels,NS));
[latency.intercept.choice_hE, latency.intercept.choice_lE,...
    AUC.intercept.choice_hE, AUC.intercept.choice_lE,...
    forcePeak.intercept.choice_hE, forcePeak.intercept.choice_lE,...
    AUC_overshoot.intercept.choice_hE, AUC_overshoot.intercept.choice_lE,...
    AUC_N.intercept.choice_hE, AUC_N.intercept.choice_lE,...
    forcePeak_N.intercept.choice_hE, forcePeak_N.intercept.choice_lE,...
    AUC_overshoot_N.intercept.choice_hE, AUC_overshoot_N.intercept.choice_lE,...
    latency.slope.choice_hE, latency.slope.choice_lE,...
    AUC.slope.choice_hE, AUC.slope.choice_lE,...
    forcePeak.slope.choice_hE, forcePeak.slope.choice_lE,...
    AUC_overshoot.slope.choice_hE, AUC_overshoot.slope.choice_lE,...
    AUC_N.slope.choice_hE, AUC_N.slope.choice_lE,...
    forcePeak_N.slope.choice_hE, forcePeak_N.slope.choice_lE,...
    AUC_overshoot_N.slope.choice_hE, AUC_overshoot_N.slope.choice_lE] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    [latency_tmp,...
        AUC_tmp, forcePeak_tmp, AUC_overshoot_tmp,...
        AUC_N_tmp, forcePeak_N_tmp, AUC_overshoot_N_tmp] = avg_Ep_perf_perSub(study_nm, sub_nm, condition);
    % extract average across all trials
    latency.avg.allTrials(iS)           = latency_tmp.allRuns.allTrials;
    AUC.avg.allTrials(iS)               = AUC_tmp.allRuns.allTrials;
    forcePeak.avg.allTrials(iS)         = forcePeak_tmp.allRuns.allTrials;
    AUC_overshoot.avg.allTrials(iS)     = AUC_overshoot_tmp.allRuns.allTrials;
    AUC_N.avg.allTrials(iS)             = AUC_N_tmp.allRuns.allTrials;
    forcePeak_N.avg.allTrials(iS)       = forcePeak_N_tmp.allRuns.allTrials;
    AUC_overshoot_N.avg.allTrials(iS)   = AUC_overshoot_N_tmp.allRuns.allTrials;
    % extract average per (effort level*choice)
    latency.avg.choice_hE(:,iS)         = latency_tmp.allRuns.per_hE.choice_highE;
    latency.avg.choice_lE(:,iS)         = latency_tmp.allRuns.per_hE.choice_lowE;
    AUC.avg.choice_hE(:,iS)             = AUC_tmp.allRuns.per_hE.choice_highE;
    AUC.avg.choice_lE(:,iS)             = AUC_tmp.allRuns.per_hE.choice_lowE;
    forcePeak.avg.choice_hE(:,iS)       = forcePeak_tmp.allRuns.per_hE.choice_highE;
    forcePeak.avg.choice_lE(:,iS)       = forcePeak_tmp.allRuns.per_hE.choice_lowE;
    AUC_overshoot.avg.choice_hE(:,iS)   = AUC_overshoot_tmp.allRuns.per_hE.choice_highE;
    AUC_overshoot.avg.choice_lE(:,iS)   = AUC_overshoot_tmp.allRuns.per_hE.choice_lowE;
    AUC_N.avg.choice_hE(:,iS)           = AUC_N_tmp.allRuns.per_hE.choice_highE;
    AUC_N.avg.choice_lE(:,iS)           = AUC_N_tmp.allRuns.per_hE.choice_lowE;
    forcePeak_N.avg.choice_hE(:,iS)     = forcePeak_N_tmp.allRuns.per_hE.choice_highE;
    forcePeak_N.avg.choice_lE(:,iS)     = forcePeak_N_tmp.allRuns.per_hE.choice_lowE;
    AUC_overshoot_N.avg.choice_hE(:,iS) = AUC_overshoot_N_tmp.allRuns.per_hE.choice_highE;
    AUC_overshoot_N.avg.choice_lE(:,iS) = AUC_overshoot_N_tmp.allRuns.per_hE.choice_lowE;
    % intercept
    latency.intercept.choice_hE(iS)         = latency_tmp.intercept.choice_highE;
    latency.intercept.choice_lE(iS)         = latency_tmp.intercept.choice_lowE;
    AUC.intercept.choice_hE(iS)             = AUC_tmp.intercept.choice_highE;
    AUC.intercept.choice_lE(iS)             = AUC_tmp.intercept.choice_lowE;
    forcePeak.intercept.choice_hE(iS)       = forcePeak_tmp.intercept.choice_highE;
    forcePeak.intercept.choice_lE(iS)       = forcePeak_tmp.intercept.choice_lowE;
    AUC_overshoot.intercept.choice_hE(iS)   = AUC_overshoot_tmp.intercept.choice_highE;
    AUC_overshoot.intercept.choice_lE(iS)   = AUC_overshoot_tmp.intercept.choice_lowE;
    AUC_N.intercept.choice_hE(iS)           = AUC_N_tmp.intercept.choice_highE;
    AUC_N.intercept.choice_lE(iS)           = AUC_N_tmp.intercept.choice_lowE;
    forcePeak_N.intercept.choice_hE(iS)     = forcePeak_N_tmp.intercept.choice_highE;
    forcePeak_N.intercept.choice_lE(iS)     = forcePeak_N_tmp.intercept.choice_lowE;
    AUC_overshoot_N.intercept.choice_hE(iS) = AUC_overshoot_N_tmp.intercept.choice_highE;
    AUC_overshoot_N.intercept.choice_lE(iS) = AUC_overshoot_N_tmp.intercept.choice_lowE;
    % slope
    latency.slope.choice_hE(iS)         = latency_tmp.slope.choice_highE;
    latency.slope.choice_lE(iS)         = latency_tmp.slope.choice_lowE;
    AUC.slope.choice_hE(iS)             = AUC_tmp.slope.choice_highE;
    AUC.slope.choice_lE(iS)             = AUC_tmp.slope.choice_lowE;
    forcePeak.slope.choice_hE(iS)       = forcePeak_tmp.slope.choice_highE;
    forcePeak.slope.choice_lE(iS)       = forcePeak_tmp.slope.choice_lowE;
    AUC_overshoot.slope.choice_hE(iS)   = AUC_overshoot_tmp.slope.choice_highE;
    AUC_overshoot.slope.choice_lE(iS)   = AUC_overshoot_tmp.slope.choice_lowE;
    AUC_N.slope.choice_hE(iS)           = AUC_N_tmp.slope.choice_highE;
    AUC_N.slope.choice_lE(iS)           = AUC_N_tmp.slope.choice_lowE;
    forcePeak_N.slope.choice_hE(iS)     = forcePeak_N_tmp.slope.choice_highE;
    forcePeak_N.slope.choice_lE(iS)     = forcePeak_N_tmp.slope.choice_lowE;
    AUC_overshoot_N.slope.choice_hE(iS) = AUC_overshoot_N_tmp.slope.choice_highE;
    AUC_overshoot_N.slope.choice_lE(iS) = AUC_overshoot_N_tmp.slope.choice_lowE;
end % subject loop

%% figure
pSize = 30;
lWidth = 3;
black = [0 0 0];
white = [1 1 1];
grey = [143 143 143]./255;
lowBinCol = [160 52 114]./255;
mediumBinCol = [59 131 189]./255;
highBinCol = [14 41 75]./255;

for iPrm = 1:nPrmToTest
    prm_nm = prmToTest{iPrm};
    prm_of_interest = prm.(prm_nm);
    goodSubs_prm = ~isnan(prm_of_interest);
    
    %% statistical tests
    % parameter vs all force variables
    goodSubs_latency = goodSubs_prm.*(~isnan(latency.avg.allTrials)) == 1;
    [b_latency,~,stats_latency] = glmfit(prm_of_interest(goodSubs_latency),...
        latency.avg.allTrials(goodSubs_latency),'normal');
    prm_ascOrder_latency = sort(prm_of_interest(goodSubs_latency));
    latency_fit = glmval(b_latency, prm_ascOrder_latency,'identity');
    betas.(prm_nm).latency = b_latency;
    pval.(prm_nm).latency = stats_latency.p;
    
    goodSubs_AUC = goodSubs_prm.*(~isnan(AUC.avg.allTrials)) == 1;
    [b_AUC,~,stats_AUC] = glmfit(prm_of_interest(goodSubs_AUC),...
        AUC.avg.allTrials(goodSubs_AUC),'normal');
    prm_ascOrder_AUC = sort(prm_of_interest(goodSubs_AUC));
    AUC_fit = glmval(b_AUC, prm_ascOrder_AUC,'identity');
    betas.(prm_nm).AUC = b_AUC;
    pval.(prm_nm).AUC = stats_AUC.p;
    
    goodSubs_forcePeak = goodSubs_prm.*(~isnan(forcePeak.avg.allTrials)) == 1;
    [b_forcePeak,~,stats_forcePeak] = glmfit(prm_of_interest(goodSubs_forcePeak),...
        forcePeak.avg.allTrials(goodSubs_forcePeak),'normal');
    prm_ascOrder_forcePeak = sort(prm_of_interest(goodSubs_forcePeak));
    forcePeak_fit = glmval(b_forcePeak, prm_ascOrder_forcePeak,'identity');
    betas.(prm_nm).forcePeak = b_forcePeak;
    pval.(prm_nm).forcePeak = stats_forcePeak.p;
    
    goodSubs_AUC_overshoot = goodSubs_prm.*(~isnan(AUC_overshoot.avg.allTrials)) == 1;
    [b_AUC_overshoot,~,stats_AUC_overshoot] = glmfit(prm_of_interest(goodSubs_AUC_overshoot),...
        AUC_overshoot.avg.allTrials(goodSubs_AUC_overshoot),'normal');
    prm_ascOrder_AUC_overshoot = sort(prm_of_interest(goodSubs_AUC_overshoot));
    AUC_overshoot_fit = glmval(b_AUC_overshoot, prm_ascOrder_AUC_overshoot,'identity');
    betas.(prm_nm).AUC_overshoot = b_AUC_overshoot;
    pval.(prm_nm).AUC_overshoot = stats_AUC_overshoot.p;
    
    goodSubs_AUC_N = goodSubs_prm.*(~isnan(AUC_N.avg.allTrials)) == 1;
    [b_AUC_N,~,stats_AUC_N] = glmfit(prm_of_interest(goodSubs_AUC_N),...
        AUC_N.avg.allTrials(goodSubs_AUC_N),'normal');
    prm_ascOrder_AUC_N = sort(prm_of_interest(goodSubs_AUC_N));
    AUC_N_fit = glmval(b_AUC_N, prm_ascOrder_AUC_N,'identity');
    betas.(prm_nm).AUC_N = b_AUC_N;
    pval.(prm_nm).AUC_N = stats_AUC_N.p;
    
    goodSubs_forcePeak_N = goodSubs_prm.*(~isnan(forcePeak_N.avg.allTrials)) == 1;
    [b_forcePeak_N,~,stats_forcePeak_N] = glmfit(prm_of_interest(goodSubs_forcePeak_N),...
        forcePeak_N.avg.allTrials(goodSubs_forcePeak_N),'normal');
    prm_ascOrder_forcePeak_N = sort(prm_of_interest(goodSubs_forcePeak_N));
    forcePeak_N_fit = glmval(b_forcePeak_N, prm_ascOrder_forcePeak_N,'identity');
    betas.(prm_nm).forcePeak_N = b_forcePeak_N;
    pval.(prm_nm).forcePeak_N = stats_forcePeak_N.p;
    
    goodSubs_AUC_overshoot_N = goodSubs_prm.*(~isnan(AUC_overshoot_N.avg.allTrials)) == 1;
    [b_AUC_overshoot_N,~,stats_AUC_overshoot_N] = glmfit(prm_of_interest(goodSubs_AUC_overshoot_N),...
        AUC_overshoot_N.avg.allTrials(goodSubs_AUC_overshoot_N),'normal');
    prm_ascOrder_AUC_overshoot = sort(prm_of_interest(goodSubs_AUC_overshoot_N));
    AUC_overshoot_N_fit = glmval(b_AUC_overshoot_N, prm_ascOrder_AUC_overshoot,'identity');
    betas.(prm_nm).AUC_overshoot_N = b_AUC_overshoot_N;
    pval.(prm_nm).AUC_overshoot_N = stats_AUC_overshoot_N.p;
    
    % parameter vs intercept and slope of effort curve
    interc_slope = {'intercept','slope'};
    choices = {'choice_lE','choice_hE'};
    for iC = 1:length(choices)
        ch_nm = choices{iC};
        [latency_bin.avg.(ch_nm),...
            AUC_bin.avg.(ch_nm),...
            forcePeak_bin.avg.(ch_nm),...
            AUC_overshoot_bin.avg.(ch_nm),...
            AUC_N_bin.avg.(ch_nm),...
            forcePeak_N_bin.avg.(ch_nm),...
            AUC_overshoot_N_bin.avg.(ch_nm)] = deal(NaN(n_hE_levels, nBins_Ech));
        
        for iIS = 1:length(interc_slope)
            IS_nm = interc_slope{iIS};
            
            goodSubs_latency_Ech = goodSubs_prm.*(~isnan(latency.(IS_nm).(ch_nm))) == 1;
            [b_latency_Ech,~,stats_latency_Ech] = glmfit(prm_of_interest(goodSubs_latency_Ech),...
                latency.(IS_nm).(ch_nm)(goodSubs_latency_Ech),'normal');
            betas.(prm_nm).latency_ch.(IS_nm).(ch_nm) = b_latency_Ech;
            pval.(prm_nm).latency_ch.(IS_nm).(ch_nm) = stats_latency_Ech.p;
            % extract bins for fit
            latency_bin.(IS_nm).(ch_nm) = do_bin2(latency.(IS_nm).(ch_nm), prm_of_interest, nBins_Ech, 0);
        
            goodSubs_AUC_Ech = goodSubs_prm.*(~isnan(AUC.(IS_nm).(ch_nm))) == 1;
            [b_AUC_Ech,~,stats_AUC_Ech] = glmfit(prm_of_interest(goodSubs_AUC_Ech),...
                AUC.(IS_nm).(ch_nm)(goodSubs_AUC_Ech),'normal');
            betas.(prm_nm).AUC_ch.(IS_nm).(ch_nm) = b_AUC_Ech;
            pval.(prm_nm).AUC_ch.(IS_nm).(ch_nm) = stats_AUC_Ech.p;
            % extract bins for fit
            AUC_bin.(IS_nm).(ch_nm) = do_bin2(AUC.(IS_nm).(ch_nm), prm_of_interest, nBins_Ech, 0);
        
            goodSubs_forcePeak_Ech = goodSubs_prm.*(~isnan(forcePeak.(IS_nm).(ch_nm))) == 1;
            [b_forcePeak_Ech,~,stats_forcePeak_Ech] = glmfit(prm_of_interest(goodSubs_forcePeak_Ech),...
                forcePeak.(IS_nm).(ch_nm)(goodSubs_forcePeak_Ech),'normal');
            betas.(prm_nm).forcePeak_ch.(IS_nm).(ch_nm) = b_forcePeak_Ech;
            pval.(prm_nm).forcePeak_ch.(IS_nm).(ch_nm) = stats_forcePeak_Ech.p;
            % extract bins for fit
            forcePeak_bin.(IS_nm).(ch_nm) = do_bin2(forcePeak.(IS_nm).(ch_nm), prm_of_interest, nBins_Ech, 0);
        
            goodSubs_AUC_overshoot_Ech = goodSubs_prm.*(~isnan(AUC_overshoot.(IS_nm).(ch_nm))) == 1;
            [b_AUC_overshoot_Ech,~,stats_AUC_overshoot_Ech] = glmfit(prm_of_interest(goodSubs_AUC_overshoot_Ech),...
                AUC_overshoot.(IS_nm).(ch_nm)(goodSubs_AUC_overshoot_Ech),'normal');
            betas.(prm_nm).AUC_overshoot_ch.(IS_nm).(ch_nm) = b_AUC_overshoot_Ech;
            pval.(prm_nm).AUC_overshoot_ch.(IS_nm).(ch_nm) = stats_AUC_overshoot_Ech.p;
            % extract bins for fit
            AUC_overshoot_bin.(IS_nm).(ch_nm) = do_bin2(AUC_overshoot.(IS_nm).(ch_nm), prm_of_interest, nBins_Ech, 0);
        
            goodSubs_AUC_N_Ech = goodSubs_prm.*(~isnan(AUC_N.(IS_nm).(ch_nm))) == 1;
            [b_AUC_N_Ech,~,stats_AUC_N_Ech] = glmfit(prm_of_interest(goodSubs_AUC_N_Ech),...
                AUC_N.(IS_nm).(ch_nm)(goodSubs_AUC_N_Ech),'normal');
            betas.(prm_nm).AUC_N_ch.(IS_nm).(ch_nm) = b_AUC_N_Ech;
            pval.(prm_nm).AUC_N_ch.(IS_nm).(ch_nm) = stats_AUC_N_Ech.p;
            % extract bins for fit
            AUC_N_bin.(IS_nm).(ch_nm) = do_bin2(AUC_N.(IS_nm).(ch_nm), prm_of_interest, nBins_Ech, 0);
        
            goodSubs_forcePeak_N_Ech = goodSubs_prm.*(~isnan(forcePeak_N.(IS_nm).(ch_nm))) == 1;
            [b_forcePeak_N_Ech,~,stats_forcePeak_N_Ech] = glmfit(prm_of_interest(goodSubs_forcePeak_N_Ech),...
                forcePeak_N.(IS_nm).(ch_nm)(goodSubs_forcePeak_N_Ech),'normal');
            betas.(prm_nm).forcePeak_N_ch.(IS_nm).(ch_nm) = b_forcePeak_N_Ech;
            pval.(prm_nm).forcePeak_N_ch.(IS_nm).(ch_nm) = stats_forcePeak_N_Ech.p;
            % extract bins for fit
            forcePeak_N_bin.(IS_nm).(ch_nm) = do_bin2(forcePeak_N.(IS_nm).(ch_nm), prm_of_interest, nBins_Ech, 0);
        
            goodSubs_AUC_overshoot_N_Ech = goodSubs_prm.*(~isnan(AUC_overshoot_N.(IS_nm).(ch_nm))) == 1;
            [b_AUC_overshoot_N_Ech,~,stats_AUC_overshoot_N_Ech] = glmfit(prm_of_interest(goodSubs_AUC_overshoot_N_Ech),...
                AUC_overshoot_N.(IS_nm).(ch_nm)(goodSubs_AUC_overshoot_N_Ech),'normal');
            betas.(prm_nm).AUC_overshoot_N_ch.(IS_nm).(ch_nm) = b_AUC_overshoot_N_Ech;
            pval.(prm_nm).AUC_overshoot_N_ch.(IS_nm).(ch_nm) = stats_AUC_overshoot_N_Ech.p;
            % extract bins for fit
            AUC_overshoot_N_bin.(IS_nm).(ch_nm) = do_bin2(AUC_overshoot_N.(IS_nm).(ch_nm), prm_of_interest, nBins_Ech, 0);
        
        end % intercept/slope
        
        % extract bins for raw values
        for iE = 1:n_hE_levels
            latency_bin.avg.(ch_nm)(iE,:)           = do_bin2(latency.avg.(ch_nm)(iE,:), prm_of_interest, nBins_Ech, 0);
            AUC_bin.avg.(ch_nm)(iE,:)               = do_bin2(AUC.avg.(ch_nm)(iE,:), prm_of_interest, nBins_Ech, 0);
            forcePeak_bin.avg.(ch_nm)(iE,:)         = do_bin2(forcePeak.avg.(ch_nm)(iE,:), prm_of_interest, nBins_Ech, 0);
            AUC_overshoot_bin.avg.(ch_nm)(iE,:)     = do_bin2(AUC_overshoot.avg.(ch_nm)(iE,:), prm_of_interest, nBins_Ech, 0);
            AUC_N_bin.avg.(ch_nm)(iE,:)             = do_bin2(AUC_N.avg.(ch_nm)(iE,:), prm_of_interest, nBins_Ech, 0);
            forcePeak_N_bin.avg.(ch_nm)(iE,:)       = do_bin2(forcePeak_N.avg.(ch_nm)(iE,:), prm_of_interest, nBins_Ech, 0);
            AUC_overshoot_N_bin.avg.(ch_nm)(iE,:)   = do_bin2(AUC_overshoot_N.avg.(ch_nm)(iE,:), prm_of_interest, nBins_Ech, 0);
        end
    end % choice loop
    
    %% figure with performance averaged
    fig;
    % latency = f(prm)
    subplot(3,3,1);
    scat_hdl = scatter(prm_of_interest(goodSubs_latency),...
        latency.avg.allTrials(goodSubs_latency));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_ascOrder_latency, latency_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('Latency (s)');
    legend_size(pSize);
    
    % AUC = f(prm)
    subplot(3,3,2);
    scat_hdl = scatter(prm_of_interest(goodSubs_AUC),...
        AUC.avg.allTrials(goodSubs_AUC));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_ascOrder_AUC, AUC_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('AUC (a.u.)');
    legend_size(pSize);
    
    % peak force (%) = f(prm)
    subplot(3,3,3);
    scat_hdl = scatter(prm_of_interest(goodSubs_forcePeak),...
        forcePeak.avg.allTrials(goodSubs_forcePeak));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_ascOrder_forcePeak, forcePeak_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('Peak force (%)');
    legend_size(pSize);
    
    % AUC overshoot = f(prm)
    subplot(3,3,4);
    scat_hdl = scatter(prm_of_interest(goodSubs_AUC_overshoot),...
        AUC_overshoot.avg.allTrials(goodSubs_AUC_overshoot));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_ascOrder_AUC_overshoot, AUC_overshoot_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('AUC overshoot (a.u.)');
    legend_size(pSize);
    
    % AUC (based in Newton force) = f(prm)
    subplot(3,3,5);
    scat_hdl = scatter(prm_of_interest(goodSubs_AUC_N),...
        AUC_N.avg.allTrials(goodSubs_AUC_N));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_ascOrder_AUC_N, AUC_N_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('AUC (a.u.)');
    legend_size(pSize);
    
    % peak force (N) = f(prm)
    subplot(3,3,6);
    scat_hdl = scatter(prm_of_interest(goodSubs_forcePeak_N),...
        forcePeak_N.avg.allTrials(goodSubs_forcePeak_N));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_ascOrder_forcePeak_N, forcePeak_N_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('Peak force (N)');
    legend_size(pSize);
    
    % AUC overshoot (based in force in N) = f(prm)
    subplot(3,3,7);
    scat_hdl = scatter(prm_of_interest(goodSubs_AUC_overshoot_N),...
        AUC_overshoot_N.avg.allTrials(goodSubs_AUC_overshoot_N));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_ascOrder_AUC_overshoot, AUC_overshoot_N_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('AUC overshoot (a.u.)');
    legend_size(pSize);
    
    %% figure with latency = f(E level)*choice depending on parameter
    fig;
    % latency = f(prm)
    for iBin = 1:nBins_Ech
        switch iBin
            case 1 % low parameter
                col = lowBinCol;
            case 2 % medium parameter
                col = mediumBinCol;
            case 3 % high parameter
                col = highBinCol;
            otherwise
                col = black;
        end
        scat_hdl_lE = scatter(1:n_hE_levels,...
            latency_bin.avg.choice_lE(:,iBin));
        scat_hdl_hE = scatter(1:n_hE_levels,...
            latency_bin.avg.choice_hE(:,iBin));
        scat_hdl_lE.LineWidth = lWidth;
        scat_hdl_lE.MarkerEdgeColor = black;
        scat_hdl_lE.MarkerFaceColor = grey;
        scat_hdl_hE.LineWidth = lWidth;
        scat_hdl_hE.MarkerEdgeColor = black;
        scat_hdl_hE.MarkerFaceColor = white;
        fit_hdl_lE = plot(1:n_hE_levels,...
            latency_bin.intercept.choice_lE(iBin) +...
            latency_bin.slope.choice_lE(iBin).*(1:n_hE_levels));
        fit_hdl_hE = plot(1:n_hE_levels,...
            latency_bin.intercept.choice_hE(iBin) +...
            latency_bin.slope.choice_hE(iBin).*(1:n_hE_levels));
        fit_hdl_lE.LineWidth = lWidth;
        fit_hdl_lE.Color = col;
        fit_hdl_lE.LineStyle = '--';
        fit_hdl_hE.LineWidth = lWidth;
        fit_hdl_hE.Color = col;
        fit_hdl_hE.LineStyle = '-';
    end
    xticks(1:n_hE_levels);
    xlabel('Effort level');
    ylabel('Latency (s)');
    legend_size(pSize);
    
    %% figure with peak force (%) = f(E level)*choice depending on parameter
    fig;
    % forcePeak = f(prm)
    for iBin = 1:nBins_Ech
        switch iBin
            case 1 % low parameter
                col = lowBinCol;
            case 2 % medium parameter
                col = mediumBinCol;
            case 3 % high parameter
                col = highBinCol;
            otherwise
                col = black;
        end
        scat_hdl_lE = scatter(1:n_hE_levels,...
            forcePeak_bin.avg.choice_lE(:,iBin));
        scat_hdl_hE = scatter(1:n_hE_levels,...
            forcePeak_bin.avg.choice_hE(:,iBin));
        scat_hdl_lE.LineWidth = lWidth;
        scat_hdl_lE.MarkerEdgeColor = black;
        scat_hdl_lE.MarkerFaceColor = grey;
        scat_hdl_hE.LineWidth = lWidth;
        scat_hdl_hE.MarkerEdgeColor = black;
        scat_hdl_hE.MarkerFaceColor = white;
        fit_hdl_lE = plot(1:n_hE_levels,...
            forcePeak_bin.intercept.choice_lE(iBin) +...
            forcePeak_bin.slope.choice_lE(iBin).*(1:n_hE_levels));
        fit_hdl_hE = plot(1:n_hE_levels,...
            forcePeak_bin.intercept.choice_hE(iBin) +...
            forcePeak_bin.slope.choice_hE(iBin).*(1:n_hE_levels));
        fit_hdl_lE.LineWidth = lWidth;
        fit_hdl_lE.Color = col;
        fit_hdl_lE.LineStyle = '--';
        fit_hdl_hE.LineWidth = lWidth;
        fit_hdl_hE.Color = col;
        fit_hdl_hE.LineStyle = '-';
    end
    xticks(1:n_hE_levels);
    xlabel('Effort level');
    ylabel('Peak force (%)');
    legend_size(pSize);
    
    %% figure with peak force (N) = f(E level)*choice depending on parameter
    fig;
    % forcePeak (N) = f(prm)
    for iBin = 1:nBins_Ech
        switch iBin
            case 1 % low parameter
                col = lowBinCol;
            case 2 % medium parameter
                col = mediumBinCol;
            case 3 % high parameter
                col = highBinCol;
            otherwise
                col = black;
        end
        scat_hdl_lE = scatter(1:n_hE_levels,...
            forcePeak_N_bin.avg.choice_lE(:,iBin));
        scat_hdl_hE = scatter(1:n_hE_levels,...
            forcePeak_N_bin.avg.choice_hE(:,iBin));
        scat_hdl_lE.LineWidth = lWidth;
        scat_hdl_lE.MarkerEdgeColor = black;
        scat_hdl_lE.MarkerFaceColor = grey;
        scat_hdl_hE.LineWidth = lWidth;
        scat_hdl_hE.MarkerEdgeColor = black;
        scat_hdl_hE.MarkerFaceColor = white;
        fit_hdl_lE = plot(1:n_hE_levels,...
            forcePeak_N_bin.intercept.choice_lE(iBin) +...
            forcePeak_N_bin.slope.choice_lE(iBin).*(1:n_hE_levels));
        fit_hdl_hE = plot(1:n_hE_levels,...
            forcePeak_N_bin.intercept.choice_hE(iBin) +...
            forcePeak_N_bin.slope.choice_hE(iBin).*(1:n_hE_levels));
        fit_hdl_lE.LineWidth = lWidth;
        fit_hdl_lE.Color = col;
        fit_hdl_lE.LineStyle = '--';
        fit_hdl_hE.LineWidth = lWidth;
        fit_hdl_hE.Color = col;
        fit_hdl_hE.LineStyle = '-';
    end
    xticks(1:n_hE_levels);
    xlabel('Effort level');
    ylabel('Peak force (N)');
    legend_size(pSize);
    
end % parameter loop