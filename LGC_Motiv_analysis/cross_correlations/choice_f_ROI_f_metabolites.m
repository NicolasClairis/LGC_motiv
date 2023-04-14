% choice_f_ROI_f_metabolites aims at looking at whether the slope of
% choices=f(ROI) is correlating with the level of metabolites.
% In other words, do the level of metabolites predict how much the brain
% area will have an impact on choices?

%% subject selection
study_nm = 'study1';
condition = subject_condition;
subject_id = LGCM_subject_selection(study_nm, condition);

%% general parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
nBins = 6;
ROI_RT_orth = 1; % orthogonalize ROI to RT
n_E_levels = 3;

%% step 1: extract the slope
% variable organized as following:
% [b_choice_f_fMRI.Ep.perElevel, b_choice_f_fMRI.Em.perElevel, b_choice_f_fMRI.EpEmPool.perElevel] = deal(NaN(2, NS, n_E_levels));
% X1 = constant; X2= slope
figDisp = 0;
[b_choice_f_fMRI, pval_choice_f_fMRI, fMRI_bins,...
    choice_hE_bins, choice_hE_fit_bins,...
    fMRI_ROI_short_nm, timePeriod_nm] = choice_f_ROI(nBins,...
    study_nm, subject_id, condition,...
    ROI_RT_orth,...
    figDisp);
b_E3_intercept.Ep = b_choice_f_fMRI.Ep.perElevel(1,:,3);
b_E3_intercept.Em = b_choice_f_fMRI.Em.perElevel(1,:,3);
b_E3_intercept.EpEm = b_choice_f_fMRI.EpEmPool.perElevel(1,:,3);
b_E3_slope.Ep = b_choice_f_fMRI.Ep.perElevel(2,:,3);
b_E3_slope.Em = b_choice_f_fMRI.Em.perElevel(2,:,3);
b_E3_slope.EpEm = b_choice_f_fMRI.EpEmPool.perElevel(2,:,3);

%% step 2: extract metabolites
[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);
metabolite_ascOrder = sort(metabolite_allSubs);

%% step 3 test the correlation between the slope and the metabolites
% test for intercept
[b_E3_intercept_metabolite.Ep,~,stats_E3_incercept_metabolite.Ep] = glmfit(metabolite_allSubs, b_E3_intercept.Ep, 'normal');
[b_E3_intercept_metabolite.Em,~,stats_E3_incercept_metabolite.Em] = glmfit(metabolite_allSubs, b_E3_intercept.Em, 'normal');
[b_E3_intercept_metabolite.EpEm,~,stats_E3_incercept_metabolite.EpEm] = glmfit(metabolite_allSubs, b_E3_intercept.EpEm, 'normal');
% test for slope
[b_E3_slope_metabolite.Ep,~,stats_E3_slope_metabolite.Ep] = glmfit(metabolite_allSubs, b_E3_slope.Ep, 'normal');
[b_E3_slope_metabolite.Em,~,stats_E3_slope_metabolite.Em] = glmfit(metabolite_allSubs, b_E3_slope.Em, 'normal');
[b_E3_slope_metabolite.EpEm,~,stats_E3_slope_metabolite.EpEm] = glmfit(metabolite_allSubs, b_E3_slope.EpEm, 'normal');

% extract fit
fit_E3_interc_mb.Ep = glmval(b_E3_intercept_metabolite.Ep,metabolite_ascOrder,'identity');
fit_E3_interc_mb.Em = glmval(b_E3_intercept_metabolite.Em,metabolite_ascOrder,'identity');
fit_E3_interc_mb.EpEm = glmval(b_E3_intercept_metabolite.EpEm,metabolite_ascOrder,'identity');
fit_E3_slope_mb.Ep = glmval(b_E3_slope_metabolite.Ep,metabolite_ascOrder,'identity');
fit_E3_slope_mb.Em = glmval(b_E3_slope_metabolite.Em,metabolite_ascOrder,'identity');
fit_E3_slope_mb.EpEm = glmval(b_E3_slope_metabolite.EpEm,metabolite_ascOrder,'identity');

% extract p.value
pval.intercept.Ep = stats_E3_incercept_metabolite.Ep.p;
pval.intercept.Em = stats_E3_incercept_metabolite.Em.p;
pval.intercept.EpEm = stats_E3_incercept_metabolite.EpEm.p;
pval.slope.Ep = stats_E3_slope_metabolite.Ep.p;
pval.slope.Em = stats_E3_slope_metabolite.Em.p;
pval.slope.EpEm = stats_E3_slope_metabolite.EpEm.p;

% extract R²
% intercept
stats_Ep_intercept = fitlm(metabolite_allSubs, b_E3_intercept.Ep);
stats_Em_intercept = fitlm(metabolite_allSubs, b_E3_intercept.Em);
stats_EpEm_intercept = fitlm(metabolite_allSubs, b_E3_intercept.EpEm);
R2.intercept.Ep = stats_Ep_intercept.Rsquared.Adjusted;
R2.intercept.Em = stats_Em_intercept.Rsquared.Adjusted;
R2.intercept.EpEm = stats_EpEm_intercept.Rsquared.Adjusted;
% slope
stats_Ep_slope = fitlm(metabolite_allSubs, b_E3_slope.Ep);
stats_Em_slope = fitlm(metabolite_allSubs, b_E3_slope.Em);
stats_EpEm_slope = fitlm(metabolite_allSubs, b_E3_slope.EpEm);
R2.slope.Ep = stats_Ep_slope.Rsquared.Adjusted;
R2.slope.Em = stats_Em_slope.Rsquared.Adjusted;
R2.slope.EpEm = stats_EpEm_slope.Rsquared.Adjusted;

%% figures
pSize = 30;
lWidth = 3;
grey = [143 143 143]./255;

%% intercept
fig;
% Ep
subplot(1,2,1); hold on;
scat_hdl = scatter(metabolite_allSubs, b_E3_intercept.Ep);
fit_hdl = plot(metabolite_ascOrder, fit_E3_interc_mb.Ep);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '-';
fit_hdl.Color = grey;
xlabel([MRS_ROI_nm,' - ',metabolite_nm]);
ylabel(['Ep E3 choice=f(',fMRI_ROI_short_nm,') intercept']);
legend_size(pSize);

% Em
subplot(1,2,2); hold on;
scat_hdl = scatter(metabolite_allSubs, b_E3_intercept.Em);
fit_hdl = plot(metabolite_ascOrder, fit_E3_interc_mb.Em);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '-';
fit_hdl.Color = grey;
xlabel([MRS_ROI_nm,' - ',metabolite_nm]);
ylabel(['Em E3 choice=f(',fMRI_ROI_short_nm,') intercept']);
legend_size(pSize);


%% slope
fig;
% Ep
subplot(1,2,1); hold on;
scat_hdl = scatter(metabolite_allSubs, b_E3_slope.Ep);
fit_hdl = plot(metabolite_ascOrder, fit_E3_slope_mb.Ep);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '-';
fit_hdl.Color = grey;
xlabel([MRS_ROI_nm,' - ',metabolite_nm]);
ylabel(['Ep E3 choice=f(',fMRI_ROI_short_nm,') slope']);
legend_size(pSize);

% Em
subplot(1,2,2); hold on;
scat_hdl = scatter(metabolite_allSubs, b_E3_slope.Em);
fit_hdl = plot(metabolite_ascOrder, fit_E3_slope_mb.Em);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '-';
fit_hdl.Color = grey;
xlabel([MRS_ROI_nm,' - ',metabolite_nm]);
ylabel(['Em E3 choice=f(',fMRI_ROI_short_nm,') slope']);
legend_size(pSize);

%% pool Ep+Em
fig;
% intercept
subplot(1,2,1); hold on;
scat_hdl = scatter(metabolite_allSubs, b_E3_intercept.EpEm);
fit_hdl = plot(metabolite_ascOrder, fit_E3_interc_mb.EpEm);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '-';
fit_hdl.Color = grey;
xlabel([MRS_ROI_nm,' - ',metabolite_nm]);
ylabel(['Ep+Em E3 choice=f(',fMRI_ROI_short_nm,') intercept']);
legend_size(pSize);
% slope
subplot(1,2,2); hold on;
scat_hdl = scatter(metabolite_allSubs, b_E3_slope.EpEm);
fit_hdl = plot(metabolite_ascOrder, fit_E3_slope_mb.EpEm);
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '-';
fit_hdl.Color = grey;
xlabel([MRS_ROI_nm,' - ',metabolite_nm]);
ylabel(['Ep+Em E3 choice=f(',fMRI_ROI_short_nm,') slope']);
legend_size(pSize);