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
[b_choice_f_fMRI, pval, fMRI_bins,...
    choice_hE_bins, choice_hE_fit_bins,...
    fMRI_ROI_short_nm, timePeriod_nm] = choice_f_ROI(nBins,...
    study_nm, subject_id, condition,...
    ROI_RT_orth,...
    figDisp);
b_E3_slope.Ep = b_choice_f_fMRI.Ep.perElevel(2,:,3);
b_E3_slope.Em = b_choice_f_fMRI.Em.perElevel(2,:,3);
b_E3_slope.EpEm = b_choice_f_fMRI.EpEmPool.perElevel(2,:,3);

%% step 2: extract metabolites
[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);
metabolite_ascOrder = sort(metabolite_allSubs);

%% step 3 test the correlation between the slope and the metabolites
[b_E3_slope_metabolite.Ep,~,stats_E3_slope_metabolite.Ep] = glmfit(metabolite_allSubs, b_E3_slope.Ep, 'normal');
[b_E3_slope_metabolite.Em,~,stats_E3_slope_metabolite.Em] = glmfit(metabolite_allSubs, b_E3_slope.Em, 'normal');
[b_E3_slope_metabolite.EpEm,~,stats_E3_slope_metabolite.EpEm] = glmfit(metabolite_allSubs, b_E3_slope.EpEm, 'normal');

% extract fit
fit_E3_mb.Ep = glmval(b_E3_slope_metabolite.Ep,metabolite_ascOrder,'identity');
fit_E3_mb.Em = glmval(b_E3_slope_metabolite.Em,metabolite_ascOrder,'identity');
fit_E3_mb.EpEm = glmval(b_E3_slope_metabolite.EpEm,metabolite_ascOrder,'identity');

% extract p.value
pval.Ep = stats_E3_slope_metabolite.Ep.p;
pval.Em = stats_E3_slope_metabolite.Em.p;
pval.EpEm = stats_E3_slope_metabolite.EpEm.p;

%% figure
pSize = 30;
lWidth = 3;
grey_ = [143 143 143]./255;

fig;
% Ep
subplot(1,2,1); hold on;
scat_hdl = scatter(metabolite_allSubs, b_E3_slope.Ep);
fit_hdl = plot(metabolite_ascOrder, fit_E3_mb.Ep);
scat_hdl.LineWidth = lWidth;
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '-';
fit_hdl.Color = grey;
xlabel([MRS_ROI_nm,' - ',metabolite_nm]);
ylabel(['Ep E3 choice=f(',fMRI_ROI_short_nm,') slope']);
legend_size(pSize);

% Em
subplot(1,2,2); hold on;
scat_hdl = scatter(metabolite_allSubs, b_E3_slope.Em);
fit_hdl = plot(metabolite_ascOrder, fit_E3_mb.Em);
scat_hdl.LineWidth = lWidth;
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '-';
fit_hdl.Color = grey;
xlabel([MRS_ROI_nm,' - ',metabolite_nm]);
ylabel(['Em E3 choice=f(',fMRI_ROI_short_nm,') slope']);
legend_size(pSize);

%% pool
fig;
scat_hdl = scatter(metabolite_allSubs, b_E3_slope.EpEm);
fit_hdl = plot(metabolite_ascOrder, fit_E3_mb.EpEm);
scat_hdl.LineWidth = lWidth;
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '-';
fit_hdl.Color = grey;
xlabel([MRS_ROI_nm,' - ',metabolite_nm]);
ylabel(['Ep+Em E3 choice=f(',fMRI_ROI_short_nm,') slope']);
legend_size(pSize);