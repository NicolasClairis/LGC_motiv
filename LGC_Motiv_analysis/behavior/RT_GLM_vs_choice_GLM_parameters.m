% function[] = RT_GLM_vs_choice_GLM_parameters()


%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% study names
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% subject selection
if ~exist('subject_id','var') || isempty(subject_id)
    condition = subject_condition;
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end

%% extract RT parameters
[betas_RT, pval_RT, betas_grp_RT, pval_grp_RT] = RT_GLM(0, computerRoot, study_nm, subject_id);

%% extract behavioral parameters
[prm_choices, mdlType, mdlN] = prm_extraction(study_nm, subject_id);

%% test correlation between RT coefficient on confidence and choice coefficient on effort(s)
bUncertainty_RT = -betas_RT.conf;
bDeltaEp_RT = betas_RT.deltaEp;
bDeltaEm_RT = betas_RT.deltaEm;
switch mdlType
    case 'bayesian'
        kEp_choice = prm_choices.kEp;
        kEm_choice = prm_choices.kEm;
    case 'simple' % be careful in that case, kE are reverted (higher kE = less sensitive to efforts => needs to be flipped)
        kEp_choice = -prm_choices.kEp;
        kEm_choice = -prm_choices.kEm;
end

% extract subjects without NaN
okSubs = ~isnan(bUncertainty_RT.*kEp_choice);

%% test for bRT-uncertainty
% test physical effort
[betas.RT_bUnc_kEp, ~, stats.RT_bUnc_kEp] = glmfit(bUncertainty_RT(okSubs), kEp_choice(okSubs),'normal');
pval.RT_bUnc_kEp = stats.RT_bUnc_kEp.p;
kEpFit_f_bUncertainty = glmval(betas.RT_bUnc_kEp, bUncertainty_RT, 'identity');
[bUncertainty_RT_asc, idx_bU_RT_asc] = sort(bUncertainty_RT,'ascend');
kEpFit_f_bUncertainty = kEpFit_f_bUncertainty(idx_bU_RT_asc);

% test mental effort
[betas.RT_bUnc_kEm, ~, stats.RT_bUnc_kEm] = glmfit(bUncertainty_RT(okSubs), kEm_choice(okSubs),'normal');
pval.RT_bUnc_kEm = stats.RT_bUnc_kEm.p;
kEmFit_f_bUncertainty = glmval(betas.RT_bUnc_kEm, bUncertainty_RT, 'identity');
kEmFit_f_bUncertainty = kEmFit_f_bUncertainty(idx_bU_RT_asc);

%% test for bRT-kE
% test physical effort
[betas.RT_bDeltaEp_kEp, ~, stats.RT_bDeltaEp_kEp] = glmfit(bDeltaEp_RT(okSubs), kEp_choice(okSubs),'normal');
pval.RT_bDeltaEp_kEp = stats.RT_bDeltaEp_kEp.p;
kEpFit_f_bDeltaEp = glmval(betas.RT_bDeltaEp_kEp, bDeltaEp_RT, 'identity');
[bDeltaEp_RT_asc, idx_bDeltaEp_RT_asc] = sort(bDeltaEp_RT,'ascend');
kEpFit_f_bDeltaEp = kEpFit_f_bDeltaEp(idx_bDeltaEp_RT_asc);

% test mental effort
[betas.RT_bDeltaEm_kEm, ~, stats.RT_bDeltaEm_kEm] = glmfit(bDeltaEm_RT(okSubs), kEm_choice(okSubs),'normal');
pval.RT_bDeltaEm_kEm = stats.RT_bDeltaEm_kEm.p;
kEmFit_f_bDeltaEm = glmval(betas.RT_bDeltaEm_kEm, bDeltaEm_RT, 'identity');
[bDeltaEm_RT_asc, idx_bDeltaEm_RT_asc] = sort(bDeltaEm_RT,'ascend');
kEmFit_f_bDeltaEm = kEmFit_f_bDeltaEm(idx_bDeltaEm_RT_asc);

%% display figure
pSize = 50;
mSize = 50;
lWidth = 1;
lWidthPlot = 3;
greyCol = [143 143 143]./255;

% uncertainty
fig;
subplot(1,2,1); hold on;
scatter_hdl = scatter(bUncertainty_RT, kEp_choice);
plotHdl = plot(bUncertainty_RT_asc, kEpFit_f_bUncertainty, 'LineStyle','--');
scatter_hdl.MarkerEdgeColor = 'k';
scatter_hdl.SizeData = mSize;
scatter_hdl.LineWidth = lWidth;
plotHdl.LineWidth = lWidthPlot;
plotHdl.Color = greyCol;
xlabel('RT b(uncertainty)');
ylabel('choice kEp');
legend_size(pSize);

% mental effort
subplot(1,2,2); hold on;
scatter_hdl = scatter(bUncertainty_RT, kEm_choice);
plotHdl = plot(bUncertainty_RT_asc, kEmFit_f_bUncertainty, 'LineStyle','--');
scatter_hdl.MarkerEdgeColor = 'k';
scatter_hdl.SizeData = mSize;
scatter_hdl.LineWidth = lWidth;
plotHdl.LineWidth = lWidthPlot;
plotHdl.Color = greyCol;
xlabel('RT b(uncertainty)');
ylabel('choice kEm');
legend_size(pSize);

% deltaE and RT
fig;
subplot(1,2,1); hold on;
scatter_hdl = scatter(bDeltaEp_RT, kEp_choice);
plotHdl = plot(bDeltaEp_RT_asc, kEpFit_f_bDeltaEp, 'LineStyle','--');
scatter_hdl.MarkerEdgeColor = 'k';
scatter_hdl.SizeData = mSize;
scatter_hdl.LineWidth = lWidth;
plotHdl.LineWidth = lWidthPlot;
plotHdl.Color = greyCol;
xlabel('RT b(delta Ep)');
ylabel('choice kEp');
legend_size(pSize);

% mental effort
subplot(1,2,2); hold on;
scatter_hdl = scatter(bDeltaEm_RT, kEm_choice);
plotHdl = plot(bDeltaEm_RT_asc, kEmFit_f_bDeltaEm, 'LineStyle','--');
scatter_hdl.MarkerEdgeColor = 'k';
scatter_hdl.SizeData = mSize;
scatter_hdl.LineWidth = lWidth;
plotHdl.LineWidth = lWidthPlot;
plotHdl.Color = greyCol;
xlabel('RT b(delta Em)');
ylabel('choice kEm');
legend_size(pSize);
% end % function