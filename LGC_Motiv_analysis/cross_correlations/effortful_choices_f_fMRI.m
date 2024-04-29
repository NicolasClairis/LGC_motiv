function[beta, pval] = effortful_choices_f_fMRI(study_nm)
% [beta, pval] = effortful_choices_f_fMRI(study_nm)
% effortful_choices_f_fMRI checks whether fMRI predicts % of efforts done
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% OUTPUTS
% beta: structure with regression estimates
%
% pval: p.value for regression estimates

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% initialize parameters of interest
[hE_percentage.EpEm,...
    hE_percentage.Ep,...
    hE_percentage.Em,...
    avg_fMRI.EpEm,...
    avg_fMRI.Ep,...
    avg_fMRI.Em] = deal(NaN(1,NS));

%% ROI selection and extraction
GLM_str = inputdlg('Which fMRI GLM?'); % suggestion: GLM122 with only onset extraction
GLM = str2double(GLM_str);

% define fMRI ROI to use
[con_vec_all,...
    ~, ~, ~,...
    con_names,...
    ROI_coords] = ROI_extraction_group('study1', GLM,...
    subject_id, condition, 0);
n_ROIs = size(con_vec_all,3);
if n_ROIs > 1
    error(['more than 1 ROI selected, script cannot work that way',...
        'please focus on one and do it again.']);
end

% define regression estimate to look for in the fMRI GLM
EpEm_choice_onset_idx = strcmp(con_names,'Ep+Em ONSET choice RP E');
avg_fMRI.EpEm(:) = con_vec_all(EpEm_choice_onset_idx, :, 1);
Ep_choice_onset_idx = strcmp(con_names,'Ep ONSET choice RP E');
avg_fMRI.Ep(:) = con_vec_all(Ep_choice_onset_idx, :, 1);
Em_choice_onset_idx = strcmp(con_names,'Em ONSET choice RP E');
avg_fMRI.Em(:) = con_vec_all(Em_choice_onset_idx, :, 1);

% ask for fMRI short name for the figure
fMRI_short_nm = inputdlg('ROI short name for figure');

%% extract total proportion of effortful choices
[choice_hE_tmp] = choice_hE_proportion(study_nm, condition, [], 0);
hE_percentage.EpEm(:) = choice_hE_tmp.EpEm;
hE_percentage.Ep(:) = choice_hE_tmp.Ep;
hE_percentage.Em(:) = choice_hE_tmp.Em;

%% perform GLM for effort = f(fMRI)
% Physical + Mental
okSubs.EpEm = ~isnan(avg_fMRI.EpEm.*hE_percentage.EpEm);
[beta.EpEm,~,stats_EpEm] = glmfit(avg_fMRI.EpEm(okSubs.EpEm),...
    hE_percentage.EpEm(okSubs.EpEm),'normal');
pval.EpEm = stats_EpEm.p;
avg_fMRI_sorted.EpEm = sort(avg_fMRI.EpEm(okSubs.EpEm));
hE_percentage_fit.EpEm = glmval(beta.EpEm, avg_fMRI_sorted.EpEm, 'identity');

% Physical
okSubs.Ep = ~isnan(avg_fMRI.Ep.*hE_percentage.Ep);
[beta.Ep,~,stats_Ep] = glmfit(avg_fMRI.Ep(okSubs.Ep),...
    hE_percentage.Ep(okSubs.Ep),'normal');
pval.Ep = stats_Ep.p;
avg_fMRI_sorted.Ep = sort(avg_fMRI.Ep(okSubs.Ep));
hE_percentage_fit.Ep = glmval(beta.Ep, avg_fMRI_sorted.Ep, 'identity');

% Mental
okSubs.Em = ~isnan(avg_fMRI.Em.*hE_percentage.Em);
[beta.Em,~,stats_Em] = glmfit(avg_fMRI.EpEm(okSubs.Em),...
    hE_percentage.Em(okSubs.Em),'normal');
pval.Em = stats_Em.p;
avg_fMRI_sorted.Em = sort(avg_fMRI.Em(okSubs.Em));
hE_percentage_fit.Em = glmval(beta.Em, avg_fMRI_sorted.Em, 'identity');

%% figure display
[pSize, lWidth, col, mSize] = general_fig_prm;
xlim_vals = [-10 10];
fig;

% Physical + Mental
subplot(2,2,1:2); hold on;
scat_hdl = scatter(avg_fMRI.EpEm, hE_percentage.EpEm);
scat_hdl.MarkerEdgeColor = col.black;
scat_hdl.SizeData = mSize;
scat_hdl.LineWidth = lWidth;
plot_hdl = plot(avg_fMRI_sorted.EpEm, hE_percentage_fit.EpEm);
plot_hdl.LineWidth = lWidth;
plot_hdl.Color = col.grey;
plot_hdl.LineStyle = '--';
xlim(xlim_vals);
ylim([0 100]);
xlabel(fMRI_short_nm);
ylabel('Choices (%)');
legend_size(pSize);

% Physical
subplot(2,2,3); hold on;
scat_hdl = scatter(avg_fMRI.Ep, hE_percentage.Ep);
scat_hdl.MarkerEdgeColor = col.black;
scat_hdl.SizeData = mSize;
scat_hdl.LineWidth = lWidth;
plot_hdl = plot(avg_fMRI_sorted.Ep, hE_percentage_fit.Ep);
plot_hdl.LineWidth = lWidth;
plot_hdl.Color = col.grey;
plot_hdl.LineStyle = '--';
xlim(xlim_vals);
ylim([0 100]);
xlabel(fMRI_short_nm);
ylabel('Physical choices (%)');
legend_size(pSize);

% Mental
subplot(2,2,4); hold on;
scat_hdl = scatter(avg_fMRI.Em, hE_percentage.Em);
scat_hdl.MarkerEdgeColor = col.black;
scat_hdl.SizeData = mSize;
scat_hdl.LineWidth = lWidth;
plot_hdl = plot(avg_fMRI_sorted.Em, hE_percentage_fit.Em);
plot_hdl.LineWidth = lWidth;
plot_hdl.Color = col.grey;
plot_hdl.LineStyle = '--';
xlim(xlim_vals);
ylim([0 100]);
xlabel(fMRI_short_nm);
ylabel('Mental choices (%)');
legend_size(pSize);
end % function