
resultsPath = fullfile('P:','boulot','postdoc_CarmenSandi',...
    'collab_with_CHUV_project');

%% load GSH
GSH_results_EP = readtable([resultsPath,filesep,'GSH_EP_forLGC.xlsx'],...
    'Sheet','EP');
NS_EP = size(GSH_results_EP,1);
GSH_results_ctrl = readtable([resultsPath,filesep,'GSH_EP_forLGC.xlsx'],...
    'Sheet','ctrl');
NS_ctrl = size(GSH_results_ctrl,1);
NS_all = NS_EP + NS_ctrl;

%% load motivational questionnaire scores
behaviorInfos = readtable([resultsPath,filesep,'DB_NClairis_30-05-2023.xlsx'],...
    'Sheet','Feuil1');

%% extract info with corresponding subjects
[GSH_EP, MADRS_EP] = deal(NaN(1,NS_EP));
subject_EP = cell(1,NS_EP);
[GSH_ctrl, MADRS_ctrl] = deal(NaN(1,NS_ctrl));
subject_ctrl = cell(1,NS_ctrl);
[GSH_all, MADRS_all, EP_state] = deal(NaN(1,NS_all));
subject_all = cell(1,NS_all);

for iS = 1:NS_EP
    % subject id
    subject_EP{iS} = strrep(GSH_results_EP.CodeLunep(iS),'LNAC','LNAC0');
    subject_all{iS} = strrep(GSH_results_EP.CodeLunep(iS),'LNAC','LNAC0');
    % extract metabolite concentration
    GSH_EP(iS) = GSH_results_EP.GSHCorrectedForTissueCompsition(iS);
    GSH_all(iS) = GSH_results_EP.GSHCorrectedForTissueCompsition(iS);
    % subject id for behavior
    bhv_idx_tmp = find(strcmp(behaviorInfos.CodeURS1,subject_EP{iS}));
    % look on the second column of subject codes if doesn't work in the
    % first one
    if isempty(bhv_idx_tmp)
        bhv_idx_tmp = find(strcmp(behaviorInfos.CodeURS2,subject_EP{iS}));
    end
    % motivational behavior info
    if ~isempty(bhv_idx_tmp)
        MADRS_EP(iS) = behaviorInfos.MADRSTOTAL(bhv_idx_tmp);
        MADRS_all(iS) = behaviorInfos.MADRSTOTAL(bhv_idx_tmp);
        EP_state(iS) = 1;
    end
end % subject loop EP

% same for controls
for iS = 1:NS_ctrl
    % subject id
    subject_ctrl{iS} = strrep(GSH_results_ctrl.CodeLunep(iS),'LNAC','LNAC0');
    subject_all{iS+NS_EP} = strrep(GSH_results_ctrl.CodeLunep(iS),'LNAC','LNAC0');
    % extract metabolite concentration
    GSH_ctrl(iS) = GSH_results_ctrl.GSHCorrectedForTissueCompsition(iS);
    GSH_all(iS+NS_EP) = GSH_results_ctrl.GSHCorrectedForTissueCompsition(iS);
    % subject id for behavior
    bhv_idx_tmp = find(strcmp(behaviorInfos.CodeURS1,subject_ctrl{iS}));
    % look on the second column of subject codes if doesn't work in the
    % first one
    if isempty(bhv_idx_tmp)
        bhv_idx_tmp = find(strcmp(behaviorInfos.CodeURS2,subject_ctrl{iS}));
    end
    % motivational behavior info
    if ~isempty(bhv_idx_tmp)
        MADRS_ctrl(iS) = behaviorInfos.MADRSTOTAL(bhv_idx_tmp);
        MADRS_all(iS+NS_EP) = behaviorInfos.MADRSTOTAL(bhv_idx_tmp);
        EP_state(iS) = 0;
    end
end % subject loop EP

%% basic correlation MADRS vs GSH
pSize = 30;
lWidth = 3;
mSize = 50;
grey = [174 174 174]./255;
MADRS_GSH_fig = fig;

% MADRS = f(GSH) early psychosis
ok_EP_subs = ~isnan(GSH_EP.*MADRS_EP);
[b_MADRS_GSH.EP,~,stats_MADRS_GSH.EP] = glmfit(GSH_EP(ok_EP_subs), MADRS_EP(ok_EP_subs),'normal');
GSH_EP_sorted = GSH_EP(ok_EP_subs);
MADRS_f_GSH_fit.EP = glmval(b_MADRS_GSH.EP, GSH_EP_sorted, 'identity');
subplot(1,3,1); hold on;
title('EP')
scat_hdl = scatter(GSH_EP(ok_EP_subs), MADRS_EP(ok_EP_subs));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
scat_hdl.SizeData = mSize;
plot_hdl = plot(GSH_EP_sorted, MADRS_f_GSH_fit.EP);
plot_hdl.LineStyle = '--';
plot_hdl.LineWidth = lWidth;
plot_hdl.Color = grey;
xlabel('GSH');
ylabel('MADRS score');
legend_size(pSize);

% MADRS = f(GSH) controls
ok_ctrl_subs = ~isnan(GSH_ctrl.*MADRS_ctrl);
[b_MADRS_GSH.ctrl,~,stats_MADRS_GSH.ctrl] = glmfit(GSH_ctrl(ok_ctrl_subs), MADRS_ctrl(ok_ctrl_subs),'normal');
GSH_ctrl_sorted = GSH_ctrl(ok_ctrl_subs);
MADRS_f_GSH_fit.ctrl = glmval(b_MADRS_GSH.ctrl, GSH_ctrl_sorted, 'identity');
subplot(1,3,2); hold on;
title('controls')
scat_hdl = scatter(GSH_ctrl(ok_ctrl_subs), MADRS_ctrl(ok_ctrl_subs));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
scat_hdl.SizeData = mSize;
plot_hdl = plot(GSH_ctrl_sorted, MADRS_f_GSH_fit.ctrl);
plot_hdl.LineStyle = '--';
plot_hdl.LineWidth = lWidth;
plot_hdl.Color = grey;
xlabel('GSH');
ylabel('MADRS score');
legend_size(pSize);

% MADRS = f(GSH) all
ok_all_subs = ~isnan(GSH_all.*MADRS_all);
[b_MADRS_GSH.all,~,stats_MADRS_GSH.all] = glmfit(GSH_all(ok_all_subs), MADRS_all(ok_all_subs),'normal');
GSH_all_sorted = GSH_all(ok_all_subs);
MADRS_f_GSH_fit.all = glmval(b_MADRS_GSH.all, GSH_all_sorted, 'identity');
subplot(1,3,3); hold on;
title('all')
scat_hdl = scatter(GSH_all(ok_all_subs), MADRS_all(ok_all_subs));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
scat_hdl.SizeData = mSize;
plot_hdl = plot(GSH_all_sorted, MADRS_f_GSH_fit.all);
plot_hdl.LineStyle = '--';
plot_hdl.LineWidth = lWidth;
plot_hdl.Color = grey;
xlabel('GSH');
ylabel('MADRS score');
legend_size(pSize);
