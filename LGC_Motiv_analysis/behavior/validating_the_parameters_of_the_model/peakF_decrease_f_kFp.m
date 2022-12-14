% peakF_decrease_f_kFp will test whether the decrease in the peak force for
% each effort chosen over time correlates, across individuals, with the
% modeled sensitivity to physical fatigue (kFp).

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
computerRoot = LGCM_root_paths;
dataRoot = [computerRoot, filesep, study_nm, filesep];

%% load parameters
[mdlType, mdlN] = behavioral_model_selection;
prm = prm_extraction(study_nm, subject_id, mdlType, mdlN);
% extract parameter of interest
kFp = prm.kFp;

%% load slope for peak force decrease with time
figDisp = 0;
[b_peakF_perEch_f_time,~,b_peakF_allE_f_time] = peakF_f_time(study_nm, subject_id, condition, figDisp);

%% correlate kFp to peak force decrease with time slope
Ech_levels = 0:3;
n_Ech_lvl = length(Ech_levels);
for iEch = Ech_levels
    Ech_nm = ['Ech',num2str(iEch)];
    jEch = iEch + 1;
    goodSubs.(Ech_nm) = (~isnan(kFp)).*(~isnan(b_peakF_perEch_f_time(jEch,:))) == 1;
    [beta.(Ech_nm),~,stats.(Ech_nm)] = glmfit(kFp(goodSubs.(Ech_nm)), b_peakF_perEch_f_time(jEch,goodSubs.(Ech_nm)), 'normal');
    kFp_sorted.(Ech_nm) = sort(kFp(goodSubs.(Ech_nm)));
    b_peakF_fit.(Ech_nm) = glmval(beta.(Ech_nm), kFp_sorted.(Ech_nm), 'identity');
    pval.(Ech_nm) = stats.(Ech_nm).p;
end % effort chosen loop

% what about across all efforts?
goodSubs.allE = (~isnan(kFp)).*(~isnan(b_peakF_allE_f_time)) == 1;
[beta.allE,~,stats.allE] = glmfit(kFp(goodSubs.allE), b_peakF_allE_f_time(goodSubs.allE), 'normal');
kFp_sorted.allE = sort(kFp(goodSubs.allE));
b_peakF_fit.allE = glmval(beta.allE, kFp_sorted.allE, 'identity');
pval.allE = stats.allE.p;

%% figure
pSize = 30;
lWidth = 3;
grey = [143 143 143]./255;

% split per effort level
fig;
for iEch = Ech_levels
    Ech_nm = ['Ech',num2str(iEch)];
    jEch = iEch + 1;
    subplot(2,2,jEch);
    hold on;
    scat_hdl = scatter(kFp(goodSubs.(Ech_nm)), b_peakF_perEch_f_time(jEch,goodSubs.(Ech_nm)));
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = 'k';
    fit_hdl = plot(kFp_sorted.(Ech_nm), b_peakF_fit.(Ech_nm));
    fit_hdl.LineWidth = lWidth;
    fit_hdl.LineStyle = '--';
    fit_hdl.Color = grey;
    xlabel('kFp');
    ylabel({'peakF decrease with time';...
        Ech_nm});
    legend_size(pSize);
end

% average across effort levels
fig;
hold on;
scat_hdl = scatter(kFp(goodSubs.allE), b_peakF_allE_f_time(jEch,goodSubs.allE));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kFp_sorted.allE, b_peakF_fit.allE);
fit_hdl.LineWidth = lWidth;
fit_hdl.LineStyle = '--';
fit_hdl.Color = grey;
xlabel('kFp');
ylabel('peakF decrease with time');
legend_size(pSize);