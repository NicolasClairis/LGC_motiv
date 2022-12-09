% check correlation between decrease of performance in max perf
% before/after each run and fatigue parameter

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
computerRoot = LGCM_root_paths;
dataRoot = [computerRoot, filesep, study_nm, filesep];
lab_or_home = 'home';
switch lab_or_home
    case 'lab'
        gitRoot = fullfile('C:','Users','clairis','Desktop');
    case 'home'
        gitRoot = fullfile('C:','Users','Loco','Documents');
end

%% load parameters
[mdlType, mdlN] = behavioral_model_selection;
[prm, mdlType, mdlN] = prm_extraction(study_nm, subject_id, mdlType, mdlN);
% extract parameter of interest
kFp = prm.kFp;

%% load fatigue rating before/after fMRI
[preIRMf_fatigue, postIRMf_fatigue,...
    deltaFatiguePrePostExp] = deal(NaN(1,NS));
% load fatigue
fatiguePath = fullfile(gitRoot,'GitHub',...
    'LGC_motiv','LGC_Motiv_analysis','behavior');
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    sub_idx = strcmp(excelReadGeneralFile.CID, sub_nm_bis);
    preIRMf_fatigue(iS)          = excelReadGeneralFile.FatiguePr__fMRI(sub_idx);
    postIRMf_fatigue(iS)          = excelReadGeneralFile.FatiguePost_fMRI(sub_idx);
    deltaFatiguePrePostExp(iS)  = postIRMf_fatigue(iS) - preIRMf_fatigue(iS);
end

%% correlate for each run and averaging across runs
goodSubs = (~isnan(kFp)).*(~isnan(deltaFatiguePrePostExp)) == 1;
[beta,~,stats] = glmfit(kFp(goodSubs), deltaFatiguePrePostExp(goodSubs),'normal');
kFp_sorted = sort(kFp(goodSubs));
deltaFrtg_fit = glmval(beta, kFp_sorted, 'identity');

%% display figure
% general figure infos
pSize = 30;
lSize = 2;
lWidth = 3;
grey = [143 143 143]./255;

fig;
hold on;
scat_hdl = scatter(kFp(goodSubs), deltaFatiguePrePostExp(goodSubs));
scat_hdl.LineWidth = lWidth;
scat_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kFp_sorted, deltaFrtg_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel('kFp');
ylabel('Fatigue rating post-IRMf - pre-IRMf');
legend_size(pSize);