% function aimed at extracting ROI in a given GLM based on the level of a
% metabolite of interest (to select) in one of the two brain areas of study
% 1 that can be checked (dmPFC/aIns)

%% define metabolite and ROI you want to focus on
% ROI
ROIs = {'dmPFC','aIns'};
nROIs = length(ROIs);
ROI_idx = spm_input('Metabolites in which brain area?',1,'m',...
    ROIs,1:nROIs,0);
ROI_nm = ROIs{ROI_idx};
% select metabolite of interest
metabolites = {'Mac','Ala','Asp','PCho','Cr','PCr','GABA',...
    'Gln','Glu','GSH','Gly','Ins','Lac','NAA','Scyllo','Tau',...
    'Asc','Glc','NAAG','GPC','PE','Ser',...
    'NAA_NAAG','Glu_Gln','GPC_PCho','Cr_PCr','Gly_Ins','Gln_div_Glu'};
n_met = length(metabolites);
metabolite_idx = spm_input('Which metabolite to focus on?',1,'m',...
    metabolites,1:n_met,0);
metabolite_nm = metabolites{metabolite_idx};
% define GLM number
GLM = spm_input('GLM number',1,'e');

%% define all subjects
[subject_id, NS] = LGCM_subject_selection('study1');

%% extract all metabolites
[metabolites] = metabolite_load(subject_id);
% focus on metabolite and brain area selected
metabolite_allSubs = metabolites.(ROI_nm).(metabolite_nm);
% perform a median split based on the metabolite selected in the ROI
% selected
med_metabolite_allSubs = median(metabolite_allSubs,'omitnan');
% extract index of participants with low or high level of metabolites
low_met_subs = metabolite_allSubs <= med_metabolite_allSubs;
high_met_subs = metabolite_allSubs > med_metabolite_allSubs;

%% extract ROI
[con_vec_all,...
    ~, ~, ~,...
    con_names,...
    ROI_coords, ttest_ROI] = ROI_extraction_group('study1', GLM, subject_id, 0);
n_cons = size(con_vec_all, 1);

%% extract ROI data split according to the metabolite levels
[con_avg_lowMet, con_avg_highMet,...
    con_sem_lowMet, con_sem_highMet,...
    con_std_lowMet, con_std_highMet,...
    ttest_pval_lowMet_vs_highMet, ttest_tval_lowMet_vs_highMet] = deal(NaN(1,n_cons));
for iCon = 1:n_cons
    [con_avg_lowMet(iCon),...
        con_sem_lowMet(iCon),...
        con_std_lowMet(iCon)] = mean_sem_sd(con_vec_all(iCon, low_met_subs),2);
    [con_avg_highMet(iCon),...
        con_sem_lowMet(iCon),...
        con_std_lowMet(iCon)] = mean_sem_sd(con_vec_all(iCon, high_met_subs),2);
    error('careful current test is a paired t.test, you need an unpaired t.test here');
    [~,ttest_lowMet_vs_highMet,~,stats_tmp] = ttest(con_vec_all(iCon, low_met_subs), con_vec_all(iCon, high_met_subs))
    ttest_tval_lowMet_vs_highMet = stats_tmp.tstat;
end % contrast loop

%% figure

