%% script to test correlation between MADRS-S scores and metabolite levels
%
% designed by N.Clairis - 2022

%% working directories
computerRoot = LGCM_root_paths;
study_nm = 'study1';

%% load MADRS-S scores
[excelReadQuestionnairesFile, MADRS_S_sub_CID_list] = load_questionnaires_data();
MADRS_S_fullList = excelReadQuestionnairesFile.MADRS_SCorrected;

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

%% extract relevant subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% extract all metabolites
[metabolites] = metabolite_load(subject_id);
% focus on metabolite and brain area selected
metabolite_allSubs = metabolites.(ROI_nm).(metabolite_nm);

[MADRS_S_score] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    MADRS_S_score(iS) = MADRS_S_fullList(strcmp(MADRS_S_sub_CID_list, sub_nm));
end % subject list
%% perform GLM and correlation between metabolites and MADRS-S
var_nm = [metabolite_nm,'_f_MADRS'];
% remove subjects with NaN values
goodSubs = ~isnan(metabolite_allSubs);
[betas.(var_nm), ~,stats_tmp] = glmfit(MADRS_S_score(goodSubs), metabolite_allSubs(goodSubs),'normal');
pval.(var_nm) = stats_tmp.p;
metabolite_f_MADRS_S_fit.(var_nm) = NaN(1,NS);
metabolite_f_MADRS_S_fit.(var_nm)(goodSubs) = glmval(betas.(var_nm),metabolite_allSubs(goodSubs),'identity');
 
%% display results
pSize = 50;
fig;
hdl=scatter(MADRS_S_score, metabolite_allSubs,...
    'MarkerEdgeColor','k','SizeData',100,'LineWidth',3);
hold on;
[MADRS_S_score_sorted, idx_sorted] = sort(MADRS_S_score);
plot(MADRS_S_score_sorted, metabolite_f_MADRS_S_fit.(var_nm)(idx_sorted),...
    'LineStyle','--','LineWidth',3);
ylim_vals = ylim();
line([4 4],ylim_vals,...
    'LineStyle','-','LineWidth',3,'Color','k');
line([6 6],ylim_vals,...
    'LineStyle','-','LineWidth',3,'Color','k');
xlabel('MADRS-S score');
ylabel(metabolite_nm);
legend_size(pSize);