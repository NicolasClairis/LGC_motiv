% script to prepare data in (NS*n) format in order to perform Structural Equation Modelling (SEM)
% NS = number of subjects and n = number of predictor variables (ie
% variables included in the SEM to predict the variable y).
%
% Please look at Usage notes 1.1.pdf from the toolbox:
% https://ch.mathworks.com/matlabcentral/fileexchange/60013-toolbox-for-structural-equation-modelling-sem
% to know how to implement the SEM properly

%% working dir
SEM_path = fullfile('P:','boulot','postdoc_CarmenSandi','results','SEM');

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% main variables
[aINS_Lac, aINS_fMRI,...
    THE,...
    PHE, kEp,...
    MHE, kEm] = deal(NaN(NS,1));

%% load brain metabolite
[metabolites] = metabolite_load(subject_id);
aINS_Lac(:) = metabolites.aINS.Lac;

%% load fMRI
% define fMRI GLM to work on
GLM_str = inputdlg('Which fMRI GLM?');
GLM = str2double(GLM_str);

% define fMRI ROI to use
[con_vec_all,...
    ~, ~, ~,...
    con_names,...
    ROI_coords, ttest_ROI] = ROI_extraction_group('study1', GLM,...
    subject_id, condition, 0);
n_cons = size(con_vec_all, 1);
n_ROIs = size(con_vec_all,3);
if n_ROIs > 1
    error(['more than 1 ROI selected, mediation script cannot work that way',...
        'please focus on one and do it again.']);
end
ROI_short_Nm = ROI_coords.ROI_nm.ROI_1_shortName;

% define regression estimate to look for in the fMRI GLM
con_idx = listdlg('promptstring','Which contrast?',...
    'ListString',con_names);
% con_nm = con_names{con_idx};
con_nm_str = inputdlg('Contrast short name?');
con_nm = con_nm_str{1};
aINS_fMRI(:) = con_vec_all(con_idx, :, 1);

%% load behavioral parameters (THE/PHE/MHE/kEp/kEm)
% THE/PHE/MHE
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, 0);
THE(:) = choice_hE.EpEm;
PHE(:) = choice_hE.Ep;
MHE(:) = choice_hE.Em;

% kEp/kEm
[prm] = prm_extraction(study_nm, subject_id);
for iS2 = 1:NS
    sub_nm2 = subject_id{iS2};
    bhvPrm_CID_idx_tmp = strcmp(prm.CID,sub_nm2);
    kEp(iS2) = prm.kEp(bhvPrm_CID_idx_tmp);
    kEm(iS2) = prm.kEm(bhvPrm_CID_idx_tmp);
end % subject loop

%% create matrix
mtrx_THE = [aINS_Lac, aINS_fMRI, THE, PHE, MHE, kEp, kEm];
% save data to load in R
save([SEM_path,filesep,'aINSLac_',ROI_short_Nm,...
    '_GLM',GLM_str{1},'_',con_nm,'_',num2str(NS),'subs.mat'],'mtrx_THE');

% %% load
% SEM_path = fullfile('P:','boulot','postdoc_CarmenSandi','results','SEM');
% ROI_short_Nm = 'MRS_aINS';
% GLM_str = inputdlg('Which fMRI GLM?');
% GLM = str2double(GLM_str);
% con_nm = 'EpEm_Ech';
% NS = 63;
% load([SEM_path,filesep,'aINSLac_',ROI_short_Nm,...
%     '_GLM',GLM_str{1},'_',con_nm,'_',num2str(NS),'subs.mat'],'mtrx_THE');
% aINS_Lac=mtrx_THE(:,1);
% aINS_fMRI=mtrx_THE(:,2);
% THE=mtrx_THE(:,3);
% PHE=mtrx_THE(:,4);
% MHE=mtrx_THE(:,5);
% kEp=mtrx_THE(:,6);
% kEm=mtrx_THE(:,7);
%% create matrix bis but after outlier removal
% remove "outliers" (><mean*3SD)
[~, ~, aINS_Lac_clean] = rmv_outliers_3sd(aINS_Lac);
[~, ~, aINS_fMRI_clean] = rmv_outliers_3sd(aINS_fMRI);
[~, ~, THE_clean] = rmv_outliers_3sd(THE);
[~, ~, MHE_clean] = rmv_outliers_3sd(MHE);
[~, ~, PHE_clean] = rmv_outliers_3sd(PHE);
[~, ~, kEp_clean] = rmv_outliers_3sd(kEp);

%% data optimized for THE
goodSubs_bis = ~isnan(aINS_Lac_clean).*~isnan(aINS_fMRI_clean).*~isnan(THE_clean) == 1;

mtrx_THE_no_outliers = [aINS_Lac(goodSubs_bis),...
    aINS_fMRI(goodSubs_bis),...
    THE(goodSubs_bis)];
NS_no_outliers = sum(goodSubs_bis);
% save data to load in R
save([SEM_path,filesep,'aINSLac_',ROI_short_Nm,...
    '_GLM',GLM_str{1},'_',con_nm,'_',num2str(NS_no_outliers),'subs_THE_no_outliers.mat'],'mtrx_THE_no_outliers');

%% data optimized for PHE
goodSubs_bis = ~isnan(aINS_Lac_clean).*~isnan(aINS_fMRI_clean).*~isnan(PHE_clean) == 1;

mtrx_PHE_no_outliers = [aINS_Lac(goodSubs_bis),...
    aINS_fMRI(goodSubs_bis),...
    PHE(goodSubs_bis)];
NS_no_outliers = sum(goodSubs_bis);
% save data to load in R
save([SEM_path,filesep,'aINSLac_',ROI_short_Nm,...
    '_GLM',GLM_str{1},'_',con_nm,'_',num2str(NS_no_outliers),'subs_PHE_no_outliers.mat'],'mtrx_PHE_no_outliers');

%% data optimized for MHE
goodSubs_bis = ~isnan(aINS_Lac_clean).*~isnan(aINS_fMRI_clean).*~isnan(MHE_clean) == 1;

mtrx_MHE_no_outliers = [aINS_Lac(goodSubs_bis),...
    aINS_fMRI(goodSubs_bis),...
    MHE(goodSubs_bis)];
NS_no_outliers = sum(goodSubs_bis);
% save data to load in R
save([SEM_path,filesep,'aINSLac_',ROI_short_Nm,...
    '_GLM',GLM_str{1},'_',con_nm,'_',num2str(NS_no_outliers),'subs_MHE_no_outliers.mat'],'mtrx_MHE_no_outliers');

%% data optimized for kEp
goodSubs_bis = ~isnan(aINS_Lac_clean).*~isnan(aINS_fMRI_clean).*~isnan(kEp_clean) == 1;
boxcox_kEp_clean = boxcox(kEp);

mtrx_kEp_no_outliers = [aINS_Lac(goodSubs_bis),...
    aINS_fMRI(goodSubs_bis),...
    kEp(goodSubs_bis),...
    boxcox_kEp_clean(goodSubs_bis)];
NS_no_outliers = sum(goodSubs_bis);
% save data to load in R
save([SEM_path,filesep,'aINSLac_',ROI_short_Nm,...
    '_GLM',GLM_str{1},'_',con_nm,'_',num2str(NS_no_outliers),'subs_kEp_no_outliers.mat'],'mtrx_kEp_no_outliers');