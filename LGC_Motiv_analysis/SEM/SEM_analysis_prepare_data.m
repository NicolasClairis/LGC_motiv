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
[blood_Lac, dmPFC_Lac, dmPFC_fMRI,...
    THE,...
    PHE, kEp,...
    MHE, kEm] = deal(NaN(NS,1));

%% load blood lactate
[Lac_tmp] = load_plasma_Lac;
for iS1 = 1:NS
    sub_nm1 = subject_id{iS1};
    Lac_CID_idx_tmp = strcmp(Lac_tmp.CID,['CID',sub_nm1]);
    blood_Lac(iS1) = Lac_tmp.Lac(Lac_CID_idx_tmp);
end % subject loop

%% load brain metabolite
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac(:) = metabolites.dmPFC.Lac;

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
dmPFC_fMRI(:) = con_vec_all(con_idx, :, 1);

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
mtrx_THE = [blood_Lac, dmPFC_Lac, dmPFC_fMRI, THE, PHE, MHE, kEp, kEm];
% save data to load in the SEM gui
save([SEM_path,filesep,'bloodLac_dmPFCLac_',ROI_short_Nm,...
    '_GLM',GLM_str{1},'_',con_nm,'_',num2str(NS),'subs.mat'],'mtrx_THE');