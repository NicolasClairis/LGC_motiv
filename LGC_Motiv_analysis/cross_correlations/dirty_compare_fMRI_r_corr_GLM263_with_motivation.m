%% load data to compare correlation between (dmPFC/dACC Ech fMRI regression 
% estimate) and HE choices to (dmPFC/dACC SVch fMRI regression estimate) 
% and (dmPFC/dACC DT fMRI regression estimate) correlations with HE choices
% => Need to first load the data carefully to ensure that the same number
% of subjects is included in all analysis before extracting the correlation
% coefficients.

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;

%% extract ROI
GLM = 263;

figure_folder = ['P:\boulot\postdoc_CarmenSandi\papers\Clairis_mediation_Lac\',...
    'figures\fig3_dmPFCdACC_aINS_fMRI\'];
dmPFC_filepath = [figure_folder,'GLM',num2str(GLM),'_',num2str(NS),'subs_dmPFC_ROI.mat'];
aIns_filepath = [figure_folder,'GLM',num2str(GLM),'_',num2str(NS),'subs_aIns_ROI.mat'];
if exist(dmPFC_filepath,'file') && exist(aIns_filepath,'file')
    dmPFC_BOLD_struct = load(dmPFC_filepath);
    dmPFC_BOLD_allCons = dmPFC_BOLD_struct.con_vec_all;
    aIns_BOLD_struct = load(aIns_filepath);
    aIns_BOLD_allCons = aIns_BOLD_struct.con_vec_all;
    con_names = dmPFC_BOLD_struct.con_names;
    close all; % in case some figures were saved as well
else
    % extract dmPFC and aIns here
    % first dmPFC
    [dmPFC_BOLD_allCons,...
        ~, ~, ~,...
        con_names,...
        ROI_coords] = ROI_extraction_group(study_nm, GLM,...
        subject_id, condition, 0);
    % then aIns
    [aIns_BOLD_allCons,...
        ~, ~, ~,...
        con_names,...
        ROI_coords] = ROI_extraction_group(study_nm, GLM,...
        subject_id, condition, 0);
end
% extract data for the contrasts of interest
Ech_con_idx = strcmp(con_names,'Ep+Em REG choice RP E: effort chosen');
SVch_con_idx = strcmp(con_names,'Ep+Em REG choice RP E: NVch-NVunch');
DT_con_idx = strcmp(con_names,'Ep+Em REG choice RP E: RT');
[dmPFC_Ech_BOLD, aIns_Ech_BOLD,...
    dmPFC_SVch_BOLD, aIns_SVch_BOLD,...
    dmPFC_DT_BOLD, aIns_DT_BOLD] = deal(NaN(1,NS));
dmPFC_Ech_BOLD(:) = dmPFC_BOLD_allCons(Ech_con_idx, :, 1);
dmPFC_SVch_BOLD(:) = dmPFC_BOLD_allCons(SVch_con_idx, :, 1);
dmPFC_DT_BOLD(:) = dmPFC_BOLD_allCons(DT_con_idx, :, 1);
aIns_Ech_BOLD(:) = aIns_BOLD_allCons(Ech_con_idx, :, 1);
aIns_SVch_BOLD(:) = aIns_BOLD_allCons(SVch_con_idx, :, 1);
aIns_DT_BOLD(:) = aIns_BOLD_allCons(DT_con_idx, :, 1);

%% extract behavioral data (choices + parameters)
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, 0);
[prm] = prm_extraction(study_nm, subject_id);
prm_names = fieldnames(prm);
prm_names(strcmp(prm_names,'CID')) = [];
nPrm = length(prm_names);
% create variables of interest
HE_ch = choice_hE.EpEm;
HPE_ch = choice_hE.Ep;
HME_ch = choice_hE.Em;
kEp = prm.kEp;
kEm = prm.kEm;
kR = prm.kR;
kP = prm.kP;
kLm = prm.kLm;
kFp = prm.kFp;
kBias = prm.kBias;

%% use 'raw' or 'filtered' data?
raw_or_corr_nm = 'filtered';

%% filter subjects on each variable
% HE <=> dmPFC
goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech = filter_fn(raw_or_corr_nm, dmPFC_Ech_BOLD, HE_ch);
goodS.(raw_or_corr_nm).HE_f_dmPFC_SVch = filter_fn(raw_or_corr_nm, dmPFC_SVch_BOLD, HE_ch);
goodS.(raw_or_corr_nm).HE_f_dmPFC_DT = filter_fn(raw_or_corr_nm, dmPFC_DT_BOLD, HE_ch);
% extract common subjects across all measures
goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch = (goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).HE_f_dmPFC_SVch) == 1;
goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT = (goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).HE_f_dmPFC_DT) == 1;

% HE <=> aIns
goodS.(raw_or_corr_nm).HE_f_aIns_Ech = filter_fn(raw_or_corr_nm, aIns_Ech_BOLD, HE_ch);
goodS.(raw_or_corr_nm).HE_f_aIns_SVch = filter_fn(raw_or_corr_nm, aIns_SVch_BOLD, HE_ch);
goodS.(raw_or_corr_nm).HE_f_aIns_DT = filter_fn(raw_or_corr_nm, aIns_DT_BOLD, HE_ch);
% extract common subjects across all measures
goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch = (goodS.(raw_or_corr_nm).HE_f_aIns_Ech.*goodS.(raw_or_corr_nm).HE_f_aIns_SVch) == 1;
goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT = (goodS.(raw_or_corr_nm).HE_f_aIns_Ech.*goodS.(raw_or_corr_nm).HE_f_aIns_DT) == 1;

% HPE <=> dmPFC
goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech = filter_fn(raw_or_corr_nm, dmPFC_Ech_BOLD, HPE_ch);
goodS.(raw_or_corr_nm).HPE_f_dmPFC_SVch = filter_fn(raw_or_corr_nm, dmPFC_SVch_BOLD, HPE_ch);
goodS.(raw_or_corr_nm).HPE_f_dmPFC_DT = filter_fn(raw_or_corr_nm, dmPFC_DT_BOLD, HPE_ch);
% extract common subjects across all measures
goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch = (goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).HPE_f_dmPFC_SVch) == 1;
goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT = (goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).HPE_f_dmPFC_DT) == 1;

% HPE <=> aIns
goodS.(raw_or_corr_nm).HPE_f_aIns_Ech = filter_fn(raw_or_corr_nm, aIns_Ech_BOLD, HPE_ch);
goodS.(raw_or_corr_nm).HPE_f_aIns_SVch = filter_fn(raw_or_corr_nm, aIns_SVch_BOLD, HPE_ch);
goodS.(raw_or_corr_nm).HPE_f_aIns_DT = filter_fn(raw_or_corr_nm, aIns_DT_BOLD, HPE_ch);
% extract common subjects across all measures
goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch = (goodS.(raw_or_corr_nm).HPE_f_aIns_Ech.*goodS.(raw_or_corr_nm).HPE_f_aIns_SVch) == 1;
goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT = (goodS.(raw_or_corr_nm).HPE_f_aIns_Ech.*goodS.(raw_or_corr_nm).HPE_f_aIns_DT) == 1;


% behavioral parameters <=> dmPFC
goodS.(raw_or_corr_nm).kEp_f_dmPFC_Ech = filter_fn(raw_or_corr_nm, dmPFC_Ech_BOLD, kEp);
goodS.(raw_or_corr_nm).kEm_f_dmPFC_Ech = filter_fn(raw_or_corr_nm, dmPFC_Ech_BOLD, kEm);
goodS.(raw_or_corr_nm).kR_f_dmPFC_Ech = filter_fn(raw_or_corr_nm, dmPFC_Ech_BOLD, kR);
goodS.(raw_or_corr_nm).kP_f_dmPFC_Ech = filter_fn(raw_or_corr_nm, dmPFC_Ech_BOLD, kP);
goodS.(raw_or_corr_nm).kLm_f_dmPFC_Ech = filter_fn(raw_or_corr_nm, dmPFC_Ech_BOLD, kLm);
goodS.(raw_or_corr_nm).kFp_f_dmPFC_Ech = filter_fn(raw_or_corr_nm, dmPFC_Ech_BOLD, kFp);
goodS.(raw_or_corr_nm).kBias_f_dmPFC_Ech = filter_fn(raw_or_corr_nm, dmPFC_Ech_BOLD, kBias);
% extract common subjects across all measures
goodS.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech = (goodS.(raw_or_corr_nm).kEp_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).kEm_f_dmPFC_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech = (goodS.(raw_or_corr_nm).kEp_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).kR_f_dmPFC_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech = (goodS.(raw_or_corr_nm).kEp_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).kP_f_dmPFC_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech = (goodS.(raw_or_corr_nm).kEp_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).kLm_f_dmPFC_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech = (goodS.(raw_or_corr_nm).kEp_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).kFp_f_dmPFC_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech = (goodS.(raw_or_corr_nm).kEp_f_dmPFC_Ech.*goodS.(raw_or_corr_nm).kBias_f_dmPFC_Ech) == 1;


% behavioral parameters <=> aIns
goodS.(raw_or_corr_nm).kEp_f_aIns_Ech = filter_fn(raw_or_corr_nm, aIns_Ech_BOLD, kEp);
goodS.(raw_or_corr_nm).kEm_f_aIns_Ech = filter_fn(raw_or_corr_nm, aIns_Ech_BOLD, kEm);
goodS.(raw_or_corr_nm).kR_f_aIns_Ech = filter_fn(raw_or_corr_nm, aIns_Ech_BOLD, kR);
goodS.(raw_or_corr_nm).kP_f_aIns_Ech = filter_fn(raw_or_corr_nm, aIns_Ech_BOLD, kP);
goodS.(raw_or_corr_nm).kLm_f_aIns_Ech = filter_fn(raw_or_corr_nm, aIns_Ech_BOLD, kLm);
goodS.(raw_or_corr_nm).kFp_f_aIns_Ech = filter_fn(raw_or_corr_nm, aIns_Ech_BOLD, kFp);
goodS.(raw_or_corr_nm).kBias_f_aIns_Ech = filter_fn(raw_or_corr_nm, aIns_Ech_BOLD, kBias);
% extract common subjects across all measures
goodS.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech = (goodS.(raw_or_corr_nm).kEp_f_aIns_Ech.*goodS.(raw_or_corr_nm).kEm_f_aIns_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech = (goodS.(raw_or_corr_nm).kEp_f_aIns_Ech.*goodS.(raw_or_corr_nm).kR_f_aIns_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech = (goodS.(raw_or_corr_nm).kEp_f_aIns_Ech.*goodS.(raw_or_corr_nm).kP_f_aIns_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech = (goodS.(raw_or_corr_nm).kEp_f_aIns_Ech.*goodS.(raw_or_corr_nm).kLm_f_aIns_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech = (goodS.(raw_or_corr_nm).kEp_f_aIns_Ech.*goodS.(raw_or_corr_nm).kFp_f_aIns_Ech) == 1;
goodS.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech = (goodS.(raw_or_corr_nm).kEp_f_aIns_Ech.*goodS.(raw_or_corr_nm).kBias_f_aIns_Ech) == 1;

%% perform correlations of interest
%% [HE <=> dmPFC/dACC Ech] vs [HE <=> dmPFC/dACC SVch]
NS_goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch = sum(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch);
% Ech fMRI
[r_corr.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch.Ech] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch)', HE_ch(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch)');
% SVch fMRI
[r_corr.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch.SVch] = corr(dmPFC_SVch_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch)', HE_ch(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch)');
% Ech vs SVch
[r_corr.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch.Ech_vs_SVch] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch)', dmPFC_SVch_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_SVch)');

%% [HE <=> dmPFC/dACC Ech] vs [HE <=> dmPFC/dACC DT]
NS_goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT = sum(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT);
% Ech fMRI
[r_corr.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT.Ech] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT)', HE_ch(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT)');
% DT fMRI
[r_corr.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT.DT] = corr(dmPFC_DT_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT)', HE_ch(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT)');
% Ech vs DT
[r_corr.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT.Ech_vs_DT] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT)', dmPFC_DT_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC_Ech_vs_DT)');


%% [HE <=> aIns Ech] vs [HE <=> aIns SVch]
NS_goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch = sum(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch);
% Ech fMRI
[r_corr.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch.Ech] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch)', HE_ch(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch)');
% SVch fMRI
[r_corr.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch.SVch] = corr(aIns_SVch_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch)', HE_ch(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch)');
% Ech vs SVch
[r_corr.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch.Ech_vs_SVch] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch)', aIns_SVch_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_SVch)');

%% [HE <=> aIns Ech] vs [HE <=> aIns DT]
NS_goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT = sum(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT);
% Ech fMRI
[r_corr.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT.Ech] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT)', HE_ch(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT)');
% DT fMRI
[r_corr.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT.DT] = corr(aIns_DT_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT)', HE_ch(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT)');
% Ech vs DT
[r_corr.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT.Ech_vs_DT] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT)', aIns_DT_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns_Ech_vs_DT)');


%% [HPE <=> dmPFC/dACC Ech] vs [HPE <=> dmPFC/dACC SVch]
NS_goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch = sum(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch);
% Ech fMRI
[r_corr.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch.Ech] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch)');
% SVch fMRI
[r_corr.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch.SVch] = corr(dmPFC_SVch_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch)');
% Ech vs SVch
[r_corr.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch.Ech_vs_SVch] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch)', dmPFC_SVch_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_SVch)');

%% [HPE <=> dmPFC/dACC Ech] vs [HPE <=> dmPFC/dACC DT]
NS_goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT = sum(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT);
% Ech fMRI
[r_corr.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT.Ech] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT)');
% DT fMRI
[r_corr.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT.DT] = corr(dmPFC_DT_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT)');
% Ech vs DT
[r_corr.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT.Ech_vs_DT] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT)', dmPFC_DT_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Ech_vs_DT)');

%% [HPE <=> aIns Ech] vs [HPE <=> aIns SVch]
NS_goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch = sum(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch);
% Ech fMRI
[r_corr.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch.Ech] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch)');
% SVch fMRI
[r_corr.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch.SVch] = corr(aIns_SVch_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch)');
% Ech vs SVch
[r_corr.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch.Ech_vs_SVch] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch)', aIns_SVch_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_SVch)');

%% [HPE <=> aIns Ech] vs [HPE <=> aIns DT]
NS_goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT = sum(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT);
% Ech fMRI
[r_corr.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT.Ech] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT)');
% DT fMRI
[r_corr.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT.DT] = corr(aIns_DT_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT)');
% Ech vs DT
[r_corr.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT.Ech_vs_DT] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT)', aIns_DT_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns_Ech_vs_DT)');


%% [kEp <=> dmPFC/dACC Ech] vs [kEm <=> dmPFC/dACC Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech.kEp] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech)');
% kEm
[r_corr.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech.kEm] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech)', kEm(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech)');
% kEp vs kEm
[r_corr.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech.kEp_vs_kEm] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech)', kEm(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_dmPFC_Ech)');

%% [kEp <=> dmPFC/dACC Ech] vs [kR <=> dmPFC/dACC Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech.kEp] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech)');
% kR
[r_corr.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech.kR] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech)', kR(goodS.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech)');
% kEp vs kR
[r_corr.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech.kEp_vs_kR] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech)', kR(goodS.(raw_or_corr_nm).kEp_vs_kR_f_dmPFC_Ech)');

%% [kEp <=> dmPFC/dACC Ech] vs [kP <=> dmPFC/dACC Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech.kEp] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech)');
% kP
[r_corr.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech.kP] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech)', kP(goodS.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech)');
% kEp vs kP
[r_corr.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech.kEp_vs_kP] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech)', kP(goodS.(raw_or_corr_nm).kEp_vs_kP_f_dmPFC_Ech)');

%% [kEp <=> dmPFC/dACC Ech] vs [kLm <=> dmPFC/dACC Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech.kEp] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech)');
% kLm
[r_corr.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech.kLm] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech)', kLm(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech)');
% kEp vs kLm
[r_corr.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech.kEp_vs_kLm] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech)', kLm(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_dmPFC_Ech)');

%% [kEp <=> dmPFC/dACC Ech] vs [kFp <=> dmPFC/dACC Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech.kEp] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech)');
% kFp
[r_corr.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech.kFp] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech)', kFp(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech)');
% kEp vs kFp
[r_corr.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech.kEp_vs_kFp] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech)', kFp(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_dmPFC_Ech)');

%% [kEp <=> dmPFC/dACC Ech] vs [kBias <=> dmPFC/dACC Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech.kEp] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech)');
% kBias
[r_corr.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech.kBias] = corr(dmPFC_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech)', kBias(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech)');
% kEp vs kBias
[r_corr.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech.kEp_vs_kBias] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech)', kBias(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_dmPFC_Ech)');



%% [kEp <=> aIns Ech] vs [kEm <=> aIns Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech.kEp] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech)');
% kEm
[r_corr.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech.kEm] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech)', kEm(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech)');
% kEp vs kEm
[r_corr.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech.kEp_vs_kEm] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech)', kEm(goodS.(raw_or_corr_nm).kEp_vs_kEm_f_aIns_Ech)');

%% [kEp <=> aIns Ech] vs [kR <=> aIns Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech.kEp] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech)');
% kR
[r_corr.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech.kR] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech)', kR(goodS.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech)');
% kEp vs kR
[r_corr.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech.kEp_vs_kR] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech)', kR(goodS.(raw_or_corr_nm).kEp_vs_kR_f_aIns_Ech)');

%% [kEp <=> aIns Ech] vs [kP <=> aIns Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech.kEp] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech)');
% kP
[r_corr.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech.kP] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech)', kP(goodS.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech)');
% kEp vs kP
[r_corr.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech.kEp_vs_kP] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech)', kP(goodS.(raw_or_corr_nm).kEp_vs_kP_f_aIns_Ech)');

%% [kEp <=> aIns Ech] vs [kLm <=> aIns Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech.kEp] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech)');
% kLm
[r_corr.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech.kLm] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech)', kLm(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech)');
% kEp vs kLm
[r_corr.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech.kEp_vs_kLm] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech)', kLm(goodS.(raw_or_corr_nm).kEp_vs_kLm_f_aIns_Ech)');

%% [kEp <=> aIns Ech] vs [kFp <=> aIns Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech.kEp] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech)');
% kFp
[r_corr.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech.kFp] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech)', kFp(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech)');
% kEp vs kFp
[r_corr.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech.kEp_vs_kFp] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech)', kFp(goodS.(raw_or_corr_nm).kEp_vs_kFp_f_aIns_Ech)');

%% [kEp <=> aIns Ech] vs [kBias <=> aIns Ech]
NS_goodS.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech = sum(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech);
% kEp
[r_corr.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech.kEp] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech)', kEp(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech)');
% kBias
[r_corr.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech.kBias] = corr(aIns_Ech_BOLD(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech)', kBias(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech)');
% kEp vs kBias
[r_corr.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech.kEp_vs_kBias] = corr(kEp(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech)', kBias(goodS.(raw_or_corr_nm).kEp_vs_kBias_f_aIns_Ech)');
