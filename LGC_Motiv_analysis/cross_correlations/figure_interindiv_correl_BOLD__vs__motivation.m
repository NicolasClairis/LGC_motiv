function[r_corr, pval, NS_goodS] = figure_interindiv_correl_BOLD__vs__motivation()
% [r_corr, pval] = figure_interindiv_correl_BOLD__vs__motivation
% figure_interindiv_correl_BOLD__vs__motivation will create a heatmap
% showing the correlation between high effort (HE), high physical effort
% (HPE), high mental effort (HME) choices and the sensitivity to physical
% effort (kEp) vs dmPFC/dACC and anterior insula (aIns) BOLD regression
% estimate for the effort chosen.
%
% OUTPUTS
% r_corr: structure with correlation coefficients for each correlation test
%
% pval: structure with p.values for each correlation test*
%
% NS_goodS: number of good subjects included in each analysis

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;

%% extract ROI
GLM = spm_input('GLM number',1,'e');

figure_folder = ['P:\boulot\postdoc_CarmenSandi\papers\Clairis_mediation_Lac\',...
    'figures\fig3_dmPFCdACC_aINS_fMRI_intra-individual\'];
dmPFC_filepath = [figure_folder,'GLM',num2str(GLM),'_prm_f_MRS_dmPFC_ROI_',num2str(NS),'subs.mat'];
aIns_filepath = [figure_folder,'GLM',num2str(GLM),'_prm_f_MRS_aINS_ROI_',num2str(NS),'subs.mat'];
if exist(dmPFC_filepath,'file') && exist(aIns_filepath,'file')
    dmPFC_BOLD_struct = load(dmPFC_filepath);
    dmPFC_BOLD_allCons = dmPFC_BOLD_struct.con_vec_all;
    aIns_BOLD_struct = load(aIns_filepath);
    aIns_BOLD_allCons = aIns_BOLD_struct.con_vec_all;
    con_names = dmPFC_BOLD_struct.con_names;
else
    % load both dmPFC and aIns here
    [con_vec_all,...
        ~, ~, ~,...
        con_names,...
        ROI_coords] = ROI_extraction_group(study_nm, GLM,...
        subject_id, condition, 0);
    dmPFC_BOLD_allCons = con_vec_all(:,:,1);
    aIns_BOLD_allCons = con_vec_all(:,:,2);
end
% extract data for the contrast of interest
con_idx = listdlg('promptstring','Which contrast?',...
    'ListString',con_names);
% con_nm = con_names{con_idx};
con_nm_str = inputdlg('Contrast short name?');
con_nm = con_nm_str{1};
[dmPFC_BOLD, aIns_BOLD] = deal(NaN(1,NS));
dmPFC_BOLD(:) = dmPFC_BOLD_allCons(con_idx, :, 1);
aIns_BOLD(:) = aIns_BOLD_allCons(con_idx, :, 1);

%% extract behavioral data (choices + parameters)
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, 0);
[prm] = prm_extraction(study_nm, subject_id);
% create variables of interest
HE_ch = choice_hE.EpEm;
HPE_ch = choice_hE.Ep;
HME_ch = choice_hE.Em;
kEp = prm.kEp;

%% perform the correlations
for iRawCorr = 1:2
    switch iRawCorr
        case 1
            raw_or_corr_nm = 'raw';
        case 2
            raw_or_corr_nm = 'filtered';
    end

    % correlations with dmPFC/dACC BOLD
    goodS.(raw_or_corr_nm).HE_f_dmPFC = ~isnan(dmPFC_BOLD.*HE_ch);
    [r_corr.(raw_or_corr_nm).HE_f_dmPFC, pval.(raw_or_corr_nm).HE_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC)', HE_ch(goodS.(raw_or_corr_nm).HE_f_dmPFC)');
    goodS.(raw_or_corr_nm).HPE_f_dmPFC = ~isnan(dmPFC_BOLD.*HPE_ch);
    [r_corr.(raw_or_corr_nm).HPE_f_dmPFC, pval.(raw_or_corr_nm).HPE_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_dmPFC)');
    goodS.(raw_or_corr_nm).HME_f_dmPFC = ~isnan(dmPFC_BOLD.*HME_ch);
    [r_corr.(raw_or_corr_nm).HME_f_dmPFC, pval.(raw_or_corr_nm).HME_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).HME_f_dmPFC)', HME_ch(goodS.(raw_or_corr_nm).HME_f_dmPFC)');
    goodS.(raw_or_corr_nm).kEp_f_dmPFC = ~isnan(dmPFC_BOLD.*kEp);
    [r_corr.(raw_or_corr_nm).kEp_f_dmPFC, pval.(raw_or_corr_nm).kEp_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).kEp_f_dmPFC)', kEp(goodS.(raw_or_corr_nm).kEp_f_dmPFC)');
    % correlations with aIns BOLD
    goodS.(raw_or_corr_nm).HE_f_aIns = ~isnan(aIns_BOLD.*HE_ch);
    [r_corr.(raw_or_corr_nm).HE_f_aIns, pval.(raw_or_corr_nm).HE_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns)', HE_ch(goodS.(raw_or_corr_nm).HE_f_aIns)');
    goodS.(raw_or_corr_nm).HPE_f_aIns = ~isnan(aIns_BOLD.*HPE_ch);
    [r_corr.(raw_or_corr_nm).HPE_f_aIns, pval.(raw_or_corr_nm).HPE_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_aIns)');
    goodS.(raw_or_corr_nm).HME_f_aIns = ~isnan(aIns_BOLD.*HME_ch);
    [r_corr.(raw_or_corr_nm).HME_f_aIns, pval.(raw_or_corr_nm).HME_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).HME_f_aIns)', HME_ch(goodS.(raw_or_corr_nm).HME_f_aIns)');
    goodS.(raw_or_corr_nm).kEp_f_aIns = ~isnan(aIns_BOLD.*kEp);
    [r_corr.(raw_or_corr_nm).kEp_f_aIns, pval.(raw_or_corr_nm).kEp_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).kEp_f_aIns)', kEp(goodS.(raw_or_corr_nm).kEp_f_aIns)');

    % extract number of good subject for each
    fields_to_check = fieldnames(goodS.(raw_or_corr_nm));
    for iF = 1:length(fields_to_check)
        field_nm = fields_to_check{iF};
        NS_goodS.(raw_or_corr_nm).(field_nm) = sum(goodS.(raw_or_corr_nm).(field_nm));
    end % loop over fields to check
end % loop over raw vs 3*SD corrected data

%% display correlation in a nice correlation matrix

end % function