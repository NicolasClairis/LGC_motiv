function[r_corr, pval, NS_goodS] = figure_interindiv_correl_Lac__vs__motivation()
% [r_corr, pval, NS_goodS] = figure_interindiv_correl_Lac__vs__motivation
% figure_interindiv_correl_Lac__vs__motivation will create a heatmap
% showing the correlation between high effort (HE), high physical effort
% (HPE), high mental effort (HME) choices and the sensitivity to physical
% effort (kEp) vs plasma, dmPFC/dACC and anterior insula (aIns) lactate.
%
% OUTPUTS
% r_corr: structure with correlation coefficients for each correlation test
%
% pval: structure with p.values for each correlation test
%
% NS_goodS: number of good subjects included in each analysis

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;

%% extract lactate levels
% load plasma lactate
[plasmaLac_struct] = load_plasma_Lac(subject_id);
plasma_Lac = plasmaLac_struct.Lac;
% load brain lactate
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac = metabolites.dmPFC.Lac;
aIns_Lac = metabolites.aIns.Lac;

%% extract behavioral data (choices + parameters)
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, 0);
[prm] = prm_extraction(study_nm, subject_id);
% create variables of interest
HE_ch = choice_hE.EpEm;
HPE_ch = choice_hE.Ep;
HME_ch = choice_hE.Em;
kEp = prm.kEp;

%% figure parameters
[pSize, ~, col] = general_fig_prm;
% define which colormap you want to use (see full list here if you are not
% happy with the selection:
% https://ch.mathworks.com/help/matlab/ref/colormap.html)
% color_range_choices = 'hot';
% color_range_choices = 'turbo';
% color_range_choices = 'jet';
color_range_choices = redblue(45);

% correlation range
corr_range = [-0.65 0.65];

%% perform the correlations
for iRawCorr = 1:2
    switch iRawCorr
        case 1
            raw_or_corr_nm = 'raw';
        case 2
            raw_or_corr_nm = 'filtered';
    end
    
    % correlations with lactate
    % plasma lactate
    goodS.(raw_or_corr_nm).HE_f_plasma_Lac = filter_fn(raw_or_corr_nm, plasma_Lac, HE_ch);
    [r_corr.(raw_or_corr_nm).HE_f_plasma_Lac,...
        pval.(raw_or_corr_nm).HE_f_plasma_Lac] = corr(plasma_Lac(goodS.(raw_or_corr_nm).HE_f_plasma_Lac)', HE_ch(goodS.(raw_or_corr_nm).HE_f_plasma_Lac)');
    goodS.(raw_or_corr_nm).HPE_f_plasma_Lac = filter_fn(raw_or_corr_nm, plasma_Lac, HPE_ch);
    [r_corr.(raw_or_corr_nm).HPE_f_plasma_Lac,...
        pval.(raw_or_corr_nm).HPE_f_plasma_Lac] = corr(plasma_Lac(goodS.(raw_or_corr_nm).HPE_f_plasma_Lac)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_plasma_Lac)');
    goodS.(raw_or_corr_nm).HME_f_plasma_Lac = filter_fn(raw_or_corr_nm, plasma_Lac, HME_ch);
    [r_corr.(raw_or_corr_nm).HME_f_plasma_Lac,...
        pval.(raw_or_corr_nm).HME_f_plasma_Lac] = corr(plasma_Lac(goodS.(raw_or_corr_nm).HME_f_plasma_Lac)', HME_ch(goodS.(raw_or_corr_nm).HME_f_plasma_Lac)');
    goodS.(raw_or_corr_nm).kEp_f_plasma_Lac = filter_fn(raw_or_corr_nm, plasma_Lac, kEp);
    [r_corr.(raw_or_corr_nm).kEp_f_plasma_Lac,...
        pval.(raw_or_corr_nm).kEp_f_plasma_Lac] = corr(plasma_Lac(goodS.(raw_or_corr_nm).kEp_f_plasma_Lac)', kEp(goodS.(raw_or_corr_nm).kEp_f_plasma_Lac)');
    
    % dmPFC/dACC lactate
    goodS.(raw_or_corr_nm).HE_f_dmPFC_Lac = filter_fn(raw_or_corr_nm, dmPFC_Lac, HE_ch);
    [r_corr.(raw_or_corr_nm).HE_f_dmPFC_Lac,...
        pval.(raw_or_corr_nm).HE_f_dmPFC_Lac] = corr(dmPFC_Lac(goodS.(raw_or_corr_nm).HE_f_dmPFC_Lac)', HE_ch(goodS.(raw_or_corr_nm).HE_f_dmPFC_Lac)');
    goodS.(raw_or_corr_nm).HPE_f_dmPFC_Lac = filter_fn(raw_or_corr_nm, dmPFC_Lac, HPE_ch);
    [r_corr.(raw_or_corr_nm).HPE_f_dmPFC_Lac,...
        pval.(raw_or_corr_nm).HPE_f_dmPFC_Lac] = corr(dmPFC_Lac(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Lac)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_dmPFC_Lac)');
    goodS.(raw_or_corr_nm).HME_f_dmPFC_Lac = filter_fn(raw_or_corr_nm, dmPFC_Lac, HME_ch);
    [r_corr.(raw_or_corr_nm).HME_f_dmPFC_Lac,...
        pval.(raw_or_corr_nm).HME_f_dmPFC_Lac] = corr(dmPFC_Lac(goodS.(raw_or_corr_nm).HME_f_dmPFC_Lac)', HME_ch(goodS.(raw_or_corr_nm).HME_f_dmPFC_Lac)');
    goodS.(raw_or_corr_nm).kEp_f_dmPFC_Lac = filter_fn(raw_or_corr_nm, dmPFC_Lac, kEp);
    [r_corr.(raw_or_corr_nm).kEp_f_dmPFC_Lac,...
        pval.(raw_or_corr_nm).kEp_f_dmPFC_Lac] = corr(dmPFC_Lac(goodS.(raw_or_corr_nm).kEp_f_dmPFC_Lac)', kEp(goodS.(raw_or_corr_nm).kEp_f_dmPFC_Lac)');
    
    % aIns lactate
    goodS.(raw_or_corr_nm).HE_f_aIns_Lac = filter_fn(raw_or_corr_nm, aIns_Lac, HE_ch);
    [r_corr.(raw_or_corr_nm).HE_f_aIns_Lac,...
        pval.(raw_or_corr_nm).HE_f_aIns_Lac] = corr(aIns_Lac(goodS.(raw_or_corr_nm).HE_f_aIns_Lac)', HE_ch(goodS.(raw_or_corr_nm).HE_f_aIns_Lac)');
    goodS.(raw_or_corr_nm).HPE_f_aIns_Lac = filter_fn(raw_or_corr_nm, aIns_Lac, HPE_ch);
    [r_corr.(raw_or_corr_nm).HPE_f_aIns_Lac,...
        pval.(raw_or_corr_nm).HPE_f_aIns_Lac] = corr(aIns_Lac(goodS.(raw_or_corr_nm).HPE_f_aIns_Lac)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_aIns_Lac)');
    goodS.(raw_or_corr_nm).HME_f_aIns_Lac = filter_fn(raw_or_corr_nm, aIns_Lac, HME_ch);
    [r_corr.(raw_or_corr_nm).HME_f_aIns_Lac,...
        pval.(raw_or_corr_nm).HME_f_aIns_Lac] = corr(aIns_Lac(goodS.(raw_or_corr_nm).HME_f_aIns_Lac)', HME_ch(goodS.(raw_or_corr_nm).HME_f_aIns_Lac)');
    goodS.(raw_or_corr_nm).kEp_f_aIns_Lac = filter_fn(raw_or_corr_nm, aIns_Lac, kEp);
    [r_corr.(raw_or_corr_nm).kEp_f_aIns_Lac,...
        pval.(raw_or_corr_nm).kEp_f_aIns_Lac] = corr(aIns_Lac(goodS.(raw_or_corr_nm).kEp_f_aIns_Lac)', kEp(goodS.(raw_or_corr_nm).kEp_f_aIns_Lac)');
    
    %% extract number of good subject for each correlation
    fields_to_check = fieldnames(goodS.(raw_or_corr_nm));
    for iF = 1:length(fields_to_check)
        field_nm = fields_to_check{iF};
        NS_goodS.(raw_or_corr_nm).(field_nm) = sum(goodS.(raw_or_corr_nm).(field_nm));
    end % loop over fields to check
    
    %% display correlation in a nice correlation matrix
    % assemble data in one correlation matrix
    plasma_r_vector = [r_corr.(raw_or_corr_nm).HE_f_plasma_Lac;...
        r_corr.(raw_or_corr_nm).HPE_f_plasma_Lac;...
        r_corr.(raw_or_corr_nm).HME_f_plasma_Lac;...
        r_corr.(raw_or_corr_nm).kEp_f_plasma_Lac];
    dmPFC_r_vector = [r_corr.(raw_or_corr_nm).HE_f_dmPFC_Lac;...
        r_corr.(raw_or_corr_nm).HPE_f_dmPFC_Lac;...
        r_corr.(raw_or_corr_nm).HME_f_dmPFC_Lac;...
        r_corr.(raw_or_corr_nm).kEp_f_dmPFC_Lac];
    aIns_r_vector = [r_corr.(raw_or_corr_nm).HE_f_aIns_Lac;...
        r_corr.(raw_or_corr_nm).HPE_f_aIns_Lac;...
        r_corr.(raw_or_corr_nm).HME_f_aIns_Lac;...
        r_corr.(raw_or_corr_nm).kEp_f_aIns_Lac];
    corr_mtrx = [plasma_r_vector, dmPFC_r_vector, aIns_r_vector];
    
    % same but for p.value
    plasma_pval_vector = [pval.(raw_or_corr_nm).HE_f_plasma_Lac;...
        pval.(raw_or_corr_nm).HPE_f_plasma_Lac;...
        pval.(raw_or_corr_nm).HME_f_plasma_Lac;...
        pval.(raw_or_corr_nm).kEp_f_plasma_Lac];
    dmPFC_pval_vector = [pval.(raw_or_corr_nm).HE_f_dmPFC_Lac;...
        pval.(raw_or_corr_nm).HPE_f_dmPFC_Lac;...
        pval.(raw_or_corr_nm).HME_f_dmPFC_Lac;...
        pval.(raw_or_corr_nm).kEp_f_dmPFC_Lac];
    aIns_pval_vector = [pval.(raw_or_corr_nm).HE_f_aIns_Lac;...
        pval.(raw_or_corr_nm).HPE_f_aIns_Lac;...
        pval.(raw_or_corr_nm).HME_f_aIns_Lac;...
        pval.(raw_or_corr_nm).kEp_f_aIns_Lac];
    pval_mtrx = [plasma_pval_vector, dmPFC_pval_vector, aIns_pval_vector];
    
    %% correlation between lactate and motivation
    fig;
    subplot_hdl = subplot(1,2,1);
    imagesc(corr_mtrx, corr_range);
    colormap(subplot_hdl, color_range_choices);
    cbar = colorbar;
    cbar.Label.String = 'r';
    xticks(1:3);
    xticklabels({'plasma','dmPFC/dACC','aIns'});
    yticks(1:size(corr_mtrx,1));
    yticklabels({'HE','HPE','HME','kEp'});
    % add stars in the graph if some correlations are significant
    for iROI = 1:size(corr_mtrx,2)
        for iBhv = 1:size(corr_mtrx,1)
            if pval_mtrx(iBhv, iROI) <= 0.05
                if pval_mtrx(iBhv, iROI) > 0.01 && pval_mtrx(iBhv, iROI) <= 0.05
                    pval_hdl = text(iROI, iBhv, '*');
                elseif pval_mtrx(iBhv, iROI) > 0.005 && pval_mtrx(iBhv, iROI) <= 0.01
                    pval_hdl = text(iROI, iBhv, '**');
                elseif pval_mtrx(iBhv, iROI) <= 0.005
                    pval_hdl = text(iROI, iBhv, '***');
                end % p.value
                % adjust p.value parameters
                pval_hdl.Color = col.white;
                pval_hdl.FontSize = 70;
                pval_hdl.FontWeight = 'bold';
                pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
                pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
            end % when p.value is significant
        end % loop over Y variables
    end % loop over X variables
    
end % loop over filter: raw vs 3*SD corrected data

end % function