function[r_corr, pval, NS_goodS] = figure_correl_Lac_measures()
% [r_corr, pval] = figure_correl_Lac_measures
% figure_correl_Lac_measures will create a heatmap
% showing the correlation between plasma, dmPFC/dACC and anterior 
% insula (aIns) lactate levels.
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
    
    % perform lactate correlations
    % plasma/dmPFC lactate
    goodS.(raw_or_corr_nm).dmPFC_f_plasma_Lac = filter_fn(raw_or_corr_nm, plasma_Lac, dmPFC_Lac);
    [r_corr.(raw_or_corr_nm).dmPFC_f_plasma_Lac,...
        pval.(raw_or_corr_nm).dmPFC_f_plasma_Lac] = corr(plasma_Lac(goodS.(raw_or_corr_nm).dmPFC_f_plasma_Lac)', dmPFC_Lac(goodS.(raw_or_corr_nm).dmPFC_f_plasma_Lac)');
    % plasma/aIns lactate
    goodS.(raw_or_corr_nm).aIns_f_plasma_Lac = filter_fn(raw_or_corr_nm, plasma_Lac, aIns_Lac);
    [r_corr.(raw_or_corr_nm).aIns_f_plasma_Lac,...
        pval.(raw_or_corr_nm).aIns_f_plasma_Lac] = corr(plasma_Lac(goodS.(raw_or_corr_nm).aIns_f_plasma_Lac)', aIns_Lac(goodS.(raw_or_corr_nm).aIns_f_plasma_Lac)');
    % dmPFC/aIns lactate
    goodS.(raw_or_corr_nm).dmPFC_f_aIns_Lac = filter_fn(raw_or_corr_nm, aIns_Lac, dmPFC_Lac);
    [r_corr.(raw_or_corr_nm).dmPFC_f_aIns_Lac,...
        pval.(raw_or_corr_nm).dmPFC_f_aIns_Lac] = corr(aIns_Lac(goodS.(raw_or_corr_nm).dmPFC_f_aIns_Lac)', dmPFC_Lac(goodS.(raw_or_corr_nm).dmPFC_f_aIns_Lac)');
    
    %% extract number of good subject for each correlation
    fields_to_check = fieldnames(goodS.(raw_or_corr_nm));
    for iF = 1:length(fields_to_check)
        field_nm = fields_to_check{iF};
        NS_goodS.(raw_or_corr_nm).(field_nm) = sum(goodS.(raw_or_corr_nm).(field_nm));
    end % loop over fields to check
    
    %% display correlation in a nice correlation matrix
    % assemble data in one correlation matrix:
    % avoid redundant information and shows diagonal as one
    dmPFC_r_vector = [1;...
        r_corr.(raw_or_corr_nm).dmPFC_f_aIns_Lac;...
        r_corr.(raw_or_corr_nm).dmPFC_f_plasma_Lac];
    aIns_r_vector = [0;...
        1;...
        r_corr.(raw_or_corr_nm).aIns_f_plasma_Lac];
    plasma_r_vector = [0;...
        0;...
        1];
    corr_mtrx = [dmPFC_r_vector, aIns_r_vector, plasma_r_vector];
    
    % same but for p.value
    dmPFC_pval_vector = [NaN;...
        pval.(raw_or_corr_nm).dmPFC_f_aIns_Lac;...
        pval.(raw_or_corr_nm).dmPFC_f_plasma_Lac];
    aIns_pval_vector = [NaN;...
        NaN;...
        pval.(raw_or_corr_nm).aIns_f_plasma_Lac];
    plasma_pval_vector = [NaN;...
        NaN;...
        NaN];
    pval_mtrx = [dmPFC_pval_vector, aIns_pval_vector, plasma_pval_vector];
    
    %% figure
    fig;
    subplot_hdl = subplot(1,2,1);
    imagesc(corr_mtrx, corr_range);
    colormap(subplot_hdl, color_range_choices);
    cbar = colorbar;
    cbar.Label.String = 'r';
    xticks(1:3);
    xticklabels({'dmPFC/dACC','aIns','plasma'});
    yticks(1:3);
    yticklabels({'dmPFC/dACC','aIns','plasma'});
    % add stars in the graph if some correlations are significant
    for iROI = 1:size(corr_mtrx,2)
        for iBhv = 1:size(corr_mtrx,1)
            if pval_mtrx(iBhv, iROI) <= 0.05
                if pval_mtrx(iBhv, iROI) > 0.01 && pval_mtrx(iBhv, iROI) <= 0.05
                    pval_hdl = text(iROI, iBhv, '*');
                elseif pval_mtrx(iBhv, iROI) > 0.001 && pval_mtrx(iBhv, iROI) <= 0.01
                    pval_hdl = text(iROI, iBhv, '**');
                elseif pval_mtrx(iBhv, iROI) <= 0.001
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