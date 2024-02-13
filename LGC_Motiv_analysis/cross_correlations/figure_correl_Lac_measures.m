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
    %% first test: avoids redundant information and shows diagonal as one
    plasma_r_vector1 = [1;...
        r_corr.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        r_corr.(raw_or_corr_nm).aIns_f_plasma_Lac];
    dmPFC_r_vector1 = [0;...
        1;...
        r_corr.(raw_or_corr_nm).dmPFC_f_aIns_Lac];
    aIns_r_vector1 = [0;...
        0;...
        1];
    corr_mtrx1 = [plasma_r_vector1, dmPFC_r_vector1, aIns_r_vector1];
    
    % same but for p.value
    plasma_pval_vector1 = [NaN;...
        pval.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        pval.(raw_or_corr_nm).aIns_f_plasma_Lac];
    dmPFC_pval_vector1 = [NaN;...
        NaN;...
        pval.(raw_or_corr_nm).dmPFC_f_aIns_Lac];
    aIns_pval_vector1 = [NaN;...
        NaN;...
        NaN];
    pval_mtrx1 = [plasma_pval_vector1, dmPFC_pval_vector1, aIns_pval_vector1];
    
    
    %% second test shows symetry (redundant)
    plasma_r_vector2 = [1;...
        r_corr.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        r_corr.(raw_or_corr_nm).aIns_f_plasma_Lac];
    dmPFC_r_vector2 = [r_corr.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        1;...
        r_corr.(raw_or_corr_nm).dmPFC_f_aIns_Lac];
    aIns_r_vector2 = [r_corr.(raw_or_corr_nm).aIns_f_plasma_Lac;...
        r_corr.(raw_or_corr_nm).dmPFC_f_aIns_Lac;...
        1];
    corr_mtrx2 = [plasma_r_vector2, dmPFC_r_vector2, aIns_r_vector2];
    
    % same but for p.value
    plasma_pval_vector2 = [0;...
        pval.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        pval.(raw_or_corr_nm).aIns_f_plasma_Lac];
    dmPFC_pval_vector2 = [pval.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        0;...
        pval.(raw_or_corr_nm).dmPFC_f_aIns_Lac];
    aIns_pval_vector2 = [pval.(raw_or_corr_nm).aIns_f_plasma_Lac;...
        pval.(raw_or_corr_nm).dmPFC_f_aIns_Lac;...
        0];
    pval_mtrx2 = [plasma_pval_vector2, dmPFC_pval_vector2, aIns_pval_vector2];
    
    %% third test: like test 1 but without diagonal
    plasma_r_vector3 = [0;...
        r_corr.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        r_corr.(raw_or_corr_nm).aIns_f_plasma_Lac];
    dmPFC_r_vector3 = [0;...
        0;...
        r_corr.(raw_or_corr_nm).dmPFC_f_aIns_Lac];
    aIns_r_vector3 = [0;...
        0;...
        0];
    corr_mtrx3 = [plasma_r_vector3, dmPFC_r_vector3, aIns_r_vector3];
    
    % same but for p.value
    plasma_pval_vector3 = [NaN;...
        pval.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        pval.(raw_or_corr_nm).aIns_f_plasma_Lac];
    dmPFC_pval_vector3 = [NaN;...
        NaN;...
        pval.(raw_or_corr_nm).dmPFC_f_aIns_Lac];
    aIns_pval_vector3 = [NaN;...
        NaN;...
        NaN];
    pval_mtrx3 = [plasma_pval_vector3, dmPFC_pval_vector3, aIns_pval_vector3];
    
    %% fourth test: like test 2 but without diagonal
    plasma_r_vector4 = [0;...
        r_corr.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        r_corr.(raw_or_corr_nm).aIns_f_plasma_Lac];
    dmPFC_r_vector4 = [r_corr.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        0;...
        r_corr.(raw_or_corr_nm).dmPFC_f_aIns_Lac];
    aIns_r_vector4 = [r_corr.(raw_or_corr_nm).aIns_f_plasma_Lac;...
        r_corr.(raw_or_corr_nm).dmPFC_f_aIns_Lac;...
        0];
    corr_mtrx4 = [plasma_r_vector4, dmPFC_r_vector4, aIns_r_vector4];
    
    % same but for p.value
    plasma_pval_vector4 = [NaN;...
        pval.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        pval.(raw_or_corr_nm).aIns_f_plasma_Lac];
    dmPFC_pval_vector4 = [pval.(raw_or_corr_nm).dmPFC_f_plasma_Lac;...
        NaN;...
        pval.(raw_or_corr_nm).dmPFC_f_aIns_Lac];
    aIns_pval_vector4 = [pval.(raw_or_corr_nm).aIns_f_plasma_Lac;...
        pval.(raw_or_corr_nm).dmPFC_f_aIns_Lac;...
        NaN];
    pval_mtrx4 = [plasma_pval_vector4, dmPFC_pval_vector4, aIns_pval_vector4];
    
    %% figure 1
    fig;
    subplot_hdl = subplot(1,2,1);
    imagesc(corr_mtrx1, corr_range);
    colormap(subplot_hdl, color_range_choices);
    cbar = colorbar;
    cbar.Label.String = 'r';
    xticks(1:3);
    xticklabels({'plasma','dmPFC/dACC','aIns'});
    yticks(1:3);
    yticklabels({'plasma','dmPFC/dACC','aIns'});
    % add stars in the graph if some correlations are significant
    for iROI = 1:size(corr_mtrx1,2)
        for iBhv = 1:size(corr_mtrx1,1)
            if pval_mtrx1(iBhv, iROI) <= 0.05
                if pval_mtrx1(iBhv, iROI) > 0.01 && pval_mtrx1(iBhv, iROI) <= 0.05
                    pval_hdl = text(iROI, iBhv, '*');
                elseif pval_mtrx1(iBhv, iROI) > 0.005 && pval_mtrx1(iBhv, iROI) <= 0.01
                    pval_hdl = text(iROI, iBhv, '**');
                elseif pval_mtrx1(iBhv, iROI) <= 0.005
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
    
    %% figure 2
    fig;
    subplot_hdl = subplot(1,2,1);
    imagesc(corr_mtrx2, corr_range);
    colormap(subplot_hdl, color_range_choices);
    cbar = colorbar;
    cbar.Label.String = 'r';
    xticks(1:3);
    xticklabels({'plasma','dmPFC/dACC','aIns'});
    yticks(1:3);
    yticklabels({'plasma','dmPFC/dACC','aIns'});
    % add stars in the graph if some correlations are significant
    for iROI = 1:size(corr_mtrx2,2)
        for iBhv = 1:size(corr_mtrx2,1)
            if pval_mtrx2(iBhv, iROI) <= 0.05
                if pval_mtrx2(iBhv, iROI) > 0.01 && pval_mtrx2(iBhv, iROI) <= 0.05
                    pval_hdl = text(iROI, iBhv, '*');
                elseif pval_mtrx2(iBhv, iROI) > 0.005 && pval_mtrx2(iBhv, iROI) <= 0.01
                    pval_hdl = text(iROI, iBhv, '**');
                elseif pval_mtrx2(iBhv, iROI) <= 0.005
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
    
    %% figure 3
    fig;
    subplot_hdl = subplot(1,2,1);
    imagesc(corr_mtrx3, corr_range);
    colormap(subplot_hdl, color_range_choices);
    cbar = colorbar;
    cbar.Label.String = 'r';
    xticks(1:3);
    xticklabels({'plasma','dmPFC/dACC','aIns'});
    yticks(1:3);
    yticklabels({'plasma','dmPFC/dACC','aIns'});
    % add stars in the graph if some correlations are significant
    for iROI = 1:size(corr_mtrx3,2)
        for iBhv = 1:size(corr_mtrx3,1)
            if pval_mtrx3(iBhv, iROI) <= 0.05
                if pval_mtrx3(iBhv, iROI) > 0.01 && pval_mtrx3(iBhv, iROI) <= 0.05
                    pval_hdl = text(iROI, iBhv, '*');
                elseif pval_mtrx3(iBhv, iROI) > 0.005 && pval_mtrx3(iBhv, iROI) <= 0.01
                    pval_hdl = text(iROI, iBhv, '**');
                elseif pval_mtrx3(iBhv, iROI) <= 0.005
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
    
    %% figure 4
    fig;
    subplot_hdl = subplot(1,2,1);
    imagesc(corr_mtrx4, corr_range);
    colormap(subplot_hdl, color_range_choices);
    cbar = colorbar;
    cbar.Label.String = 'r';
    xticks(1:3);
    xticklabels({'plasma','dmPFC/dACC','aIns'});
    yticks(1:3);
    yticklabels({'plasma','dmPFC/dACC','aIns'});
    % add stars in the graph if some correlations are significant
    for iROI = 1:size(corr_mtrx4,2)
        for iBhv = 1:size(corr_mtrx4,1)
            if pval_mtrx4(iBhv, iROI) <= 0.05
                if pval_mtrx4(iBhv, iROI) > 0.01 && pval_mtrx4(iBhv, iROI) <= 0.05
                    pval_hdl = text(iROI, iBhv, '*');
                elseif pval_mtrx4(iBhv, iROI) > 0.005 && pval_mtrx4(iBhv, iROI) <= 0.01
                    pval_hdl = text(iROI, iBhv, '**');
                elseif pval_mtrx4(iBhv, iROI) <= 0.005
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