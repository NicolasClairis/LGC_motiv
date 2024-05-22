function[r_corr, pval, NS_goodS] = figure_interindiv_correl_BOLD__vs__motivation
% [r_corr, pval, NS_goodS] = figure_interindiv_correl_BOLD__vs__motivation
% figure_interindiv_correl_BOLD__vs__motivation will create a heatmap
% showing the correlation between high effort (HE), high physical effort
% (HPE), high mental effort (HME) choices and the sensitivity to physical
% effort (kEp) vs dmPFC/dACC and anterior insula (aIns) BOLD regression
% estimate for the effort chosen.
%
% OUTPUTS
% r_corr: structure with correlation coefficients for each correlation test
%
% pval: structure with p.values for each correlation test
%
% NS_goodS: number of good subjects included in each analysis

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;

%% extract ROI
GLM = spm_input('GLM number',1,'e');

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
% extract data for the contrast of interest
con_idx = listdlg('promptstring','Which contrast?',...
    'ListString',con_names);
[dmPFC_BOLD, aIns_BOLD] = deal(NaN(1,NS));
dmPFC_BOLD(:) = dmPFC_BOLD_allCons(con_idx, :, 1);
aIns_BOLD(:) = aIns_BOLD_allCons(con_idx, :, 1);

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

%% figure parameters
[~, ~, col] = general_fig_prm;
pSize = 21;
% define which colormap you want to use (see full list here if you are not
% happy with the selection:
% https://ch.mathworks.com/help/matlab/ref/colormap.html)
% color_range_choices = 'hot';
% color_range_choices = 'turbo';
% color_range_choices = 'jet';
color_range_choices = redblue(45);

% correlation range
corr_range = [-0.65 0.65];

% x/y-axis ratio
ax_ratio = [1.5 1 1];

% include kEm in the graph or not?
include_kEm = 1;

%% perform the correlations
for iRawCorr = 1:2
    switch iRawCorr
        case 1
            raw_or_corr_nm = 'raw';
        case 2
            raw_or_corr_nm = 'filtered';
    end
    
    % correlations with dmPFC/dACC BOLD
    % HE = f(dmPFC/dACC)
    goodS.(raw_or_corr_nm).HE_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, HE_ch);
    [r_corr.(raw_or_corr_nm).HE_f_dmPFC, pval.(raw_or_corr_nm).HE_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC)', HE_ch(goodS.(raw_or_corr_nm).HE_f_dmPFC)');
    [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).HE_f_dmPFC,...
        HE_ch_fit_dmPFC_sorted.(raw_or_corr_nm).HE_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC)', HE_ch(goodS.(raw_or_corr_nm).HE_f_dmPFC)','normal');
    
    % HPE = f(dmPFC/dACC)
    goodS.(raw_or_corr_nm).HPE_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, HPE_ch);
    [r_corr.(raw_or_corr_nm).HPE_f_dmPFC, pval.(raw_or_corr_nm).HPE_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_dmPFC)');
    [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).HPE_f_dmPFC,...
        HPE_ch_fit_dmPFC_sorted.(raw_or_corr_nm).HPE_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_dmPFC)','normal');
    
    % HME = f(dmPFC/dACC)
    goodS.(raw_or_corr_nm).HME_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, HME_ch);
    [r_corr.(raw_or_corr_nm).HME_f_dmPFC, pval.(raw_or_corr_nm).HME_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).HME_f_dmPFC)', HME_ch(goodS.(raw_or_corr_nm).HME_f_dmPFC)');
    [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).HME_f_dmPFC,...
        HME_ch_fit_dmPFC_sorted.(raw_or_corr_nm).HME_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).HME_f_dmPFC)', HME_ch(goodS.(raw_or_corr_nm).HME_f_dmPFC)','normal');
    
    % kEp = f(dmPFC/dACC)
    goodS.(raw_or_corr_nm).kEp_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, kEp);
    [r_corr.(raw_or_corr_nm).kEp_f_dmPFC, pval.(raw_or_corr_nm).kEp_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).kEp_f_dmPFC)', kEp(goodS.(raw_or_corr_nm).kEp_f_dmPFC)');
    [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).kEp_f_dmPFC,...
        kEp_fit_dmPFC_sorted.(raw_or_corr_nm).kEp_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).kEp_f_dmPFC)', kEp(goodS.(raw_or_corr_nm).kEp_f_dmPFC)','normal');
    
    % kEm = f(dmPFC/dACC)
    if include_kEm == 1
        goodS.(raw_or_corr_nm).kEm_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, kEm);
        [r_corr.(raw_or_corr_nm).kEm_f_dmPFC, pval.(raw_or_corr_nm).kEm_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).kEm_f_dmPFC)', kEm(goodS.(raw_or_corr_nm).kEm_f_dmPFC)');
        [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).kEm_f_dmPFC,...
            kEm_fit_dmPFC_sorted.(raw_or_corr_nm).kEm_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).kEm_f_dmPFC)', kEm(goodS.(raw_or_corr_nm).kEm_f_dmPFC)','normal');
    end
    
    % correlations with aIns BOLD
    % HE = f(aIns)
    goodS.(raw_or_corr_nm).HE_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, HE_ch);
    [r_corr.(raw_or_corr_nm).HE_f_aIns, pval.(raw_or_corr_nm).HE_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns)', HE_ch(goodS.(raw_or_corr_nm).HE_f_aIns)');
    [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).HE_f_aIns,...
        HE_ch_fit_aIns_sorted.(raw_or_corr_nm).HE_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns)', HE_ch(goodS.(raw_or_corr_nm).HE_f_aIns)','normal');
    
    % HPE = f(aIns)
    goodS.(raw_or_corr_nm).HPE_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, HPE_ch);
    [r_corr.(raw_or_corr_nm).HPE_f_aIns, pval.(raw_or_corr_nm).HPE_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_aIns)');
    [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).HPE_f_aIns,...
        HPE_ch_fit_aIns_sorted.(raw_or_corr_nm).HPE_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns)', HPE_ch(goodS.(raw_or_corr_nm).HPE_f_aIns)','normal');
    
    % HME = f(aIns)
    goodS.(raw_or_corr_nm).HME_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, HME_ch);
    [r_corr.(raw_or_corr_nm).HME_f_aIns, pval.(raw_or_corr_nm).HME_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).HME_f_aIns)', HME_ch(goodS.(raw_or_corr_nm).HME_f_aIns)');
    [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).HME_f_aIns,...
        HME_ch_fit_aIns_sorted.(raw_or_corr_nm).HME_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).HME_f_aIns)', HME_ch(goodS.(raw_or_corr_nm).HME_f_aIns)','normal');
    
    % kEp = f(aIns)
    goodS.(raw_or_corr_nm).kEp_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, kEp);
    [r_corr.(raw_or_corr_nm).kEp_f_aIns, pval.(raw_or_corr_nm).kEp_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).kEp_f_aIns)', kEp(goodS.(raw_or_corr_nm).kEp_f_aIns)');
    [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).kEp_f_aIns,...
        kEp_fit_aIns_sorted.(raw_or_corr_nm).kEp_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).kEp_f_aIns)', kEp(goodS.(raw_or_corr_nm).kEp_f_aIns)','normal');
    
    % kEm = f(aIns)
    if include_kEm == 1
        goodS.(raw_or_corr_nm).kEm_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, kEm);
        [r_corr.(raw_or_corr_nm).kEm_f_aIns, pval.(raw_or_corr_nm).kEm_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).kEm_f_aIns)', kEm(goodS.(raw_or_corr_nm).kEm_f_aIns)');
        [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).kEm_f_aIns,...
        kEm_fit_aIns_sorted.(raw_or_corr_nm).kEm_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).kEm_f_aIns)', kEm(goodS.(raw_or_corr_nm).kEm_f_aIns)','normal');
    end
    
    % also check the other behavioral parameters
    for iPrm = 1:nPrm
        prm_nm = prm_names{iPrm};
        % dmPFC/dACC
        dmPFC_prm_nm = [prm_nm,'_f_dmPFC'];
        goodS.(raw_or_corr_nm).(dmPFC_prm_nm) = filter_fn(raw_or_corr_nm, dmPFC_BOLD, prm.(prm_nm));
        [r_corr.(raw_or_corr_nm).prm.(dmPFC_prm_nm),...
            pval.(raw_or_corr_nm).prm.(dmPFC_prm_nm)] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).(dmPFC_prm_nm))',...
            prm.(prm_nm)(goodS.(raw_or_corr_nm).(dmPFC_prm_nm))');
        % same for aIns
        aIns_prm_nm = [prm_nm,'_f_aIns'];
        goodS.(raw_or_corr_nm).(aIns_prm_nm) = filter_fn(raw_or_corr_nm, aIns_BOLD, prm.(prm_nm));
        [r_corr.(raw_or_corr_nm).prm.(aIns_prm_nm),...
            pval.(raw_or_corr_nm).prm.(aIns_prm_nm)] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).(aIns_prm_nm))',...
            prm.(prm_nm)(goodS.(raw_or_corr_nm).(aIns_prm_nm))');
    end % prm loop
    %% extract number of good subject for each correlation
    fields_to_check = fieldnames(goodS.(raw_or_corr_nm));
    for iF = 1:length(fields_to_check)
        field_nm = fields_to_check{iF};
        NS_goodS.(raw_or_corr_nm).(field_nm) = sum(goodS.(raw_or_corr_nm).(field_nm));
    end % loop over fields to check
    
    %% display correlation in a nice correlation matrix
    
    switch include_kEm
        case 0 % correlation matrix HE/HPE/HME/kEp
            
            % assemble data in one correlation matrix
            dmPFC_r_vector = [r_corr.(raw_or_corr_nm).HE_f_dmPFC;...
                r_corr.(raw_or_corr_nm).HPE_f_dmPFC;...
                r_corr.(raw_or_corr_nm).HME_f_dmPFC;...
                r_corr.(raw_or_corr_nm).kEp_f_dmPFC];
            aIns_r_vector = [r_corr.(raw_or_corr_nm).HE_f_aIns;...
                r_corr.(raw_or_corr_nm).HPE_f_aIns;...
                r_corr.(raw_or_corr_nm).HME_f_aIns;...
                r_corr.(raw_or_corr_nm).kEp_f_aIns];
            corr_mtrx = [dmPFC_r_vector, aIns_r_vector];
            % same but for p.value
            dmPFC_pval_vector = [pval.(raw_or_corr_nm).HE_f_dmPFC;...
                pval.(raw_or_corr_nm).HPE_f_dmPFC;...
                pval.(raw_or_corr_nm).HME_f_dmPFC;...
                pval.(raw_or_corr_nm).kEp_f_dmPFC];
            aIns_pval_vector = [pval.(raw_or_corr_nm).HE_f_aIns;...
                pval.(raw_or_corr_nm).HPE_f_aIns;...
                pval.(raw_or_corr_nm).HME_f_aIns;...
                pval.(raw_or_corr_nm).kEp_f_aIns];
            pval_mtrx = [dmPFC_pval_vector, aIns_pval_vector];
            
            
        case 1 % correlation matrix HE/HPE/HME/kEp/kEm
            
            % assemble data in one correlation matrix
            dmPFC_r_vector = [r_corr.(raw_or_corr_nm).HE_f_dmPFC;...
                r_corr.(raw_or_corr_nm).HPE_f_dmPFC;...
                r_corr.(raw_or_corr_nm).HME_f_dmPFC;...
                r_corr.(raw_or_corr_nm).kEp_f_dmPFC;...
                r_corr.(raw_or_corr_nm).kEm_f_dmPFC];
            aIns_r_vector = [r_corr.(raw_or_corr_nm).HE_f_aIns;...
                r_corr.(raw_or_corr_nm).HPE_f_aIns;...
                r_corr.(raw_or_corr_nm).HME_f_aIns;...
                r_corr.(raw_or_corr_nm).kEp_f_aIns;...
                r_corr.(raw_or_corr_nm).kEm_f_aIns];
            corr_mtrx = [dmPFC_r_vector, aIns_r_vector];
            % same but for p.value
            dmPFC_pval_vector = [pval.(raw_or_corr_nm).HE_f_dmPFC;...
                pval.(raw_or_corr_nm).HPE_f_dmPFC;...
                pval.(raw_or_corr_nm).HME_f_dmPFC;...
                pval.(raw_or_corr_nm).kEp_f_dmPFC;...
                pval.(raw_or_corr_nm).kEm_f_dmPFC];
            aIns_pval_vector = [pval.(raw_or_corr_nm).HE_f_aIns;...
                pval.(raw_or_corr_nm).HPE_f_aIns;...
                pval.(raw_or_corr_nm).HME_f_aIns;...
                pval.(raw_or_corr_nm).kEp_f_aIns;...
                pval.(raw_or_corr_nm).kEm_f_aIns];
            pval_mtrx = [dmPFC_pval_vector, aIns_pval_vector];
    end
    %% figure
    switch include_kEm
        case 0 % correlation matrix HE/HPE/HME/kEp
            fig;
            subplot_hdl = subplot(1,2,1);
            imagesc(corr_mtrx, corr_range);
            colormap(subplot_hdl, color_range_choices);
            cbar = colorbar;
            cbar.Label.String = 'r';
            xticks(1:2);
            xticklabels({'dmPFC/dACC','aIns'});
            yticks(1:size(corr_mtrx,1));
            yticklabels({'HE','HPE','HME','kEp'});
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
            
            
        case 1 % correlation matrix HE/HPE/HME/kEp/kEm
            fig;
            subplot_hdl = subplot(1,2,1);
            imagesc(corr_mtrx, corr_range);
            colormap(subplot_hdl, color_range_choices);
            cbar = colorbar;
            cbar.Label.String = 'r';
            xticks(1:2);
            xticklabels({'dmPFC/dACC','aIns'});
            yticks(1:size(corr_mtrx,1));
            yticklabels({'HE','HPE','HME','kEp','kEm'});
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
            
    end % include kEm or not?
    
    %% display figures with correlation plots
    
    %% significant (HE, HPE, kEp) contrasts
    % dmPFC/dACC column figure
    fig;
    
    % HE = f(dmPFC/dACC)
    subplot(3,1,1); hold on;
    pbaspect(ax_ratio);
    HE_dmPFC_hdl = scatter(dmPFC_BOLD(goodS.(raw_or_corr_nm).HE_f_dmPFC)',...
        HE_ch(goodS.(raw_or_corr_nm).HE_f_dmPFC));
    HE_dmPFC_fit_hdl = plot(dmPFC_sorted.(raw_or_corr_nm).HE_f_dmPFC,...
        HE_ch_fit_dmPFC_sorted.(raw_or_corr_nm).HE_f_dmPFC);
    scat_hdl_upgrade(HE_dmPFC_hdl, col.grey);
    fit_hdl_upgrade(HE_dmPFC_fit_hdl, col.black);
    xlabel('dmPFC/dACC regression estimate');
    ylabel('HE choices (%)');
    legend_size(pSize);
    
    % HPE = f(dmPFC/dACC)
    subplot(3,1,2); hold on;
    pbaspect(ax_ratio);
    HPE_dmPFC_hdl = scatter(dmPFC_BOLD(goodS.(raw_or_corr_nm).HPE_f_dmPFC)',...
        HPE_ch(goodS.(raw_or_corr_nm).HPE_f_dmPFC));
    HPE_dmPFC_fit_hdl = plot(dmPFC_sorted.(raw_or_corr_nm).HPE_f_dmPFC,...
        HPE_ch_fit_dmPFC_sorted.(raw_or_corr_nm).HPE_f_dmPFC);
    scat_hdl_upgrade(HPE_dmPFC_hdl, col.grey);
    fit_hdl_upgrade(HPE_dmPFC_fit_hdl, col.black);
    xlabel('dmPFC/dACC regression estimate');
    ylabel('HPE choices (%)');
    legend_size(pSize);
    
    % kEp = f(dmPFC/dACC)
    subplot(3,1,3); hold on;
    pbaspect(ax_ratio);
    kEp_dmPFC_hdl = scatter(dmPFC_BOLD(goodS.(raw_or_corr_nm).kEp_f_dmPFC)',...
        kEp(goodS.(raw_or_corr_nm).kEp_f_dmPFC));
    kEp_dmPFC_fit_hdl = plot(dmPFC_sorted.(raw_or_corr_nm).kEp_f_dmPFC,...
        kEp_fit_dmPFC_sorted.(raw_or_corr_nm).kEp_f_dmPFC);
    scat_hdl_upgrade(kEp_dmPFC_hdl, col.grey);
    fit_hdl_upgrade(kEp_dmPFC_fit_hdl, col.black);
    xlabel('dmPFC/dACC regression estimate');
    ylabel('kEp');
    legend_size(pSize);
    
    
    % aIns column figure
    fig;
    % HE = f(aIns)
    subplot(3,1,1); hold on;
    pbaspect(ax_ratio);
    HE_aIns_hdl = scatter(aIns_BOLD(goodS.(raw_or_corr_nm).HE_f_aIns)',...
        HE_ch(goodS.(raw_or_corr_nm).HE_f_aIns));
    HE_aIns_fit_hdl = plot(aIns_sorted.(raw_or_corr_nm).HE_f_aIns,...
        HE_ch_fit_aIns_sorted.(raw_or_corr_nm).HE_f_aIns);
    scat_hdl_upgrade(HE_aIns_hdl, col.grey);
    fit_hdl_upgrade(HE_aIns_fit_hdl, col.black);
    xlabel('aIns regression estimate');
    ylabel('HE choices (%)');
    legend_size(pSize);
    
    % HPE = f(aIns)
    subplot(3,1,2); hold on;
    pbaspect(ax_ratio);
    HPE_aIns_hdl = scatter(aIns_BOLD(goodS.(raw_or_corr_nm).HPE_f_aIns)',...
        HPE_ch(goodS.(raw_or_corr_nm).HPE_f_aIns));
    HPE_aIns_fit_hdl = plot(aIns_sorted.(raw_or_corr_nm).HPE_f_aIns,...
        HPE_ch_fit_aIns_sorted.(raw_or_corr_nm).HPE_f_aIns);
    scat_hdl_upgrade(HPE_aIns_hdl, col.grey);
    fit_hdl_upgrade(HPE_aIns_fit_hdl, col.black);
    xlabel('aIns regression estimate');
    ylabel('HPE choices (%)');
    legend_size(pSize);
    
    
    
    % kEp = f(aIns)
    subplot(3,1,3); hold on;
    pbaspect(ax_ratio);
    kEp_aIns_hdl = scatter(aIns_BOLD(goodS.(raw_or_corr_nm).kEp_f_aIns)',...
        kEp(goodS.(raw_or_corr_nm).kEp_f_aIns));
    kEp_aIns_fit_hdl = plot(aIns_sorted.(raw_or_corr_nm).kEp_f_aIns,...
        kEp_fit_aIns_sorted.(raw_or_corr_nm).kEp_f_aIns);
    scat_hdl_upgrade(kEp_aIns_hdl, col.grey);
    fit_hdl_upgrade(kEp_aIns_fit_hdl, col.black);
    xlabel('aIns regression estimate');
    ylabel('kEp');
    legend_size(pSize);
    
    
    %% non-significant (HME, kEm) contrasts
    % dmPFC/dACC column figure
    fig;
    
    % HME = f(dmPFC/dACC)
    subplot(3,1,1); hold on;
    pbaspect(ax_ratio);
    HME_dmPFC_hdl = scatter(dmPFC_BOLD(goodS.(raw_or_corr_nm).HME_f_dmPFC)',...
        HME_ch(goodS.(raw_or_corr_nm).HME_f_dmPFC));
    HME_dmPFC_fit_hdl = plot(dmPFC_sorted.(raw_or_corr_nm).HME_f_dmPFC,...
        HME_ch_fit_dmPFC_sorted.(raw_or_corr_nm).HME_f_dmPFC);
    scat_hdl_upgrade(HME_dmPFC_hdl, col.grey);
    fit_hdl_upgrade(HME_dmPFC_fit_hdl, col.black);
    xlabel('dmPFC/dACC regression estimate');
    ylabel('HME choices (%)');
    legend_size(pSize);
    
    % kEm = f(dmPFC/dACC)
    subplot(3,1,2); hold on;
    pbaspect(ax_ratio);
    kEm_dmPFC_hdl = scatter(dmPFC_BOLD(goodS.(raw_or_corr_nm).kEm_f_dmPFC)',...
        kEm(goodS.(raw_or_corr_nm).kEm_f_dmPFC));
    kEm_dmPFC_fit_hdl = plot(dmPFC_sorted.(raw_or_corr_nm).kEm_f_dmPFC,...
        kEm_fit_dmPFC_sorted.(raw_or_corr_nm).kEm_f_dmPFC);
    scat_hdl_upgrade(kEm_dmPFC_hdl, col.grey);
    fit_hdl_upgrade(kEm_dmPFC_fit_hdl, col.black);
    xlabel('dmPFC/dACC regression estimate');
    ylabel('kEm');
    legend_size(pSize);
    
    
    % aIns column figure
    fig;
    
    % HME = f(aIns)
    subplot(3,1,1); hold on;
    pbaspect(ax_ratio);
    HME_aIns_hdl = scatter(aIns_BOLD(goodS.(raw_or_corr_nm).HME_f_aIns)',...
        HME_ch(goodS.(raw_or_corr_nm).HME_f_aIns));
    HME_aIns_fit_hdl = plot(aIns_sorted.(raw_or_corr_nm).HME_f_aIns,...
        HME_ch_fit_aIns_sorted.(raw_or_corr_nm).HME_f_aIns);
    scat_hdl_upgrade(HME_aIns_hdl, col.grey);
    fit_hdl_upgrade(HME_aIns_fit_hdl, col.black);
    xlabel('aIns regression estimate');
    ylabel('HME choices (%)');
    legend_size(pSize);
    
    % kEm = f(aIns)
    subplot(3,1,2); hold on;
    pbaspect(ax_ratio);
    kEm_aIns_hdl = scatter(aIns_BOLD(goodS.(raw_or_corr_nm).kEm_f_aIns)',...
        kEm(goodS.(raw_or_corr_nm).kEm_f_aIns));
    kEm_aIns_fit_hdl = plot(aIns_sorted.(raw_or_corr_nm).kEm_f_aIns,...
        kEm_fit_aIns_sorted.(raw_or_corr_nm).kEm_f_aIns);
    scat_hdl_upgrade(kEm_aIns_hdl, col.grey);
    fit_hdl_upgrade(kEm_aIns_fit_hdl, col.black);
    xlabel('aIns regression estimate');
    ylabel('kEm');
    legend_size(pSize);
    
end % loop over filter: raw vs 3*SD corrected data

end % function