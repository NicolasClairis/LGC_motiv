function[r_corr, pval, NS_goodS] = figure_interindiv_correl_BOLD__vs__motivation_bis
% [r_corr, pval, NS_goodS] = figure_interindiv_correl_BOLD__vs__motivation_bis
% figure_interindiv_correl_BOLD__vs__motivation_bis will create a heatmap
% showing the correlation between all the behavioral parameters not displayed
% in figure_interindiv_correl_BOLD__vs__motivation.m (kFp, kLm, kR, kP, kBias)
% and the dmPFC/dACC and anterior insula (aIns) BOLD regression
% estimate selected.
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
[prm] = prm_extraction(study_nm, subject_id);
prm_names = fieldnames(prm);
prm_names(strcmp(prm_names,'CID')) = [];
nPrm = length(prm_names);
% create variables of interest
kR = prm.kR;
kP = prm.kP;
kFp = prm.kFp;
kLm = prm.kLm;
kBias = prm.kBias;

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

%% perform the correlations
for iRawCorr = 1:2
    switch iRawCorr
        case 1
            raw_or_corr_nm = 'raw';
        case 2
            raw_or_corr_nm = 'filtered';
    end
    
    % correlations with dmPFC/dACC BOLD
    % kFp = f(dmPFC/dACC)
    goodS.(raw_or_corr_nm).kFp_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, kFp);
    [r_corr.(raw_or_corr_nm).kFp_f_dmPFC, pval.(raw_or_corr_nm).kFp_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).kFp_f_dmPFC)', kFp(goodS.(raw_or_corr_nm).kFp_f_dmPFC)');
    [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).kFp_f_dmPFC,...
        kFp_fit_dmPFC_sorted.(raw_or_corr_nm).kFp_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).kFp_f_dmPFC)', kFp(goodS.(raw_or_corr_nm).kFp_f_dmPFC)','normal');
    
    % kLm = f(dmPFC/dACC)
    goodS.(raw_or_corr_nm).kLm_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, kLm);
    [r_corr.(raw_or_corr_nm).kLm_f_dmPFC, pval.(raw_or_corr_nm).kLm_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).kLm_f_dmPFC)', kLm(goodS.(raw_or_corr_nm).kLm_f_dmPFC)');
    [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).kLm_f_dmPFC,...
        kLm_fit_dmPFC_sorted.(raw_or_corr_nm).kLm_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).kLm_f_dmPFC)', kLm(goodS.(raw_or_corr_nm).kLm_f_dmPFC)','normal');
    
    % kR = f(dmPFC/dACC)
    goodS.(raw_or_corr_nm).kR_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, kR);
    [r_corr.(raw_or_corr_nm).kR_f_dmPFC, pval.(raw_or_corr_nm).kR_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).kR_f_dmPFC)', kR(goodS.(raw_or_corr_nm).kR_f_dmPFC)');
    [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).kR_f_dmPFC,...
        kR_fit_dmPFC_sorted.(raw_or_corr_nm).kR_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).kR_f_dmPFC)', kR(goodS.(raw_or_corr_nm).kR_f_dmPFC)','normal');
    
    % kP = f(dmPFC/dACC)
    goodS.(raw_or_corr_nm).kP_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, kP);
    [r_corr.(raw_or_corr_nm).kP_f_dmPFC, pval.(raw_or_corr_nm).kP_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).kP_f_dmPFC)', kP(goodS.(raw_or_corr_nm).kP_f_dmPFC)');
    [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).kP_f_dmPFC,...
        kP_fit_dmPFC_sorted.(raw_or_corr_nm).kP_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).kP_f_dmPFC)', kP(goodS.(raw_or_corr_nm).kP_f_dmPFC)','normal');
    
    % kBias = f(dmPFC/dACC)
    goodS.(raw_or_corr_nm).kBias_f_dmPFC = filter_fn(raw_or_corr_nm, dmPFC_BOLD, kBias);
    [r_corr.(raw_or_corr_nm).kBias_f_dmPFC, pval.(raw_or_corr_nm).kBias_f_dmPFC] = corr(dmPFC_BOLD(goodS.(raw_or_corr_nm).kBias_f_dmPFC)', kBias(goodS.(raw_or_corr_nm).kBias_f_dmPFC)');
    [~, ~, ~, ~, dmPFC_sorted.(raw_or_corr_nm).kBias_f_dmPFC,...
        kBias_fit_dmPFC_sorted.(raw_or_corr_nm).kBias_f_dmPFC] = glm_package(dmPFC_BOLD(goodS.(raw_or_corr_nm).kBias_f_dmPFC)', kBias(goodS.(raw_or_corr_nm).kBias_f_dmPFC)','normal');
    
    % correlations with aIns BOLD
    % kFp = f(aIns)
    goodS.(raw_or_corr_nm).kFp_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, kFp);
    [r_corr.(raw_or_corr_nm).kFp_f_aIns, pval.(raw_or_corr_nm).kFp_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).kFp_f_aIns)', kFp(goodS.(raw_or_corr_nm).kFp_f_aIns)');
    [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).kFp_f_aIns,...
        kFp_fit_aIns_sorted.(raw_or_corr_nm).kFp_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).kFp_f_aIns)', kFp(goodS.(raw_or_corr_nm).kFp_f_aIns)','normal');
    
    % kLm = f(aIns)
    goodS.(raw_or_corr_nm).kLm_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, kLm);
    [r_corr.(raw_or_corr_nm).kLm_f_aIns, pval.(raw_or_corr_nm).kLm_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).kLm_f_aIns)', kLm(goodS.(raw_or_corr_nm).kLm_f_aIns)');
    [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).kLm_f_aIns,...
        kLm_fit_aIns_sorted.(raw_or_corr_nm).kLm_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).kLm_f_aIns)', kLm(goodS.(raw_or_corr_nm).kLm_f_aIns)','normal');
    
    % kR = f(aIns)
    goodS.(raw_or_corr_nm).kR_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, kR);
    [r_corr.(raw_or_corr_nm).kR_f_aIns, pval.(raw_or_corr_nm).kR_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).kR_f_aIns)', kR(goodS.(raw_or_corr_nm).kR_f_aIns)');
    [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).kR_f_aIns,...
        kR_fit_aIns_sorted.(raw_or_corr_nm).kR_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).kR_f_aIns)', kR(goodS.(raw_or_corr_nm).kR_f_aIns)','normal');
    
    % kP = f(aIns)
    goodS.(raw_or_corr_nm).kP_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, kP);
    [r_corr.(raw_or_corr_nm).kP_f_aIns, pval.(raw_or_corr_nm).kP_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).kP_f_aIns)', kP(goodS.(raw_or_corr_nm).kP_f_aIns)');
    [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).kP_f_aIns,...
        kP_fit_aIns_sorted.(raw_or_corr_nm).kP_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).kP_f_aIns)', kP(goodS.(raw_or_corr_nm).kP_f_aIns)','normal');
    
    % kBias = f(aIns)
    goodS.(raw_or_corr_nm).kBias_f_aIns = filter_fn(raw_or_corr_nm, aIns_BOLD, kBias);
    [r_corr.(raw_or_corr_nm).kBias_f_aIns, pval.(raw_or_corr_nm).kBias_f_aIns] = corr(aIns_BOLD(goodS.(raw_or_corr_nm).kBias_f_aIns)', kBias(goodS.(raw_or_corr_nm).kBias_f_aIns)');
    [~, ~, ~, ~, aIns_sorted.(raw_or_corr_nm).kBias_f_aIns,...
        kBias_fit_aIns_sorted.(raw_or_corr_nm).kBias_f_aIns] = glm_package(aIns_BOLD(goodS.(raw_or_corr_nm).kBias_f_aIns)', kBias(goodS.(raw_or_corr_nm).kBias_f_aIns)','normal');
    
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
    % assemble data in one correlation matrix
    dmPFC_r_vector = [r_corr.(raw_or_corr_nm).kR_f_dmPFC;...
        r_corr.(raw_or_corr_nm).kP_f_dmPFC;...
        r_corr.(raw_or_corr_nm).kFp_f_dmPFC;...
        r_corr.(raw_or_corr_nm).kLm_f_dmPFC;...
        r_corr.(raw_or_corr_nm).kBias_f_dmPFC];
    aIns_r_vector = [r_corr.(raw_or_corr_nm).kR_f_aIns;...
        r_corr.(raw_or_corr_nm).kP_f_aIns;...
        r_corr.(raw_or_corr_nm).kFp_f_aIns;...
        r_corr.(raw_or_corr_nm).kLm_f_aIns;...
        r_corr.(raw_or_corr_nm).kBias_f_aIns];
    corr_mtrx = [dmPFC_r_vector, aIns_r_vector];
    % same but for p.value
    dmPFC_pval_vector = [pval.(raw_or_corr_nm).kR_f_dmPFC;...
        pval.(raw_or_corr_nm).kP_f_dmPFC;...
        pval.(raw_or_corr_nm).kFp_f_dmPFC;...
        pval.(raw_or_corr_nm).kLm_f_dmPFC;...
        pval.(raw_or_corr_nm).kBias_f_dmPFC];
    aIns_pval_vector = [pval.(raw_or_corr_nm).kR_f_aIns;...
        pval.(raw_or_corr_nm).kP_f_aIns;...
        pval.(raw_or_corr_nm).kFp_f_aIns;...
        pval.(raw_or_corr_nm).kLm_f_aIns;...
        pval.(raw_or_corr_nm).kBias_f_aIns];
    pval_mtrx = [dmPFC_pval_vector, aIns_pval_vector];
    %% figure
    fig;
    subplot_hdl = subplot(1,2,1);
    imagesc(corr_mtrx, corr_range);
    colormap(subplot_hdl, color_range_choices);
    cbar = colorbar;
    cbar.Label.String = 'r';
    xticks(1:2);
    xticklabels({'dmPFC/dACC','aIns'});
    yticks(1:size(corr_mtrx,1));
    yticklabels({'kR','kP','kFp','kLm','kBias'});
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
    
%     %% display figures with correlation plots
%     % dmPFC/dACC column figure
%     fig;
%     
%     % kFp = f(dmPFC/dACC)
%     subplot(3,1,3); hold on;
%     pbaspect(ax_ratio);
%     kFp_dmPFC_hdl = scatter(dmPFC_BOLD(goodS.(raw_or_corr_nm).kFp_f_dmPFC)',...
%         kFp(goodS.(raw_or_corr_nm).kFp_f_dmPFC));
%     kFp_dmPFC_fit_hdl = plot(dmPFC_sorted.(raw_or_corr_nm).kFp_f_dmPFC,...
%         kFp_fit_dmPFC_sorted.(raw_or_corr_nm).kFp_f_dmPFC);
%     scat_hdl_upgrade(kFp_dmPFC_hdl, col.grey);
%     fit_hdl_upgrade(kFp_dmPFC_fit_hdl, col.black);
%     xlabel('dmPFC/dACC regression estimate');
%     ylabel('kFp');
%     legend_size(pSize);
%     
%     
%     % aIns column figure
%     fig;
%     
%     % kFp = f(aIns)
%     subplot(3,1,3); hold on;
%     pbaspect(ax_ratio);
%     kFp_aIns_hdl = scatter(aIns_BOLD(goodS.(raw_or_corr_nm).kFp_f_aIns)',...
%         kFp(goodS.(raw_or_corr_nm).kFp_f_aIns));
%     kFp_aIns_fit_hdl = plot(aIns_sorted.(raw_or_corr_nm).kFp_f_aIns,...
%         kFp_fit_aIns_sorted.(raw_or_corr_nm).kFp_f_aIns);
%     scat_hdl_upgrade(kFp_aIns_hdl, col.grey);
%     fit_hdl_upgrade(kFp_aIns_fit_hdl, col.black);
%     xlabel('aIns regression estimate');
%     ylabel('kFp');
%     legend_size(pSize);
    
end % loop over filter: raw vs 3*SD corrected data

end % function