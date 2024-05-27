% correlation between plasma and brain metabolites

%% subject identification
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load plasma metabolites
[plasmaM, plasma_mb_names, n_plasma_mb] = load_plasma_metabolites(subject_id);

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
brain_mb_names = fieldnames(metabolites.dmPFC);
ignore_metabolites_ratios = {'Scyllo','Glu_div_GSH','Glu_div_Tau','Glu_div_antiox','Glu_div_z_antiox','Glu_div_GABA_div_GSH','Glu_div_GABA','Asp_plus_Lac',...
    'zAsp_plus_zLac','Glu_Ushape','GSH_Ushape','Glu_Ushape_div_GSH','Glu_div_GSH_Ushape','Glu_Ushape_div_GSH_Ushape','Glu_div_GSH_ratio_Ushape',...
    'Glu_div_Asp','Lac_div_GSH','Lac_div_antiox','Glu_plus_Lac_div_GSH','Glu_plus_Lac_div_antiox'};
brain_mb_names(ismember(brain_mb_names, ignore_metabolites_ratios)) = [];
n_brain_mb = length(brain_mb_names);

%% loop through metabolites to extract all correlations
[dmPFC_corr_mtrx, dmPFC_pval_mtrx,...
    aIns_corr_mtrx, aIns_pval_mtrx] = deal(NaN(n_brain_mb, n_plasma_mb));
plasma_mb_names_short = plasma_mb_names;
brain_mb_names_short = brain_mb_names;
for iPlasma = 1:n_plasma_mb
    plasma_mb_nm = plasma_mb_names{iPlasma};
    plasma_mb_names_short{iPlasma} = plasma_mb_names_short{iPlasma}(1:3); % shorten name for final display
    for iBrain = 1:n_brain_mb
        brain_mb_nm = brain_mb_names{iBrain};
        % shorten name for final display (if name is long)
        if length(brain_mb_names{iBrain}) > 3
            brain_mb_names_short{iBrain} = brain_mb_names{iBrain}(1:3);
        end
        
        %% extract relevant data
        plasma_mb_tmp = plasmaM.(plasma_mb_nm);
        dmPFC_mb_tmp = metabolites.dmPFC.(brain_mb_nm);
        aIns_mb_tmp = metabolites.aIns.(brain_mb_nm);
        
        %% check number of subjects
        NS_goodS.(plasma_mb_nm).dmPFC.(brain_mb_nm) = sum(~isnan(plasma_mb_tmp.*dmPFC_mb_tmp'));
        NS_goodS.(plasma_mb_nm).aIns.(brain_mb_nm) = sum(~isnan(plasma_mb_tmp.*aIns_mb_tmp'));
        
        %% test significance
        [r_corr.(plasma_mb_nm).dmPFC.(brain_mb_nm),...
            betas.(plasma_mb_nm).dmPFC.(brain_mb_nm),...
            pval.(plasma_mb_nm).dmPFC.(brain_mb_nm)] = glm_package(plasma_mb_tmp', dmPFC_mb_tmp, 'normal', 'on');
        [r_corr.(plasma_mb_nm).aIns.(brain_mb_nm),...
            betas.(plasma_mb_nm).aIns.(brain_mb_nm),...
            pval.(plasma_mb_nm).aIns.(brain_mb_nm)] = glm_package(plasma_mb_tmp', aIns_mb_tmp, 'normal', 'on');
        [r_corr.(plasma_mb_nm).dmPFC_vs_aIns.(brain_mb_nm),...
            betas.(plasma_mb_nm).dmPFC_vs_aIns.(brain_mb_nm),...
            pval.(plasma_mb_nm).dmPFC_vs_aIns.(brain_mb_nm)] = glm_package(plasma_mb_tmp', dmPFC_mb_tmp-aIns_mb_tmp, 'normal', 'on');
        
        % load data in big matrix for correlation plot
        % dmPFC/dACC
        if abs(r_corr.(plasma_mb_nm).dmPFC.(brain_mb_nm)) > 0.2
            dmPFC_corr_mtrx(iBrain, iPlasma) = r_corr.(plasma_mb_nm).dmPFC.(brain_mb_nm);
        else
            dmPFC_corr_mtrx(iBrain, iPlasma) = 0;
        end
        dmPFC_pval_mtrx(iBrain, iPlasma) = pval.(plasma_mb_nm).dmPFC.(brain_mb_nm)(2);
        
        % anterior insula
        if abs(r_corr.(plasma_mb_nm).aIns.(brain_mb_nm)) > 0.2
            aIns_corr_mtrx(iBrain, iPlasma) = r_corr.(plasma_mb_nm).aIns.(brain_mb_nm);
        else
            aIns_corr_mtrx(iBrain, iPlasma) = 0;
        end
        aIns_pval_mtrx(iBrain, iPlasma) = pval.(plasma_mb_nm).aIns.(brain_mb_nm)(2);
    end % brain metabolite loop
end % plasma metabolite loop

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
corr_range = [-1 1];

%% plot correlation matrix for dmPFC/dACC and aIns
fig;

% dmPFC/dACC
dmPFC_subplot_hdl = subplot(1,2,1);
imagesc(dmPFC_corr_mtrx, corr_range);
colormap(dmPFC_subplot_hdl, color_range_choices);
cbar = colorbar;
cbar.Label.String = 'r';
xticks(1:n_plasma_mb);
xticklabels(plasma_mb_names_short);
xlabel('plasma metabolites');
yticks(1:n_brain_mb);
yticklabels(brain_mb_names_short);
ylabel('dmPFC/dACC metabolites');
% add stars in the graph if some correlations are significant
for iPlasma = 1:size(dmPFC_corr_mtrx,2)
    for idmPFC_mb = 1:size(dmPFC_corr_mtrx,1)
        if dmPFC_pval_mtrx(idmPFC_mb, iPlasma) <= 0.05
            if dmPFC_pval_mtrx(idmPFC_mb, iPlasma) > 0.01 && dmPFC_pval_mtrx(idmPFC_mb, iPlasma) <= 0.05
                pval_hdl = text(iPlasma, idmPFC_mb, '*');
            elseif dmPFC_pval_mtrx(idmPFC_mb, iPlasma) > 0.001 && dmPFC_pval_mtrx(idmPFC_mb, iPlasma) <= 0.01
                pval_hdl = text(iPlasma, idmPFC_mb, '**');
            elseif dmPFC_pval_mtrx(idmPFC_mb, iPlasma) <= 0.001
                pval_hdl = text(iPlasma, idmPFC_mb, '***');
            end % p.value
            % adjust p.value parameters
            pval_hdl.Color = col.black;
            pval_hdl.FontSize = 10;
            pval_hdl.FontWeight = 'bold';
            pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
            pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
        end % when p.value is significant
    end % loop over Y variables
end % loop over X variables

% anterior insula
aIns_subplot_hdl = subplot(1,2,2);
imagesc(aIns_corr_mtrx, corr_range);
colormap(aIns_subplot_hdl, color_range_choices);
cbar = colorbar;
cbar.Label.String = 'r';
xticks(1:n_plasma_mb);
xticklabels(plasma_mb_names_short);
xlabel('plasma metabolites');
yticks(1:n_brain_mb);
yticklabels(brain_mb_names_short);
ylabel('aIns metabolites');
% add stars in the graph if some correlations are significant
for iPlasma = 1:size(aIns_corr_mtrx,2)
    for iaIns_mb = 1:size(aIns_corr_mtrx,1)
        if aIns_pval_mtrx(iaIns_mb, iPlasma) <= 0.05
            if aIns_pval_mtrx(iaIns_mb, iPlasma) > 0.01 && aIns_pval_mtrx(iaIns_mb, iPlasma) <= 0.05
                pval_hdl = text(iPlasma, iaIns_mb, '*');
            elseif aIns_pval_mtrx(iaIns_mb, iPlasma) > 0.001 && aIns_pval_mtrx(iaIns_mb, iPlasma) <= 0.01
                pval_hdl = text(iPlasma, iaIns_mb, '**');
            elseif aIns_pval_mtrx(iaIns_mb, iPlasma) <= 0.001
                pval_hdl = text(iPlasma, iaIns_mb, '***');
            end % p.value
            % adjust p.value parameters
            pval_hdl.Color = col.black;
            pval_hdl.FontSize = 10;
            pval_hdl.FontWeight = 'bold';
            pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
            pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
        end % when p.value is significant
    end % loop over Y variables
end % loop over X variables
legend_size(15);