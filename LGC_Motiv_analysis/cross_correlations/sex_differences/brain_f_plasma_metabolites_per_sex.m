% correlation between plasma and brain metabolites

%% subject identification
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

%% load plasma metabolites
[plasmaM_males, plasma_mb_names, n_plasma_mb] = load_plasma_metabolites(male_CIDS);
[plasmaM_females, ~, ~] = load_plasma_metabolites(female_CIDS);

%% load brain metabolites
[metabolites_males] = metabolite_load(male_CIDS);
[metabolites_females] = metabolite_load(female_CIDS);
brain_mb_names = fieldnames(metabolites_males.dmPFC);
ignore_metabolites_ratios = {'Scyllo','Glu_div_GSH','Glu_div_Tau','Glu_div_antiox','Glu_div_z_antiox','Glu_div_GABA_div_GSH','Glu_div_GABA','Asp_plus_Lac',...
    'zAsp_plus_zLac','Glu_Ushape','GSH_Ushape','Glu_Ushape_div_GSH','Glu_div_GSH_Ushape','Glu_Ushape_div_GSH_Ushape','Glu_div_GSH_ratio_Ushape',...
    'Glu_div_Asp','Lac_div_GSH','Lac_div_antiox','Glu_plus_Lac_div_GSH','Glu_plus_Lac_div_antiox'};
brain_mb_names(ismember(brain_mb_names, ignore_metabolites_ratios)) = [];
n_brain_mb = length(brain_mb_names);

%% define r threshold for display
rThresh = 0.2;

%% loop through metabolites to extract all correlations
[dmPFC_corr_males_mtrx, dmPFC_pval_males_mtrx,...
    aIns_corr_males_mtrx, aIns_pval_males_mtrx,...
    dmPFC_corr_females_mtrx, dmPFC_pval_females_mtrx,...
    aIns_corr_females_mtrx, aIns_pval_females_mtrx] = deal(NaN(n_brain_mb, n_plasma_mb));
plasma_mb_names_short = plasma_mb_names;
brain_mb_names_short = brain_mb_names;
for iPlasma = 1:n_plasma_mb
    plasma_mb_nm = plasma_mb_names{iPlasma};
    plasma_mb_names_short{iPlasma} = plasma_mb_names_short{iPlasma}(1:3); % shorten name for final display
    for iBrain = 1:n_brain_mb
        brain_mb_nm = brain_mb_names{iBrain};
        % shorten name for final display (if name is long)
        if strcmp(brain_mb_names{iBrain},'Gln_div_Glu')
            brain_mb_names_short{iBrain} = 'Gln/Glu';
        elseif length(brain_mb_names{iBrain}) > 3
            brain_mb_names_short{iBrain} = brain_mb_names{iBrain}(1:3);
        end
        
        %% extract relevant data
        % males
        plasma_mb_males_tmp = plasmaM_males.(plasma_mb_nm);
        dmPFC_mb_males_tmp = metabolites_males.dmPFC.(brain_mb_nm);
        aIns_mb_males_tmp = metabolites_males.aIns.(brain_mb_nm);
        % females
        plasma_mb_females_tmp = plasmaM_females.(plasma_mb_nm);
        dmPFC_mb_females_tmp = metabolites_females.dmPFC.(brain_mb_nm);
        aIns_mb_females_tmp = metabolites_females.aIns.(brain_mb_nm);
        
        %% check number of subjects
        NS_goodS.(plasma_mb_nm).dmPFC.(brain_mb_nm).males = sum(~isnan(plasma_mb_males_tmp.*dmPFC_mb_males_tmp'));
        NS_goodS.(plasma_mb_nm).aIns.(brain_mb_nm).females = sum(~isnan(plasma_mb_males_tmp.*aIns_mb_males_tmp'));
        
        %% test significance
        % correlation dmPFC/dACC in males
        [r_corr.(plasma_mb_nm).dmPFC.(brain_mb_nm).males,...
            betas.(plasma_mb_nm).dmPFC.(brain_mb_nm).males,...
            pval.(plasma_mb_nm).dmPFC.(brain_mb_nm).males,...
            ~,...
            plasma.(plasma_mb_nm).sorted.males,...
            dmPFC.([brain_mb_nm,'_f_',plasma_mb_nm]).males] = glm_package(plasma_mb_males_tmp', dmPFC_mb_males_tmp, 'normal', 'on');
        % correlation dmPFC/dACC in females
        [r_corr.(plasma_mb_nm).dmPFC.(brain_mb_nm).females,...
            betas.(plasma_mb_nm).dmPFC.(brain_mb_nm).females,...
            pval.(plasma_mb_nm).dmPFC.(brain_mb_nm).females,...
            ~,...
            plasma.(plasma_mb_nm).sorted.females,...
            dmPFC.([brain_mb_nm,'_f_',plasma_mb_nm]).females] = glm_package(plasma_mb_females_tmp', dmPFC_mb_females_tmp, 'normal', 'on');
        
        % correlation aIns in males
        [r_corr.(plasma_mb_nm).aIns.(brain_mb_nm).males,...
            betas.(plasma_mb_nm).aIns.(brain_mb_nm).males,...
            pval.(plasma_mb_nm).aIns.(brain_mb_nm).males,...
            ~,...
            ~,...
            aIns.([brain_mb_nm,'_f_',plasma_mb_nm]).males] = glm_package(plasma_mb_males_tmp', aIns_mb_males_tmp, 'normal', 'on');
        % correlation aIns in females
        [r_corr.(plasma_mb_nm).aIns.(brain_mb_nm).females,...
            betas.(plasma_mb_nm).aIns.(brain_mb_nm).females,...
            pval.(plasma_mb_nm).aIns.(brain_mb_nm).females,...
            ~,...
            ~,...
            aIns.([brain_mb_nm,'_f_',plasma_mb_nm]).females] = glm_package(plasma_mb_females_tmp', aIns_mb_females_tmp, 'normal', 'on');
        
        % correlation dmPFC/dACC <=> aIns in males
        [r_corr.(plasma_mb_nm).dmPFC_vs_aIns.(brain_mb_nm).males,...
            betas.(plasma_mb_nm).dmPFC_vs_aIns.(brain_mb_nm).males,...
            pval.(plasma_mb_nm).dmPFC_vs_aIns.(brain_mb_nm).males,...
            ~,...
            plasma.(plasma_mb_nm).sorted.males,...
            dmPFC_vs_aIns.([brain_mb_nm,'_f_',plasma_mb_nm]).males] = glm_package(plasma_mb_males_tmp', dmPFC_mb_males_tmp-aIns_mb_males_tmp, 'normal', 'on');
        [r_corr.(plasma_mb_nm).dmPFC_vs_aIns.(brain_mb_nm).females,...
            betas.(plasma_mb_nm).dmPFC_vs_aIns.(brain_mb_nm).females,...
            pval.(plasma_mb_nm).dmPFC_vs_aIns.(brain_mb_nm).females,...
            ~,...
            plasma.(plasma_mb_nm).sorted.females,...
            dmPFC_vs_aIns.([brain_mb_nm,'_f_',plasma_mb_nm]).females] = glm_package(plasma_mb_females_tmp', dmPFC_mb_females_tmp-aIns_mb_females_tmp, 'normal', 'on');
        
        % load data in big matrix for correlation plot
        % dmPFC/dACC
        % males
        if abs(r_corr.(plasma_mb_nm).dmPFC.(brain_mb_nm).males) > rThresh
            dmPFC_corr_males_mtrx(iBrain, iPlasma) = r_corr.(plasma_mb_nm).dmPFC.(brain_mb_nm).males;
        else
            dmPFC_corr_males_mtrx(iBrain, iPlasma) = 0;
        end
        dmPFC_pval_males_mtrx(iBrain, iPlasma) = pval.(plasma_mb_nm).dmPFC.(brain_mb_nm).males(2);
        % females
        if abs(r_corr.(plasma_mb_nm).dmPFC.(brain_mb_nm).females) > rThresh
            dmPFC_corr_females_mtrx(iBrain, iPlasma) = r_corr.(plasma_mb_nm).dmPFC.(brain_mb_nm).females;
        else
            dmPFC_corr_females_mtrx(iBrain, iPlasma) = 0;
        end
        dmPFC_pval_females_mtrx(iBrain, iPlasma) = pval.(plasma_mb_nm).dmPFC.(brain_mb_nm).females(2);
        
        % anterior insula
        % males
        if abs(r_corr.(plasma_mb_nm).aIns.(brain_mb_nm).males) > rThresh
            aIns_corr_males_mtrx(iBrain, iPlasma) = r_corr.(plasma_mb_nm).aIns.(brain_mb_nm).males;
        else
            aIns_corr_males_mtrx(iBrain, iPlasma) = 0;
        end
        aIns_pval_males_mtrx(iBrain, iPlasma) = pval.(plasma_mb_nm).aIns.(brain_mb_nm).males(2);
        % females
        if abs(r_corr.(plasma_mb_nm).aIns.(brain_mb_nm).females) > rThresh
            aIns_corr_females_mtrx(iBrain, iPlasma) = r_corr.(plasma_mb_nm).aIns.(brain_mb_nm).females;
        else
            aIns_corr_females_mtrx(iBrain, iPlasma) = 0;
        end
        aIns_pval_females_mtrx(iBrain, iPlasma) = pval.(plasma_mb_nm).aIns.(brain_mb_nm).females(2);
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

% color for each sex
[col_m1, col_f1, col_m2, col_f2] = col_per_sex;

%% plot correlation matrix for dmPFC/dACC and aIns - MALES
fig;

% dmPFC/dACC
dmPFC_males_subplot_hdl = subplot(1,2,1);
imagesc(dmPFC_corr_males_mtrx, corr_range);
colormap(dmPFC_males_subplot_hdl, color_range_choices);
cbar = colorbar;
cbar.Label.String = 'r';
xticks(1:n_plasma_mb);
xticklabels(plasma_mb_names_short);
xlabel('plasma metabolites');
yticks(1:n_brain_mb);
yticklabels(brain_mb_names_short);
ylabel('dmPFC/dACC metabolites');
% add stars in the graph if some correlations are significant
for iPlasma = 1:size(dmPFC_corr_males_mtrx,2)
    for idmPFC_mb = 1:size(dmPFC_corr_males_mtrx,1)
        if dmPFC_pval_males_mtrx(idmPFC_mb, iPlasma) <= 0.05
            if dmPFC_pval_males_mtrx(idmPFC_mb, iPlasma) > 0.01 && dmPFC_pval_males_mtrx(idmPFC_mb, iPlasma) <= 0.05
                pval_hdl = text(iPlasma, idmPFC_mb, '*');
            elseif dmPFC_pval_males_mtrx(idmPFC_mb, iPlasma) > 0.001 && dmPFC_pval_males_mtrx(idmPFC_mb, iPlasma) <= 0.01
                pval_hdl = text(iPlasma, idmPFC_mb, '**');
            elseif dmPFC_pval_males_mtrx(idmPFC_mb, iPlasma) <= 0.001
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
legend_size(15);

% anterior insula
aIns_males_subplot_hdl = subplot(1,2,2);
imagesc(aIns_corr_males_mtrx, corr_range);
colormap(aIns_males_subplot_hdl, color_range_choices);
cbar = colorbar;
cbar.Label.String = 'r';
xticks(1:n_plasma_mb);
xticklabels(plasma_mb_names_short);
xlabel('plasma metabolites');
yticks(1:n_brain_mb);
yticklabels(brain_mb_names_short);
ylabel('aIns metabolites');
% add stars in the graph if some correlations are significant
for iPlasma = 1:size(aIns_corr_males_mtrx,2)
    for iaIns_mb = 1:size(aIns_corr_males_mtrx,1)
        if aIns_pval_males_mtrx(iaIns_mb, iPlasma) <= 0.05
            if aIns_pval_males_mtrx(iaIns_mb, iPlasma) > 0.01 && aIns_pval_males_mtrx(iaIns_mb, iPlasma) <= 0.05
                pval_hdl = text(iPlasma, iaIns_mb, '*');
            elseif aIns_pval_males_mtrx(iaIns_mb, iPlasma) > 0.001 && aIns_pval_males_mtrx(iaIns_mb, iPlasma) <= 0.01
                pval_hdl = text(iPlasma, iaIns_mb, '**');
            elseif aIns_pval_males_mtrx(iaIns_mb, iPlasma) <= 0.001
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

%% plot correlation matrix for dmPFC/dACC and aIns - FEMALES
fig;

% dmPFC/dACC
dmPFC_females_subplot_hdl = subplot(1,2,1);
imagesc(dmPFC_corr_females_mtrx, corr_range);
colormap(dmPFC_females_subplot_hdl, color_range_choices);
cbar = colorbar;
cbar.Label.String = 'r';
xticks(1:n_plasma_mb);
xticklabels(plasma_mb_names_short);
xlabel('plasma metabolites');
yticks(1:n_brain_mb);
yticklabels(brain_mb_names_short);
ylabel('dmPFC/dACC metabolites');
% add stars in the graph if some correlations are significant
for iPlasma = 1:size(dmPFC_corr_females_mtrx,2)
    for idmPFC_mb = 1:size(dmPFC_corr_females_mtrx,1)
        if dmPFC_pval_females_mtrx(idmPFC_mb, iPlasma) <= 0.05
            if dmPFC_pval_females_mtrx(idmPFC_mb, iPlasma) > 0.01 && dmPFC_pval_females_mtrx(idmPFC_mb, iPlasma) <= 0.05
                pval_hdl = text(iPlasma, idmPFC_mb, '*');
            elseif dmPFC_pval_females_mtrx(idmPFC_mb, iPlasma) > 0.001 && dmPFC_pval_females_mtrx(idmPFC_mb, iPlasma) <= 0.01
                pval_hdl = text(iPlasma, idmPFC_mb, '**');
            elseif dmPFC_pval_females_mtrx(idmPFC_mb, iPlasma) <= 0.001
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
legend_size(15);

% anterior insula
aIns_females_subplot_hdl = subplot(1,2,2);
imagesc(aIns_corr_females_mtrx, corr_range);
colormap(aIns_females_subplot_hdl, color_range_choices);
cbar = colorbar;
cbar.Label.String = 'r';
xticks(1:n_plasma_mb);
xticklabels(plasma_mb_names_short);
xlabel('plasma metabolites');
yticks(1:n_brain_mb);
yticklabels(brain_mb_names_short);
ylabel('aIns metabolites');
% add stars in the graph if some correlations are significant
for iPlasma = 1:size(aIns_corr_females_mtrx,2)
    for iaIns_mb = 1:size(aIns_corr_females_mtrx,1)
        if aIns_pval_females_mtrx(iaIns_mb, iPlasma) <= 0.05
            if aIns_pval_females_mtrx(iaIns_mb, iPlasma) > 0.01 && aIns_pval_females_mtrx(iaIns_mb, iPlasma) <= 0.05
                pval_hdl = text(iPlasma, iaIns_mb, '*');
            elseif aIns_pval_females_mtrx(iaIns_mb, iPlasma) > 0.001 && aIns_pval_females_mtrx(iaIns_mb, iPlasma) <= 0.01
                pval_hdl = text(iPlasma, iaIns_mb, '**');
            elseif aIns_pval_females_mtrx(iaIns_mb, iPlasma) <= 0.001
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

%% display correlations for some specific metabolites between periphery and brain
metabolites_to_disp = {'Gln','Glu'};
n_mb_to_disp = length(metabolites_to_disp);

for iMb = 1:n_mb_to_disp
    mb_nm = metabolites_to_disp{iMb};

    fig;
    % dmPFC/dACC
    subplot(1,2,1); hold on;
    scat_dmPFC_males_hdl = scatter(plasmaM_males.(mb_nm), metabolites_males.dmPFC.(mb_nm));
    scat_dmPFC_females_hdl = scatter(plasmaM_females.(mb_nm), metabolites_females.dmPFC.(mb_nm));
    scat_hdl_upgrade(scat_dmPFC_males_hdl);
    scat_hdl_upgrade(scat_dmPFC_females_hdl);
    scat_dmPFC_males_hdl.MarkerFaceColor = col_m1;
    scat_dmPFC_females_hdl.MarkerFaceColor = col_f1;
    % add fit
    dmPFC_males_fit_hdl = plot(plasma.(mb_nm).sorted.males, dmPFC.([mb_nm,'_f_',mb_nm]).males);
    dmPFC_females_fit_hdl = plot(plasma.(mb_nm).sorted.females, dmPFC.([mb_nm,'_f_',mb_nm]).females);
    fit_hdl_upgrade(dmPFC_males_fit_hdl);
    fit_hdl_upgrade(dmPFC_females_fit_hdl);
    dmPFC_males_fit_hdl.Color = col_m2;
    dmPFC_females_fit_hdl.Color = col_f2;
    % add correlation
    place_r_and_pval_2(r_corr.(mb_nm).dmPFC.(mb_nm).males,...
        pval.(mb_nm).dmPFC.(mb_nm).males(2),...
        col_m2,...
        r_corr.(mb_nm).dmPFC.(mb_nm).females,...
        pval.(mb_nm).dmPFC.(mb_nm).females(2),...
        col_f2);
    xlabel(['Plasma ',mb_nm,' (μM)']);
    ylabel(['dmPFC/dACC ',mb_nm,' (mM)']);

    % aIns
    subplot(1,2,2); hold on;
    scat_aIns_males_hdl = scatter(plasmaM_males.(mb_nm), metabolites_males.aIns.(mb_nm));
    scat_aIns_females_hdl = scatter(plasmaM_females.(mb_nm), metabolites_females.aIns.(mb_nm));
    scat_hdl_upgrade(scat_aIns_males_hdl);
    scat_hdl_upgrade(scat_aIns_females_hdl);
    scat_aIns_males_hdl.MarkerFaceColor = col_m1;
    scat_aIns_females_hdl.MarkerFaceColor = col_f1;
    % add fit
    aIns_males_fit_hdl = plot(plasma.(mb_nm).sorted.males, aIns.([mb_nm,'_f_',mb_nm]).males);
    aIns_females_fit_hdl = plot(plasma.(mb_nm).sorted.females, aIns.([mb_nm,'_f_',mb_nm]).females);
    fit_hdl_upgrade(aIns_males_fit_hdl);
    fit_hdl_upgrade(aIns_females_fit_hdl);
    aIns_males_fit_hdl.Color = col_m2;
    aIns_females_fit_hdl.Color = col_f2;
    % add correlation
    place_r_and_pval_2(r_corr.(mb_nm).aIns.(mb_nm).males,...
        pval.(mb_nm).aIns.(mb_nm).males(2),...
        col_m2,...
        r_corr.(mb_nm).aIns.(mb_nm).females,...
        pval.(mb_nm).aIns.(mb_nm).females(2),...
        col_f2);
    xlabel(['Plasma ',mb_nm,' (μM)']);
    ylabel(['aIns ',mb_nm,' (mM)']);
end % metabolite loop
