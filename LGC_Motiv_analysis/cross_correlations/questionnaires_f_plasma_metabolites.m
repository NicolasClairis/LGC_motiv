% function[rho, pval, signif,...
%     corr_mtrx_plasma, corr_mtrx_dmPFC, corr_mtrx_aIns,...
%     pval_unc_plasma, pval_unc_dmPFC, pval_unc_aIns]= questionnaires_f_plasma_metabolites(fig_disp, rmv_outliers_yn)
% [rho, pval, signif,...
%     corr_mtrx_plasma, corr_mtrx_dmPFC, corr_mtrx_aIns,...
%     pval_unc_plasma, pval_unc_dmPFC, pval_unc_aIns]= questionnaires_f_plasma_metabolites(fig_disp, rmv_outliers_yn)
% questionnaires_f_plasma_metabolites will test whether the scores from the questionnaires
% correlate with lactate measures in the plasma or brain,
% after grouping them by category (stress/motivation/etc.).
%
% INPUTS
% fig_disp: display figures? (1) yes (0) no
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% rho: structure with correlation coefficients
% 
% pval: structure with p.values
% 
% signif: structure with significant results
%
% corr_mtrx_plasma: nQuestionnaires*nPlasma matrix of correlation
% coefficients
% 
% corr_mtrx_dmPFC: nQuestionnaires*nBrainMetabolites matrix of correlation
% coefficients for dmPFC/dACC
% 
% corr_mtrx_aIns: nQuestionnaires*nBrainMetabolites matrix of correlation
% coefficients for anterior insula
% 
% pval_unc_plasma: nQuestionnaires*nPlasma matrix of p.values for plasma vs
% questionnaires correlations
% 
% pval_unc_dmPFC: nQuestionnaires*nPlasma matrix of p.values for dmPFC vs
% questionnaires correlations
% 
% pval_unc_aIns: nQuestionnaires*nPlasma matrix of p.values for aIns vs
% questionnaires correlations

%% define inputs by default
% display figure by default
if ~exist('fig_disp','var') || isempty(fig_disp) || ~ismember(fig_disp,[0,1])
    fig_disp = 1;
end
% remove outliers by default
if ~exist('rmv_outliers_yn','var') || isempty(rmv_outliers_yn) || ~ismember(rmv_outliers_yn,[0,1])
    rmv_outliers_yn = 1;
end

%% subject selection
[study_nm, condition, gender, subject_id, NS] = sub_id;

%% extract lactate measures
% plasma metabolites
[plasmaM, plasmaM_names, n_plasmaM] = load_plasma_metabolites(subject_id);
% change values from μM to mM
for iPlasmaM = 1:n_plasmaM
    plasmaM_nm = plasmaM_names{iPlasmaM};
    plasmaM.(plasmaM_nm) = plasmaM.(plasmaM_nm)./1000;
end

% brain metabolites
[~,~,brainMetabolites] = metabolite_load(subject_id);
dmPFC_mb = brainMetabolites.dmPFC;
aIns_mb = brainMetabolites.aIns;
brainM_names = fieldnames(dmPFC_mb);
n_brainM = length(brainM_names);

%% extract questionnaires
[questionnaires, categ_quests, n_categ] = extract_questionnaires(study_nm, subject_id, NS);

%% remove outliers
if rmv_outliers_yn == 1
    % plasma metabolites
    for iPlasmaM = 1:n_plasmaM
        plasmaM_nm = plasmaM_names{iPlasmaM};
        [~,~,plasmaM.(plasmaM_nm)] = rmv_outliers_3sd(plasmaM.(plasmaM_nm));
    end

    % brain metabolites
    for iBrainM = 1:n_brainM
        brainM_nm = brainM_names{iBrainM};
        [~,~,dmPFC_mb.(brainM_nm)] = rmv_outliers_3sd(dmPFC_mb.(brainM_nm));
        [~,~,aIns_mb.(brainM_nm)] = rmv_outliers_3sd(aIns_mb.(brainM_nm));
    end
    
    % questionnaires
    for iC = 1:n_categ
        categ_nm = categ_quests{iC};
        
        quest_names = fieldnames(questionnaires.(categ_nm));
        n_Q = length(quest_names);
        for iQ = 1:n_Q
            quest_nm = quest_names{iQ};
            [~,~,questionnaires.(categ_nm).(quest_nm)] = rmv_outliers_3sd(questionnaires.(categ_nm).(quest_nm));
        end % questionnaire loop
    end % questionnaire categories
end % outliers

%% figure parameters
if fig_disp == 1
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
    corr_range = [-1 1];

    % x/y-axis ratio
    ax_ratio = [1.5 1 1];
end

%% perform correlations for each questionnaire
for iCateg = 1:n_categ
    categ_nm = categ_quests{iCateg};
    
    % category subfields
    quest_names = fieldnames(questionnaires.(categ_nm));
    n_quests = length(quest_names);
    
    [corr_mtrx_plasma, pval_unc_plasma] = deal(NaN(n_quests, n_plasmaM));
    [corr_mtrx_dmPFC, corr_mtrx_aIns,...
        pval_unc_dmPFC, pval_unc_aIns] = deal(NaN(n_quests, n_brainM));
    for iQ = 1:n_quests
        quest_nm = quest_names{iQ};
        
        % plasma metabolites
        for iPM = 1:n_plasmaM
            plasmaM_nm = plasmaM_names{iPlasmaM};
            goodS = ~isnan(plasmaM.(plasmaM_nm).*questionnaires.(categ_nm).(quest_nm));
            % perform the correlation
            [rho.(categ_nm).(quest_nm).plasma.(plasmaM_nm), pval.uncorrected.(categ_nm).(quest_nm).plasma.(plasmaM_nm)] = corr(...
                plasmaM.(plasmaM_nm)(goodS)',...
                questionnaires.(categ_nm).(quest_nm)(goodS)');
            % extract fitted data
            [~,~,~,~,plasmaM_sorted.(plasmaM_nm).(categ_nm).(quest_nm),...
                questionnairess_fit_sorted_f_plasma.(plasmaM_nm).(categ_nm).(quest_nm)] = glm_package(plasmaM.(plasmaM_nm)', questionnaires.(categ_nm).(quest_nm)','normal');
            % extract the data in a big matrix for later display
            corr_mtrx_plasma(iQ,iPM) = rho.(categ_nm).(quest_nm).plasma.(plasmaM_nm);
            pval_unc_plasma(iQ,iPM) = pval.uncorrected.(categ_nm).(quest_nm).plasma.(plasmaM_nm);
            % store significant results in a different subfield
            if pval.uncorrected.(categ_nm).(quest_nm).plasma.(plasmaM_nm) < 0.05
                % store coefficient
                signif.unc.(categ_nm).(quest_nm).plasma.(plasmaM_nm).r_corr = rho.(categ_nm).(quest_nm).plasma.(plasmaM_nm);
                % store p.value
                signif.unc.(categ_nm).(quest_nm).plasma.(plasmaM_nm).pval = pval.uncorrected.(categ_nm).(quest_nm).plasma.(plasmaM_nm);
            end
        end % plasma loop

        % brain metabolites
        for iBM = 1:n_brainM
            brainM_nm = brainM_names{iBM};
            
            % dmPFC
            goodS = ~isnan(dmPFC_mb.(brainM_nm).*questionnaires.(categ_nm).(quest_nm));
            % perform the correlation
            [rho.(categ_nm).(quest_nm).dmPFC.(brainM_nm), pval.uncorrected.(categ_nm).(quest_nm).dmPFC.(brainM_nm)] = corr(...
                dmPFC_mb.(brainM_nm)(goodS)',...
                questionnaires.(categ_nm).(quest_nm)(goodS)');
            % extract fitted data
            [~,~,~,~,dmPFC_sorted.(brainM_nm).(categ_nm).(quest_nm),...
                questionnairess_fit_sorted_f_dmPFC.(brainM_nm).(categ_nm).(quest_nm)] = glm_package(dmPFC_mb.(brainM_nm)', questionnaires.(categ_nm).(quest_nm)','normal');
            % extract the data in a big matrix for later display
            corr_mtrx_dmPFC(iQ,iBM) = rho.(categ_nm).(quest_nm).dmPFC.(brainM_nm);
            pval_unc_dmPFC(iQ,iBM) = pval.uncorrected.(categ_nm).(quest_nm).dmPFC.(brainM_nm);
            % store significant results in a different subfield
            if pval.uncorrected.(categ_nm).(quest_nm).dmPFC.(brainM_nm) < 0.05
                % store coefficient
                signif.unc.(categ_nm).(quest_nm).dmPFC.(brainM_nm).r_corr = rho.(categ_nm).(quest_nm).dmPFC.(brainM_nm);
                % store p.value
                signif.unc.(categ_nm).(quest_nm).dmPFC.(brainM_nm).pval = pval.uncorrected.(categ_nm).(quest_nm).dmPFC.(brainM_nm);
            end

            % aIns
            goodS = ~isnan(aIns_mb.(brainM_nm).*questionnaires.(categ_nm).(quest_nm));
            % perform the correlation
            [rho.(categ_nm).(quest_nm).aIns.(brainM_nm), pval.uncorrected.(categ_nm).(quest_nm).aIns.(brainM_nm)] = corr(...
                aIns_mb.(brainM_nm)(goodS)',...
                questionnaires.(categ_nm).(quest_nm)(goodS)');
            % extract fitted data
            [~,~,~,~,aIns_sorted.(brainM_nm).(categ_nm).(quest_nm),...
                questionnairess_fit_sorted_f_aIns.(brainM_nm).(categ_nm).(quest_nm)] = glm_package(aIns_mb.(brainM_nm)', questionnaires.(categ_nm).(quest_nm)','normal');
            % extract the data in a big matrix for later display
            corr_mtrx_aIns(iQ,iBM) = rho.(categ_nm).(quest_nm).aIns.(brainM_nm);
            pval_unc_aIns(iQ,iBM) = pval.uncorrected.(categ_nm).(quest_nm).aIns.(brainM_nm);
            % store significant results in a different subfield
            if pval.uncorrected.(categ_nm).(quest_nm).aIns.(brainM_nm) < 0.05
                % store coefficient
                signif.unc.(categ_nm).(quest_nm).aIns.(brainM_nm).r_corr = rho.(categ_nm).(quest_nm).aIns.(brainM_nm);
                % store p.value
                signif.unc.(categ_nm).(quest_nm).aIns.(brainM_nm).pval = pval.uncorrected.(categ_nm).(quest_nm).aIns.(brainM_nm);
            end
        end % brain metabolites
    end % loop over questionnaires
    
    %% show results for current in a graph (+ if results are significant or not)
    if fig_disp == 1
        % general figure parameters
        [pSize, lW, col, mSize] = general_fig_prm;

        %% display correlation matrices

        %% plasma
        plasma_corr_fig = fig;
        imagesc(corr_mtrx_plasma, corr_range);
        colormap(plasma_corr_fig, color_range_choices);
        cbar = colorbar;
        cbar.Label.String = 'r';
        xticks(1:n_plasmaM);
        xticklabels(plasmaM_names);
        xlabel('Plasma metabolites');
        yticks(1:n_quests);
        yticklabels(quest_names);
        ylabel('Questionnaires');
        % add stars in the graph if some correlations are significant
        for iPM = 1:n_plasmaM
            for iQ = 1:n_quests
                if pval_unc_plasma(iQ, iPM) <= 0.05
                    if pval_unc_plasma(iQ, iPM) > 0.01 && pval_unc_plasma(iQ, iPM) <= 0.05
                        pval_hdl = text(iPM, iQ, '*');
                    elseif pval_unc_plasma(iQ, iPM) > 0.001 && pval_unc_plasma(iQ, iPM) <= 0.01
                        pval_hdl = text(iPM, iQ, '**');
                    elseif pval_unc_plasma(iQ, iPM) <= 0.001
                        pval_hdl = text(iPM, iQ, '***');
                    end % p.value
                    % adjust p.value parameters
                    pval_hdl.Color = col.white;
                    pval_hdl.FontSize = 70;
                    pval_hdl.FontWeight = 'bold';
                    pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
                    pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
                end % when p.value is significant
            end % loop over questionnaires
        end % loop over plasma metabolites

        %% dmPFC/dACC
        dmPFC_corr_fig = fig;
        imagesc(corr_mtrx_dmPFC, corr_range);
        colormap(dmPFC_corr_fig, color_range_choices);
        cbar = colorbar;
        cbar.Label.String = 'r';
        xticks(1:n_brainM);
        xticklabels(brainM_names);
        xlabel('dmPFC/dACC metabolites');
        yticks(1:n_quests);
        yticklabels(quest_names);
        ylabel('Questionnaires');
        % add stars in the graph if some correlations are significant
        for i_dmPFC_mb = 1:n_brainM
            for iQ = 1:n_quests
                if pval_unc_dmPFC(iQ, i_dmPFC_mb) <= 0.05
                    if pval_unc_dmPFC(iQ, i_dmPFC_mb) > 0.01 && pval_unc_dmPFC(iQ, i_dmPFC_mb) <= 0.05
                        pval_hdl = text(i_dmPFC_mb, iQ, '*');
                    elseif pval_unc_dmPFC(iQ, i_dmPFC_mb) > 0.001 && pval_unc_dmPFC(iQ, i_dmPFC_mb) <= 0.01
                        pval_hdl = text(i_dmPFC_mb, iQ, '**');
                    elseif pval_unc_dmPFC(iQ, i_dmPFC_mb) <= 0.001
                        pval_hdl = text(i_dmPFC_mb, iQ, '***');
                    end % p.value
                    % adjust p.value parameters
                    pval_hdl.Color = col.white;
                    pval_hdl.FontSize = 70;
                    pval_hdl.FontWeight = 'bold';
                    pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
                    pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
                end % when p.value is significant
            end % loop over questionnaires
        end % loop over dmPFC metabolites

        %% aIns
        aIns_corr_fig = fig;
        imagesc(corr_mtrx_aIns, corr_range);
        colormap(aIns_corr_fig, color_range_choices);
        cbar = colorbar;
        cbar.Label.String = 'r';
        xticks(1:n_brainM);
        xticklabels(brainM_names);
        xlabel('aIns metabolites');
        yticks(1:n_quests);
        yticklabels(quest_names);
        ylabel('Questionnaires');
        % add stars in the graph if some correlations are significant
        for i_aIns_mb = 1:n_brainM
            for iQ = 1:n_quests
                if pval_unc_aIns(iQ, i_aIns_mb) <= 0.05
                    if pval_unc_aIns(iQ, i_aIns_mb) > 0.01 && pval_unc_aIns(iQ, i_aIns_mb) <= 0.05
                        pval_hdl = text(i_aIns_mb, iQ, '*');
                    elseif pval_unc_aIns(iQ, i_aIns_mb) > 0.001 && pval_unc_aIns(iQ, i_aIns_mb) <= 0.01
                        pval_hdl = text(i_aIns_mb, iQ, '**');
                    elseif pval_unc_aIns(iQ, i_aIns_mb) <= 0.001
                        pval_hdl = text(i_aIns_mb, iQ, '***');
                    end % p.value
                    % adjust p.value parameters
                    pval_hdl.Color = col.white;
                    pval_hdl.FontSize = 70;
                    pval_hdl.FontWeight = 'bold';
                    pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
                    pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
                end % when p.value is significant
            end % loop over questionnaires
        end % loop over aIns metabolites

        % for iQ = 1:n_quests
        %     quest_nm = quest_names{iQ};
        % %% plasma
        % fig;
        %
        % scat_hdl = scatter(plasma_Lac, questionnaires.(categ_nm).(quest_nm));
        % scat_hdl_upgrade(scat_hdl);
        %
        % % add fit
        % fit_hdl = plot(plasma_Lac_sorted.(categ_nm).(quest_nm),...
        %     questionnairess_fit_sorted_f_plasma.(categ_nm).(quest_nm));
        % fit_hdl_upgrade(fit_hdl);
        % % add pvalue and coeff correl
        % place_r_and_pval(rho.(categ_nm).(quest_nm).plasma,...
        %     pval.uncorrected.(categ_nm).(quest_nm).plasma);
        % % x/y labels
        % xlabel('Plasma lactate (mM)');
        % switch quest_nm
        %     case 'money'
        %         ylabel('Money (kCHF/year)');
        %     otherwise
        %         ylabel(quest_nm);
        % end
        % % figure title
        % switch categ_nm
        %     case 'dominance_compet'
            %         title('Dominance/Competitiveness');
            %     case 'stress_anxiety'
            %         title('Stress/Anxiety');
            %     otherwise
            %         title(categ_nm);
            % end
            % 
            % %% dmPFC
            % fig;
            % 
            % scat_hdl = scatter(dmPFC_Lac, questionnaires.(categ_nm).(quest_nm));
            % scat_hdl_upgrade(scat_hdl);
            % 
            % % add fit
            % fit_hdl = plot(dmPFC_Lac_sorted.(categ_nm).(quest_nm),...
            %     questionnairess_fit_sorted_f_dmPFC.(categ_nm).(quest_nm));
            % fit_hdl_upgrade(fit_hdl);
            % % add pvalue and coeff correl
            % place_r_and_pval(rho.(categ_nm).(quest_nm).dmPFC,...
            %     pval.uncorrected.(categ_nm).(quest_nm).dmPFC);
            % % x/y labels
            % xlabel('dmPFC/dACC lactate (mM)');
            % switch quest_nm
            %     case 'money'
            %         ylabel('Money (kCHF/year)');
            %     otherwise
            %         ylabel(quest_nm);
            % end
            % % figure title
            % switch categ_nm
            %     case 'dominance_compet'
            %         title('Dominance/Competitiveness');
            %     case 'stress_anxiety'
            %         title('Stress/Anxiety');
            %     otherwise
            %         title(categ_nm);
            % end
            % 
            % %% aIns
            % fig;
            % 
            % scat_hdl = scatter(aIns_Lac, questionnaires.(categ_nm).(quest_nm));
            % scat_hdl_upgrade(scat_hdl);
            % 
            % % add fit
            % fit_hdl = plot(aIns_Lac_sorted.(categ_nm).(quest_nm),...
            %     questionnairess_fit_sorted_f_aIns.(categ_nm).(quest_nm));
            % fit_hdl_upgrade(fit_hdl);
            % % add pvalue and coeff correl
            % place_r_and_pval(rho.(categ_nm).(quest_nm).aIns,...
            %     pval.uncorrected.(categ_nm).(quest_nm).aIns);
            % % x/y labels
            % xlabel('aIns lactate (mM)');
            % switch quest_nm
            %     case 'money'
            %         ylabel('Money (kCHF/year)');
            %     otherwise
            %         ylabel(quest_nm);
            % end
            % % figure title
            % switch categ_nm
            %     case 'dominance_compet'
            %         title('Dominance/Competitiveness');
            %     case 'stress_anxiety'
            %         title('Stress/Anxiety');
            %     otherwise
            %         title(categ_nm);
            % end
        % end % questionnaire loop
    end % figure
end % loop over questionnaire categories
% end % function