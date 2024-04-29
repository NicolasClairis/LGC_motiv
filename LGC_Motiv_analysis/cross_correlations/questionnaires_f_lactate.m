% function[]= questionnaires_f_lactate(fig_disp, rmv_outliers_yn)
% questionnaires_f_lactate will test whether the scores from the questionnaires
% correlate with lactate measures in the plasma or brain,
% after grouping them by category (stress/motivation/etc.).
%
% INPUTS
% fig_disp: display figures? (1) yes (0) no
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%

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
% plasma
plasmaM = load_plasma_metabolites(subject_id);
plasma_Lac = plasmaM.Lac./1000;
% brain
brainMetabolites = metabolite_load(subject_id);
dmPFC_Lac = brainMetabolites.dmPFC.Lac;
aIns_Lac = brainMetabolites.aIns.Lac;

%% extract questionnaires
[questionnaires, categ_quests, n_categ] = extract_questionnaires(study_nm, subject_id, NS);

%% remove outliers
if rmv_outliers_yn == 1
    % plasma
    [~,~,plasma_Lac] = rmv_outliers_3sd(plasma_Lac);
    % brain
    [~,~,dmPFC_Lac] = rmv_outliers_3sd(dmPFC_Lac);
    [~,~,aIns_Lac] = rmv_outliers_3sd(aIns_Lac);
    
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


%% method for multiple comparisons correction
corr_method = 'bonferroni';
corr_method_nm = ['corr_',corr_method];

%% perform correlations for each questionnaire
for iCateg = 1:n_categ
    categ_nm = categ_quests{iCateg};
    
    % category subfields
    quest_names = fieldnames(questionnaires.(categ_nm));
    n_quests = length(quest_names);
    
    [pval_unc_plasma, pval_unc_dmPFC, pval_unc_aIns] = deal([]);
    for iQ = 1:n_quests
        quest_nm = quest_names{iQ};
        
        % plasma
        goodS = ~isnan(plasma_Lac.*questionnaires.(categ_nm).(quest_nm));
        [rho.(categ_nm).(quest_nm).plasma, pval.uncorrected.(categ_nm).(quest_nm).plasma] = corr(...
            plasma_Lac(goodS)',...
            questionnaires.(categ_nm).(quest_nm)(goodS)');
        [~,~,~,~,plasma_Lac_sorted.(categ_nm).(quest_nm),...
            questionnairess_fit_sorted_f_plasma.(categ_nm).(quest_nm)] = glm_package(plasma_Lac', questionnaires.(categ_nm).(quest_nm)','normal');
        pval_unc_plasma = [pval_unc_plasma, pval.uncorrected.(categ_nm).(quest_nm).plasma];
        % dmPFC
        goodS = ~isnan(dmPFC_Lac.*questionnaires.(categ_nm).(quest_nm));
        [rho.(categ_nm).(quest_nm).dmPFC, pval.uncorrected.(categ_nm).(quest_nm).dmPFC] = corr(...
            dmPFC_Lac(goodS)',...
            questionnaires.(categ_nm).(quest_nm)(goodS)');
        [~,~,~,~,dmPFC_Lac_sorted.(categ_nm).(quest_nm),...
            questionnairess_fit_sorted_f_dmPFC.(categ_nm).(quest_nm)] = glm_package(dmPFC_Lac', questionnaires.(categ_nm).(quest_nm)','normal');
        pval_unc_dmPFC = [pval_unc_dmPFC, pval.uncorrected.(categ_nm).(quest_nm).dmPFC];
        % aIns
        goodS = ~isnan(aIns_Lac.*questionnaires.(categ_nm).(quest_nm));
        [rho.(categ_nm).(quest_nm).aIns, pval.uncorrected.(categ_nm).(quest_nm).aIns] = corr(...
            aIns_Lac(goodS)',...
            questionnaires.(categ_nm).(quest_nm)(goodS)');
        [~,~,~,~,aIns_Lac_sorted.(categ_nm).(quest_nm),...
            questionnairess_fit_sorted_f_aIns.(categ_nm).(quest_nm)] = glm_package(aIns_Lac', questionnaires.(categ_nm).(quest_nm)','normal');
        pval_unc_aIns = [pval_unc_aIns, pval.uncorrected.(categ_nm).(quest_nm).aIns];
        
        %% store significant results in a different subfield
        % plasma
        if pval.uncorrected.(categ_nm).(quest_nm).plasma < 0.05
            % store coefficient
            signif.unc.(categ_nm).(quest_nm).plasma.r_corr = rho.(categ_nm).(quest_nm).plasma;
            % store p.value
            signif.unc.(categ_nm).(quest_nm).plasma.pval = pval.uncorrected.(categ_nm).(quest_nm).plasma;
        end
        % dmPFC
        if pval.uncorrected.(categ_nm).(quest_nm).dmPFC < 0.05
            % store coefficient
            signif.unc.(categ_nm).(quest_nm).dmPFC.r_corr = rho.(categ_nm).(quest_nm).dmPFC;
            % store p.value
            signif.unc.(categ_nm).(quest_nm).dmPFC.pval = pval.uncorrected.(categ_nm).(quest_nm).dmPFC;
        end
        % aIns
        if pval.uncorrected.(categ_nm).(quest_nm).aIns < 0.05
            % store coefficient
            signif.unc.(categ_nm).(quest_nm).aIns.r_corr = rho.(categ_nm).(quest_nm).aIns;
            % store p.value
            signif.unc.(categ_nm).(quest_nm).aIns.pval = pval.uncorrected.(categ_nm).(quest_nm).aIns;
        end
    end % loop over questionnaires
    
    %% attempt to correct for multiple comparisons
    % use pval_adjust if you want to correct for multiple comparisons
    % temporarily decided to correct for the number of dimensions
    % tested (considering that the tests should not be considered to be independent within the same dimension)
    % but note this is quite arbitrary... Another way would be to
    % correct for the number of questionnaires or for the number of
    % scores extracted, but again the independence between tests is
    % highly debatable...
    [pval_corr_plasma] = pval_adjust(pval_unc_plasma,corr_method);
    [pval_corr_dmPFC] = pval_adjust(pval_unc_dmPFC,corr_method);
    [pval_corr_aIns] = pval_adjust(pval_unc_aIns,corr_method);
    % extract individual corrected p.values
    for iQ = 1:n_quests
        quest_nm = quest_names{iQ};
        pval.(corr_method_nm).(categ_nm).(quest_nm).plasma = pval_corr_plasma(iQ);
        pval.(corr_method_nm).(categ_nm).(quest_nm).dmPFC = pval_corr_dmPFC(iQ);
        pval.(corr_method_nm).(categ_nm).(quest_nm).aIns = pval_corr_aIns(iQ);
        
        % store significant results in a different subfield
        % plasma
        if pval.(corr_method_nm).(categ_nm).(quest_nm).plasma < 0.05
            % store coefficient
            signif.(corr_method_nm).(categ_nm).(quest_nm).plasma.r_corr = rho.(categ_nm).(quest_nm).plasma;
            % store p.value
            signif.(corr_method_nm).(categ_nm).(quest_nm).plasma.pval = pval.(corr_method_nm).(categ_nm).(quest_nm).plasma;
        end
        % dmPFC
        if pval.(corr_method_nm).(categ_nm).(quest_nm).dmPFC < 0.05
            % store coefficient
            signif.(corr_method_nm).(categ_nm).(quest_nm).dmPFC.r_corr = rho.(categ_nm).(quest_nm).dmPFC;
            % store p.value
            signif.(corr_method_nm).(categ_nm).(quest_nm).dmPFC.pval = pval.(corr_method_nm).(categ_nm).(quest_nm).dmPFC;
        end
        % aIns
        if pval.(corr_method_nm).(categ_nm).(quest_nm).aIns < 0.05
            % store coefficient
            signif.(corr_method_nm).(categ_nm).(quest_nm).aIns.r_corr = rho.(categ_nm).(quest_nm).aIns;
            % store p.value
            signif.(corr_method_nm).(categ_nm).(quest_nm).aIns.pval = pval.(corr_method_nm).(categ_nm).(quest_nm).aIns;
        end
    end % loop over questionnaires
    
    %% show results for current in a graph (+ if results are significant or not)
    if fig_disp == 1
        % general figure parameters
        [pSize, lW, col, mSize] = general_fig_prm;
        
        for iQ = 1:n_quests
            quest_nm = quest_names{iQ};
            
            %% plasma
            fig;
            
            scat_hdl = scatter(plasma_Lac, questionnaires.(categ_nm).(quest_nm));
            scat_hdl_upgrade(scat_hdl);
            
            % add fit
            fit_hdl = plot(plasma_Lac_sorted.(categ_nm).(quest_nm),...
                questionnairess_fit_sorted_f_plasma.(categ_nm).(quest_nm));
            fit_hdl_upgrade(fit_hdl);
            % add pvalue and coeff correl
            place_r_and_pval(rho.(categ_nm).(quest_nm).plasma,...
                pval.uncorrected.(categ_nm).(quest_nm).plasma);
            % x/y labels
            xlabel('Plasma lactate (mM)');
            switch quest_nm
                case 'money'
                    ylabel('Money (kCHF/year)');
                otherwise
                    ylabel(quest_nm);
            end
            % figure title
            switch categ_nm
                case 'dominance_compet'
                    title('Dominance/Competitiveness');
                case 'stress_anxiety'
                    title('Stress/Anxiety');
                otherwise
                    title(categ_nm);
            end
            
            %% dmPFC
            fig;
            
            scat_hdl = scatter(dmPFC_Lac, questionnaires.(categ_nm).(quest_nm));
            scat_hdl_upgrade(scat_hdl);
            
            % add fit
            fit_hdl = plot(dmPFC_Lac_sorted.(categ_nm).(quest_nm),...
                questionnairess_fit_sorted_f_dmPFC.(categ_nm).(quest_nm));
            fit_hdl_upgrade(fit_hdl);
            % add pvalue and coeff correl
            place_r_and_pval(rho.(categ_nm).(quest_nm).dmPFC,...
                pval.uncorrected.(categ_nm).(quest_nm).dmPFC);
            % x/y labels
            xlabel('dmPFC/dACC lactate (mM)');
            switch quest_nm
                case 'money'
                    ylabel('Money (kCHF/year)');
                otherwise
                    ylabel(quest_nm);
            end
            % figure title
            switch categ_nm
                case 'dominance_compet'
                    title('Dominance/Competitiveness');
                case 'stress_anxiety'
                    title('Stress/Anxiety');
                otherwise
                    title(categ_nm);
            end
            
            %% aIns
            fig;
            
            scat_hdl = scatter(aIns_Lac, questionnaires.(categ_nm).(quest_nm));
            scat_hdl_upgrade(scat_hdl);
            
            % add fit
            fit_hdl = plot(aIns_Lac_sorted.(categ_nm).(quest_nm),...
                questionnairess_fit_sorted_f_aIns.(categ_nm).(quest_nm));
            fit_hdl_upgrade(fit_hdl);
            % add pvalue and coeff correl
            place_r_and_pval(rho.(categ_nm).(quest_nm).aIns,...
                pval.uncorrected.(categ_nm).(quest_nm).aIns);
            % x/y labels
            xlabel('aIns lactate (mM)');
            switch quest_nm
                case 'money'
                    ylabel('Money (kCHF/year)');
                otherwise
                    ylabel(quest_nm);
            end
            % figure title
            switch categ_nm
                case 'dominance_compet'
                    title('Dominance/Competitiveness');
                case 'stress_anxiety'
                    title('Stress/Anxiety');
                otherwise
                    title(categ_nm);
            end
        end % questionnaire loop
    end % figure
end % loop over questionnaire categories
% end % function