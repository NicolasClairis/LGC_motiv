function[betas, pval] = brain_GSH_f_nutrition_precursors(figDisp)
% [betas, pval] = brain_GSH_f_nutrition_precursors(figDisp)
%% brain_GSH_f_nutrition_precursors will test whether there is any 
% correlation between nutrition precursors of GSH (namely Gly, Cys and Glu)
% and brain levels of GSH in the dmPFC or aINS (for study 1)
%
% INPUTS
% figDisp: display figure (1) or not (0) ? Will be equal to 1 by default
% if not enterrd
%
% OUTPUTS
% betas: structure with betas for each test
%
% pval: structure with p.value for each test

%% main parameters
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
end
%% working directories
root = LGCM_root_paths;
% study
study_nm = 'study1';
studyPath = [root, filesep, study_nm, filesep];
% nutrition data
switch root
    case ['E:',filesep]
        gitPath = fullfile('C:','Users','clairis','Desktop');
    case {[fullfile('C:','Users','Loco','Downloads'),filesep],...
            [fullfile('L:','human_data_private','raw_data_subject'),filesep]}
        gitPath = fullfile('C:','Users','Loco','Downloads');
    otherwise
        error('case not ready yet');
end
nutritionPath = [fullfile(gitPath,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'nutrition'),filesep];

% condition
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% initialize variables of interest
[dmPFC_GSH, aINS_GSH,...
    nutri.Gly, nutri.Cys, nutri.Glu,...
    nutri.calories,...
    nutri.GlyCysGlu_sum,...
    nutri.Gly_div_totalCal,  nutri.Cys_div_totalCal, nutri.Glu_div_totalCal,...
    nutri.GlyCysGlu_sum_div_totalCal] = deal(NaN(1,NS));
nutri_vars = fieldnames(nutri);
% load 
Gly_table = readtable([nutritionPath, 'Glycine_scoring.xlsx'],...
    'Sheet','Sheet1');
Cys_table = readtable([nutritionPath, 'cystein_scoring.xlsx'],...
    'Sheet','Sheet1');
Glu_table = readtable([nutritionPath, 'Glutamate_scoring.xlsx'],...
    'Sheet','Sheet1');
calories_table = readtable([nutritionPath, 'calories_scoring.xlsx'],...
    'Sheet','Sheet1');
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    %% load brain metabolites
    [metabolite_tmp] = metabolite_load({sub_nm});
    dmPFC_GSH(iS) = metabolite_tmp.dmPFC.GSH;
    aINS_GSH(iS) = metabolite_tmp.aIns.GSH;
    
    %% load nutrition score
    sub_Gly_idx = find(strcmp(Gly_table.CID, sub_nm));
    sub_Cys_idx = find(strcmp(Cys_table.CID, sub_nm));
    sub_Glu_idx = find(strcmp(Glu_table.CID, sub_nm));
    sub_calories_idx = find(strcmp(calories_table.CID, sub_nm));
    % extract nutrition intake values
    if ~isempty(sub_Gly_idx) && size(sub_Gly_idx,1) == 1
        nutri.Gly(iS) = Gly_table.GlyParSemaine_mg_(sub_Gly_idx);
    end
    if ~isempty(sub_Cys_idx) && size(sub_Cys_idx,1) == 1
        nutri.Cys(iS) = Cys_table.CysteineParSemaine_mg_(sub_Cys_idx);
    end
    if ~isempty(sub_Glu_idx) && size(sub_Glu_idx,1) == 1
        nutri.Glu(iS) = Glu_table.GlutamateParSemaine_mg_(sub_Glu_idx);
    end
    if ~isempty(sub_calories_idx) && size(sub_calories_idx,1) == 1
        nutri.calories(iS) = calories_table.CaloriesParSemaine_kcal_(sub_calories_idx);
    end
    nutri.GlyCysGlu_sum(iS) = nutri.Gly(iS) + nutri.Cys(iS) + nutri.Glu(iS);
    
    % divide by calories
    if ~isnan(nutri.calories(iS))
        nutri.Gly_div_totalCal(iS) = nutri.Gly(iS)/nutri.calories(iS);
        nutri.Cys_div_totalCal(iS) = nutri.Cys(iS)/nutri.calories(iS);
        nutri.Glu_div_totalCal(iS) = nutri.Gly(iS)/nutri.calories(iS);
        nutri.GlyCysGlu_sum_div_totalCal(iS) = nutri.GlyCysGlu_sum(iS)/nutri.calories(iS);
    end
end % subject loop

%% test correlations
for iROI = 1:2
    switch iROI
        case 1
            ROI_nm = 'dmPFC_GSH';
            ROI_var = dmPFC_GSH;
        case 2
            ROI_nm = 'aINS_GSH';
            ROI_var = aINS_GSH;
    end
    for iNutri = 1:length(nutri_vars)
        nutri_nm = nutri_vars{iNutri};
        % filter bad subjects
        subs_ok = (~isnan(nutri.(nutri_nm)).*~isnan(ROI_var)) == 1;
        if sum(subs_ok) > 0
            % perform glm brain GSH = b0 + b1*nutrition variable
            [betas.(ROI_nm).(['f_',nutri_nm]),~,stats_tmp] = glmfit(nutri.(nutri_nm)(subs_ok), ROI_var(subs_ok),'normal');
            % extract p.value
            pval.(ROI_nm).(['f_',nutri_nm]) = stats_tmp.p;
            % extract the fit
            nutri_bis.(ROI_nm).(nutri_nm) = sort(nutri.(nutri_nm)(subs_ok));
            fittedData.(ROI_nm).(['f_',nutri_nm]) = glmval(betas.(ROI_nm).(['f_',nutri_nm]), nutri_bis.(ROI_nm).(nutri_nm), 'identity');
        end
    end % loop through nutritional intake variables
end % loop through ROIs

%% figure
if figDisp == 1
    
    dmPFC_col = 'm';
    aINS_col = 'b';
    lWidth = 3;
    pSize = 25;
    
    %% show result Glu/Cys/Gly alone
    fig;
    for iNutri = 1:3
        switch iNutri
            case 1
                nutri_nm = 'Glu';
                nutri_nm_bis = 'Glu_div_totalCal';
            case 2
                nutri_nm = 'Cys';
                nutri_nm_bis = 'Cys_div_totalCal';
            case 3
                nutri_nm = 'Gly';
                nutri_nm_bis = 'Gly_div_totalCal';
        end
        
        % raw metabolite
        subplot(2,3,iNutri);
        hold on;
        % show raw data
        dmPFC_hdl = scatter(nutri.(nutri_nm), dmPFC_GSH);
        aINS_hdl = scatter(nutri.(nutri_nm), aINS_GSH);
        dmPFC_hdl.LineWidth = lWidth;
        aINS_hdl.LineWidth = lWidth;
        dmPFC_hdl.MarkerEdgeColor = dmPFC_col;
        aINS_hdl.MarkerEdgeColor = aINS_col;
        % add the fit
        dmPFC_fit_hdl = plot(nutri_bis.dmPFC_GSH.(nutri_nm),...
            fittedData.dmPFC_GSH.(['f_',(nutri_nm)]));
        aINS_fit_hdl = plot(nutri_bis.aINS_GSH.(nutri_nm),...
            fittedData.aINS_GSH.(['f_',(nutri_nm)]));
        dmPFC_fit_hdl.LineStyle = '--';
        dmPFC_fit_hdl.LineWidth = lWidth;
        dmPFC_fit_hdl.Color = dmPFC_col;
        aINS_fit_hdl.LineStyle = '--';
        aINS_fit_hdl.LineWidth = lWidth;
        aINS_fit_hdl.Color = aINS_col;
        % add legend
        legend([dmPFC_hdl, aINS_hdl],{'dmPFC','aINS'});
        legend('boxoff');
        xlabel([nutri_nm,' (mg/week)']);
        ylabel('GSH (μmol/g)');
        legend_size(pSize);
        
        % Nutrition/total calories
        subplot(2,3,iNutri+3);
        hold on;
        % show data normalized
        dmPFC_hdl = scatter(nutri.(nutri_nm_bis), dmPFC_GSH);
        aINS_hdl = scatter(nutri.(nutri_nm_bis), aINS_GSH);
        dmPFC_hdl.LineWidth = lWidth;
        aINS_hdl.LineWidth = lWidth;
        dmPFC_hdl.MarkerEdgeColor = dmPFC_col;
        aINS_hdl.MarkerEdgeColor = aINS_col;
        % add the fit
        dmPFC_fit_hdl = plot(nutri_bis.dmPFC_GSH.(nutri_nm_bis),...
            fittedData.dmPFC_GSH.(['f_',(nutri_nm_bis)]));
        aINS_fit_hdl = plot(nutri_bis.aINS_GSH.(nutri_nm_bis),...
            fittedData.aINS_GSH.(['f_',(nutri_nm_bis)]));
        dmPFC_fit_hdl.LineStyle = '--';
        dmPFC_fit_hdl.LineWidth = lWidth;
        dmPFC_fit_hdl.Color = dmPFC_col;
        aINS_fit_hdl.LineStyle = '--';
        aINS_fit_hdl.LineWidth = lWidth;
        aINS_fit_hdl.Color = aINS_col;
        % add legend
        legend([dmPFC_hdl, aINS_hdl],{'dmPFC','aINS'});
        legend('boxoff');
        xlabel([nutri_nm,'/calories']);
        ylabel('GSH (μmol/g)');
        legend_size(pSize);
    end % nutrition metabolite
    
    %% show sum
    fig;
    % Gly + Cys + Glu
    subplot(2,1,1);
    hold on;
    dmPFC_hdl = scatter(nutri.GlyCysGlu_sum, dmPFC_GSH);
    aINS_hdl = scatter(nutri.GlyCysGlu_sum, aINS_GSH);
    dmPFC_hdl.LineWidth = lWidth;
    aINS_hdl.LineWidth = lWidth;
    dmPFC_hdl.MarkerEdgeColor = dmPFC_col;
    aINS_hdl.MarkerEdgeColor = aINS_col;
    % add the fit
    dmPFC_fit_hdl = plot(nutri_bis.dmPFC_GSH.GlyCysGlu_sum,...
        fittedData.dmPFC_GSH.f_GlyCysGlu_sum);
    aINS_fit_hdl = plot(nutri_bis.aINS_GSH.GlyCysGlu_sum,...
        fittedData.aINS_GSH.f_GlyCysGlu_sum);
    dmPFC_fit_hdl.LineStyle = '--';
    dmPFC_fit_hdl.LineWidth = lWidth;
    dmPFC_fit_hdl.Color = dmPFC_col;
    aINS_fit_hdl.LineStyle = '--';
    aINS_fit_hdl.LineWidth = lWidth;
    aINS_fit_hdl.Color = aINS_col;
    % add legend
    legend([dmPFC_hdl, aINS_hdl],{'dmPFC','aINS'});
    legend('boxoff');
    xlabel('Glu + Cys + Gly (mg/week)');
    ylabel('GSH (μmol/g)');
    legend_size(pSize);
    
    % (Gly + Cys + Glu)/(total calories)
    subplot(2,1,2);
    hold on;
    dmPFC_hdl = scatter(nutri.GlyCysGlu_sum_div_totalCal, dmPFC_GSH);
    aINS_hdl = scatter(nutri.GlyCysGlu_sum_div_totalCal, aINS_GSH);
    dmPFC_hdl.LineWidth = lWidth;
    aINS_hdl.LineWidth = lWidth;
    dmPFC_hdl.MarkerEdgeColor = dmPFC_col;
    aINS_hdl.MarkerEdgeColor = aINS_col;
    % add the fit
    dmPFC_fit_hdl = plot(nutri_bis.dmPFC_GSH.GlyCysGlu_sum_div_totalCal,...
        fittedData.dmPFC_GSH.f_GlyCysGlu_sum_div_totalCal);
    aINS_fit_hdl = plot(nutri_bis.aINS_GSH.GlyCysGlu_sum_div_totalCal,...
        fittedData.aINS_GSH.f_GlyCysGlu_sum_div_totalCal);
    dmPFC_fit_hdl.LineStyle = '--';
    dmPFC_fit_hdl.LineWidth = lWidth;
    dmPFC_fit_hdl.Color = dmPFC_col;
    aINS_fit_hdl.LineStyle = '--';
    aINS_fit_hdl.LineWidth = lWidth;
    aINS_fit_hdl.Color = aINS_col;
    % add legend
    legend([dmPFC_hdl, aINS_hdl],{'dmPFC','aINS'});
    legend('boxoff');
    xlabel('(Glu + Cys + Gly)/calories');
    ylabel('GSH (μmol/g)');
    legend_size(pSize);
end % figure display
end % function