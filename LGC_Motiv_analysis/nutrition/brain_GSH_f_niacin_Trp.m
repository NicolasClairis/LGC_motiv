%% test correlation between brain GSH and niacin


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
    'LGC_Motiv_results','nutrition'),filesep];
% condition
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);
%% initialize variables of interest
[dmPFC_GSH, aINS_GSH,...
    nutri.niacin, nutri.Trp, nutri.niacinPlusTrp,...
    nutri.calories,...
    nutri.niacin_div_totalCal,  nutri.Trp_div_totalCal,...
    nutri.niacinPlusTrp_div_totalCal] = deal(NaN(1,NS));
nutri_vars = fieldnames(nutri);
%% extract all calories data
calories_table = readtable([nutritionPath, 'calories_scoring.xlsx'],...
    'Sheet','Sheet1');
%% extract all niacin data
niacineFilePath = [nutritionPath,'niacine_scoring.xlsx'];
niacin_table = readtable(niacineFilePath,...
    'Sheet','Sheet1');
%% extract all Tryptophane data
TrpFilePath = [nutritionPath,'Tryptophan_scoring.xlsx'];
Trp_table = readtable(TrpFilePath,...
    'Sheet','Sheet1');

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    %% load brain metabolites
    [metabolite_tmp] = metabolite_load({sub_nm});
    dmPFC_GSH(iS) = metabolite_tmp.dmPFC.GSH;
    aINS_GSH(iS) = metabolite_tmp.aIns.GSH;
    
    %% load nutrition score
    sub_niacin_idx = find(strcmp(niacin_table.CID, sub_nm));
    sub_Trp_idx = find(strcmp(Trp_table.CID, sub_nm));
    sub_calories_idx = find(strcmp(calories_table.CID, sub_nm));
    % extract nutrition intake values
    % extract niacin values
    if ~isempty(sub_niacin_idx) && size(sub_niacin_idx,1) == 1
        nutri.niacin(iS) = niacin_table.NiacineParSemaine_ug_(sub_niacin_idx);
    end
    % extract Tryptophan values
    if ~isempty(sub_Trp_idx) && size(sub_Trp_idx,1) == 1
        nutri.Trp(iS) = Trp_table.TryptophaneParSemaine_ug_(sub_Trp_idx);
    end
    % extract calories values
    if ~isempty(sub_calories_idx) && size(sub_calories_idx,1) == 1
        nutri.calories(iS) = calories_table.CaloriesParSemaine_kcal_(sub_calories_idx);
    end
    nutri.niacinPlusTrp(iS) = nutri.niacin(iS) + nutri.Trp(iS);
    
    % divide by calories
    if ~isnan(nutri.calories(iS))
        nutri.niacin_div_totalCal(iS) = nutri.niacin(iS)/nutri.calories(iS);
        nutri.Trp_div_totalCal(iS) = nutri.Trp(iS)/nutri.calories(iS);
        nutri.niacinPlusTrp_div_totalCal(iS) = nutri.niacinPlusTrp(iS)/nutri.calories(iS);
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
    
    %% show results
    fig;
    for iNutri = 1:3
        switch iNutri
            case 1
                nutri_nm = 'niacin';
                nutri_nm_bis = 'niacin_div_totalCal';
            case 2
                nutri_nm = 'Trp';
                nutri_nm_bis = 'Trp_div_totalCal';
            case 3
                nutri_nm = 'niacinPlusTrp';
                nutri_nm_bis = 'niacinPlusTrp_div_totalCal';
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
        if ismember(nutri_nm,{'niacin','Trp'})
            xlabel([nutri_nm,' (μg/week)']);
        elseif strcmp(nutri_nm,'niacinPlusTrp')
            xlabel('niacin + Trp (μg/week)');
        end
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
        
        if ismember(nutri_nm,{'niacin','Trp'})
            xlabel([nutri_nm,'/calories']);
        elseif strcmp(nutri_nm,'niacinPlusTrp')
            xlabel('(niacin + Trp)/calories');
        end
        ylabel('GSH (μmol/g)');
        legend_size(pSize);
    end % nutrition metabolite
    
end % figure display