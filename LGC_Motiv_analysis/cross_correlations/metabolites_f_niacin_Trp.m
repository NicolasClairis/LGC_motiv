%% test correlation between brain GSH and niacin


if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 0;
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
[nutri.niacin, nutri.Trp, nutri.niacinEquiv, nutri.niacinPlusTrp,...
    nutri.calories,...
    nutri.niacin_div_totalCal,  nutri.Trp_div_totalCal,...
    nutri.niacinEquiv,...
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
%% extract niacin equivalents data
niacineEquivFilePath = [nutritionPath,'niacine_equivalents_scoring.xlsx'];
niacinEquiv_table = readtable(niacineEquivFilePath,...
    'Sheet','Sheet1');

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
switch study_nm
    case 'study1'
        MRS_ROIs = {'dmPFC','aIns'};
    case 'study2'
        error('not ready yet');
end
nROIs = length(MRS_ROIs);
for iMRS_ROI = 1:nROIs
    MRS_ROI_nm = MRS_ROIs{iMRS_ROI};
    metabolite_names.(MRS_ROI_nm) = fieldnames(metabolites.(MRS_ROI_nm));
    n_metabolites.(MRS_ROI_nm) = length(metabolite_names.(MRS_ROI_nm));
end % roi loop

%% loop through subjects to extract nutrition
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    %% load nutrition score
    sub_niacin_idx = find(strcmp(niacin_table.CID, sub_nm));
    sub_Trp_idx = find(strcmp(Trp_table.CID, sub_nm));
    sub_niacinEquiv_idx = find(strcmp(niacinEquiv_table.CID, sub_nm));
    sub_calories_idx = find(strcmp(calories_table.CID, sub_nm));
    % extract nutrition intake values
    % extract niacin values
    if ~isempty(sub_niacin_idx) && size(sub_niacin_idx,1) == 1
        nutri.niacin(iS) = niacin_table.NiacineParSemaine__g_(sub_niacin_idx);
    end
    % extract Tryptophan values
    if ~isempty(sub_Trp_idx) && size(sub_Trp_idx,1) == 1
        nutri.Trp(iS) = Trp_table.TryptophaneParSemaine_mg_(sub_Trp_idx);
    end
    % extract niacin equivalents values
    if ~isempty(sub_niacinEquiv_idx) && size(sub_niacinEquiv_idx,1) == 1
        nutri.niacinEquiv(iS) = niacinEquiv_table.x_quivalentsDeNiacineParSemaine__g_(sub_niacinEquiv_idx);
    end
    % extract calories values
    if ~isempty(sub_calories_idx) && size(sub_calories_idx,1) == 1
        nutri.calories(iS) = calories_table.CaloriesParSemaine_kcal_(sub_calories_idx);
    end
    nutri.niacinPlusTrp(iS) = nutri.niacin(iS)/1000 + nutri.Trp(iS);
    
    % divide by calories
    if ~isnan(nutri.calories(iS))
        nutri.niacin_div_totalCal(iS) = nutri.niacin(iS)/nutri.calories(iS);
        nutri.niacinEquiv_div_totalCal(iS) = nutri.niacinEquiv(iS)/nutri.calories(iS);
        nutri.Trp_div_totalCal(iS) = nutri.Trp(iS)/nutri.calories(iS);
        nutri.niacinPlusTrp_div_totalCal(iS) = nutri.niacinPlusTrp(iS)/nutri.calories(iS);
    end
end % subject loop

%% test correlations
for iMRS_ROI = 1:nROIs
    MRS_ROI_nm = MRS_ROIs{iMRS_ROI};
    
    for iMb = 1:n_metabolites.(MRS_ROI_nm)
        mb_nm = metabolite_names.(MRS_ROI_nm){iMb};
        ROI_mb_nm = [MRS_ROI_nm,'_',mb_nm];
        ROI_mb_var = metabolites.(MRS_ROI_nm).(mb_nm);
        for iNutri = 1:length(nutri_vars)
            nutri_nm = nutri_vars{iNutri};
            % filter bad subjects
            subs_ok = (~isnan(nutri.(nutri_nm)).*~isnan(ROI_mb_var)) == 1;
            if sum(subs_ok) > 0
                % perform glm brain GSH = b0 + b1*nutrition variable
                [betas.(ROI_mb_nm).(['f_',nutri_nm]),~,stats_tmp] = glmfit(nutri.(nutri_nm)(subs_ok), ROI_mb_var(subs_ok),'normal');
                [r_corr.(ROI_mb_nm).(['f_',nutri_nm])] = glmfit(zscore(nutri.(nutri_nm)(subs_ok)), zscore(ROI_mb_var(subs_ok)),'normal');
                % extract p.value
                pval.(ROI_mb_nm).(['f_',nutri_nm]) = stats_tmp.p;
                if stats_tmp.p(2) < 0.05
                    pval_signif.(ROI_mb_nm).(['f_',nutri_nm]).pval = stats_tmp.p(2);
                    pval_signif.(ROI_mb_nm).(['f_',nutri_nm]).r_corr = r_corr.(ROI_mb_nm).(['f_',nutri_nm])(2);
                end
                % extract the fit
                nutri_bis.(ROI_mb_nm).(nutri_nm) = sort(nutri.(nutri_nm)(subs_ok));
                fittedData.(ROI_mb_nm).(['f_',nutri_nm]) = glmval(betas.(ROI_mb_nm).(['f_',nutri_nm]), nutri_bis.(ROI_mb_nm).(nutri_nm), 'identity');
            end
        end % loop through nutritional intake variables
    end % loop through metabolites
end % loop through ROIs

%% figure
if figDisp == 1
    
    dmPFC_col = 'm';
    aINS_col = 'b';
    lWidth = 3;
    pSize = 25;
    
    for iMb = 1:n_metabolites.(MRS_ROIs{1})
        mb_nm = metabolite_names.(MRS_ROIs{1}){iMb};
        
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
                    nutri_nm = 'niacinEquiv';
                    nutri_nm_bis = 'niacinEquiv_div_totalCal';
%                     nutri_nm = 'niacinPlusTrp';
%                     nutri_nm_bis = 'niacinPlusTrp_div_totalCal';
            end
            
            % raw metabolite
            subplot(2,3,iNutri);
            hold on;
            % show raw data
            dmPFC_hdl = scatter(nutri.(nutri_nm), metabolites.dmPFC.(mb_nm));
            aINS_hdl = scatter(nutri.(nutri_nm), metabolites.aIns.(mb_nm));
            dmPFC_hdl.LineWidth = lWidth;
            aINS_hdl.LineWidth = lWidth;
            dmPFC_hdl.MarkerEdgeColor = dmPFC_col;
            aINS_hdl.MarkerEdgeColor = aINS_col;
            % add the fit
            dmPFC_fit_hdl = plot(nutri_bis.(['dmPFC_',mb_nm]).(nutri_nm),...
                fittedData.(['dmPFC_',mb_nm]).(['f_',(nutri_nm)]));
            aINS_fit_hdl = plot(nutri_bis.(['aIns_',mb_nm]).(nutri_nm),...
                fittedData.(['aIns_',mb_nm]).(['f_',(nutri_nm)]));
            dmPFC_fit_hdl.LineStyle = '--';
            dmPFC_fit_hdl.LineWidth = lWidth;
            dmPFC_fit_hdl.Color = dmPFC_col;
            aINS_fit_hdl.LineStyle = '--';
            aINS_fit_hdl.LineWidth = lWidth;
            aINS_fit_hdl.Color = aINS_col;
            % add legend
            legend([dmPFC_hdl, aINS_hdl],{'dmPFC','aINS'});
            legend('boxoff');
            switch nutri_nm
                case 'niacin'
                    xlabel([nutri_nm,' (μg/week)']);
                case 'Trp'
                    xlabel([nutri_nm,' (mg/week)']);
                case 'niacinEquiv'
                    xlabel('niacin equivalents (mg/week)');
                case 'niacinPlusTrp'
                    xlabel('niacin + Trp (mg/week)');
            end
            ylabel([mb_nm,' (μmol/g)']);
            legend_size(pSize);
            
            % Nutrition/total calories
            subplot(2,3,iNutri+3);
            hold on;
            % show data normalized
            dmPFC_hdl = scatter(nutri.(nutri_nm_bis), metabolites.dmPFC.(mb_nm));
            aINS_hdl = scatter(nutri.(nutri_nm_bis), metabolites.aIns.(mb_nm));
            dmPFC_hdl.LineWidth = lWidth;
            aINS_hdl.LineWidth = lWidth;
            dmPFC_hdl.MarkerEdgeColor = dmPFC_col;
            aINS_hdl.MarkerEdgeColor = aINS_col;
            % add the fit
            dmPFC_fit_hdl = plot(nutri_bis.(['dmPFC_',mb_nm]).(nutri_nm_bis),...
                fittedData.(['dmPFC_',mb_nm]).(['f_',(nutri_nm_bis)]));
            aINS_fit_hdl = plot(nutri_bis.(['aIns_',mb_nm]).(nutri_nm_bis),...
                fittedData.(['aIns_',mb_nm]).(['f_',(nutri_nm_bis)]));
            dmPFC_fit_hdl.LineStyle = '--';
            dmPFC_fit_hdl.LineWidth = lWidth;
            dmPFC_fit_hdl.Color = dmPFC_col;
            aINS_fit_hdl.LineStyle = '--';
            aINS_fit_hdl.LineWidth = lWidth;
            aINS_fit_hdl.Color = aINS_col;
            % add legend
            legend([dmPFC_hdl, aINS_hdl],{'dmPFC','aINS'});
            legend('boxoff');
            
            switch nutri_nm
                case {'niacin','Trp'}
                    xlabel([nutri_nm,'/calories']);
                case 'niacinPlusTrp'
                    xlabel('(niacin + Trp)/calories');
                case 'niacinEquiv'
                    xlabel('niacin equivalents/calories (mg/week)');
            end
            ylabel([mb_nm,' (μmol/g)']);
            legend_size(pSize);
        end % nutrition metabolite
        
    end % loop through metabolites
end % figure display