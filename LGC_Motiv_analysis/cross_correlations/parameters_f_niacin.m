%% test whether parameters from the bayesian model are correlated to the nutritional intake of niacin/Trp/niacin equivalents

%% by default display the figures but can be set to 0 if you only want the p.value
figDisp = 1;

%% define subjects
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

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
        gitPath = fullfile('C:','Users','Loco','Documents');
    otherwise
        error('case not ready yet');
end
nutritionPath = [fullfile(gitPath,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'nutrition'),filesep];

%% load nutrition
[nutri.niacin, nutri.Trp, nutri.niacinEquiv,...
    nutri.calories,...
    nutri.niacin_div_totalCal,  nutri.Trp_div_totalCal,...
    nutri.niacinEquiv_div_totalCal] = deal(NaN(1,NS));
nutri_vars = fieldnames(nutri);
% extract all calories data
calories_table = readtable([nutritionPath, 'calories_scoring.xlsx'],...
    'Sheet','Sheet1');
% extract all niacin data
niacineFilePath = [nutritionPath,'niacine_scoring.xlsx'];
niacin_table = readtable(niacineFilePath,...
    'Sheet','Sheet1');
% extract all Tryptophane data
TrpFilePath = [nutritionPath,'Tryptophan_scoring.xlsx'];
Trp_table = readtable(TrpFilePath,...
    'Sheet','Sheet1');
% extract all niacin equivalents data
niacineEquivFilePath = [nutritionPath,'niacine_equivalents_scoring.xlsx'];
niacinEquiv_table = readtable(niacineEquivFilePath,...
    'Sheet','Sheet1');

%% load parameters
[prm_all, mdlType, mdlN] = prm_extraction(study_nm, subject_id);
parameters = fieldnames(prm_all);
n_prm = length(parameters);
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    if ~strcmp(prm_nm,'CID')
        prm.(prm_nm) = NaN(1,NS);
    end
end

%% loop through subjects
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
        nutri.niacin(iS) = niacin_table.NiacineParSemaine__g_(sub_niacin_idx)./1000;
    end
    % extract Tryptophan values
    if ~isempty(sub_Trp_idx) && size(sub_Trp_idx,1) == 1
        nutri.Trp(iS) = Trp_table.TryptophaneParSemaine_mg_(sub_Trp_idx);
    end
    % extract niacin equivalents values
    if ~isempty(sub_niacinEquiv_idx) && size(sub_niacinEquiv_idx,1) == 1
        nutri.niacinEquiv(iS) = niacinEquiv_table.x_quivalentsDeNiacineParSemaine__g_(sub_niacinEquiv_idx)./1000;
    end
    % extract calories values
    if ~isempty(sub_calories_idx) && size(sub_calories_idx,1) == 1
        nutri.calories(iS) = calories_table.CaloriesParSemaine_kcal_(sub_calories_idx);
    end
    
    % divide by calories
    if ~isnan(nutri.calories(iS))
        nutri.niacin_div_totalCal(iS) = nutri.niacin(iS)/nutri.calories(iS);
        nutri.Trp_div_totalCal(iS) = nutri.Trp(iS)/nutri.calories(iS);
        nutri.niacinEquiv_div_totalCal(iS) = nutri.niacinEquiv(iS)/nutri.calories(iS);
    end
    
    %% load parameters
    jS = strcmp(sub_nm,prm_all.CID);
    for iPrm = 1:n_prm
        prm_nm = parameters{iPrm};
        if ~strcmp(prm_nm,'CID')
            prm.(prm_nm)(iS) = prm_all.(prm_nm)(jS);
        end
    end
end % subject loop

%% correlate each other
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    if ~strcmp(prm_nm,'CID')
        for iNutri = 1:length(nutri_vars)
            nutri_nm = nutri_vars{iNutri};
            % filter bad subjects
            subs_ok = (~isnan(nutri.(nutri_nm)).*~isnan(prm.(prm_nm))) == 1;
            if sum(subs_ok) > 0
                % perform glm brain GSH = b0 + b1*nutrition variable
                [betas.(prm_nm).(['f_',nutri_nm]),~,stats_tmp] = glmfit(nutri.(nutri_nm)(subs_ok), prm.(prm_nm)(subs_ok),'normal');
                % extract p.value
                pval.(prm_nm).(['f_',nutri_nm]) = stats_tmp.p;
                % note which are significant
                if stats_tmp.p(2) < 0.05
                    pval.signif.([prm_nm,'_f_',nutri_nm]) = stats_tmp.p(2);
                end
                % extract the fit
                nutri_bis.(prm_nm).(nutri_nm) = sort(nutri.(nutri_nm)(subs_ok));
                fittedData.(prm_nm).(['f_',nutri_nm]) = glmval(betas.(prm_nm).(['f_',nutri_nm]), nutri_bis.(prm_nm).(nutri_nm), 'identity');
            end
        end % nutrition component
    end % CID field filter
end % parameter

%% figure
if figDisp == 1
    
    lWidth = 3;
    pSize = 25;
    rawDataCol = 'k';
    fitCol = [255 143 143]./255;
    
    %% show results
    for iPrm = 1:n_prm
        prm_nm = parameters{iPrm};
        
        if ~strcmp(prm_nm,'CID')
            % initiliaze figures
            fig1 = fig;
            fig2 = fig;
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
                end
                
                %% raw food
                figure(fig1);
                subplot(1,3,iNutri);
                hold on;
                % show raw data
                curve_hdl = scatter(nutri.(nutri_nm), prm.(prm_nm));
                curve_hdl.LineWidth = lWidth;
                curve_hdl.MarkerEdgeColor = rawDataCol;
                % add the fit
                curve_fit_hdl = plot(nutri_bis.(prm_nm).(nutri_nm),...
                    fittedData.(prm_nm).(['f_',(nutri_nm)]));
                curve_fit_hdl.LineStyle = '--';
                curve_fit_hdl.LineWidth = lWidth;
                curve_fit_hdl.Color = fitCol;
                switch nutri_nm
                    case {'niacin','Trp'}
                        xlabel([nutri_nm,' (mg/week)']);
                    case 'niacinEquiv'
                        xlabel('niacin eq. (mg/week)');
                end
                ylabel(prm_nm);
                legend_size(pSize);
                
                %% Nutrition/total calories
                figure(fig2);
                subplot(1,3,iNutri);
                hold on;
                % show data normalized
                curve_hdl = scatter(nutri.(nutri_nm_bis), prm.(prm_nm));
                curve_hdl.LineWidth = lWidth;
                curve_hdl.MarkerEdgeColor = rawDataCol;
                % add the fit
                curve_fit_hdl = plot(nutri_bis.(prm_nm).(nutri_nm_bis),...
                    fittedData.(prm_nm).(['f_',(nutri_nm_bis)]));
                curve_fit_hdl.LineStyle = '--';
                curve_fit_hdl.LineWidth = lWidth;
                curve_fit_hdl.Color = fitCol;
                if ismember(nutri_nm,{'niacin','Trp'})
                    xlabel([nutri_nm,'/calories']);
                elseif strcmp(nutri_nm,'niacinEquiv')
                    xlabel('niacin eq./calories');
                end
                ylabel(prm_nm);
                legend_size(pSize);
            end % nutrition metabolite
            
            %% add plot with calories only
            fig3 = fig;
            % show raw data
            curve_hdl = scatter(nutri.calories, prm.(prm_nm));
            curve_hdl.LineWidth = lWidth;
            curve_hdl.MarkerEdgeColor = rawDataCol;
            % add the fit
            curve_fit_hdl = plot(nutri_bis.(prm_nm).calories,...
                fittedData.(prm_nm).f_calories);
            curve_fit_hdl.LineStyle = '--';
            curve_fit_hdl.LineWidth = lWidth;
            curve_fit_hdl.Color = fitCol;
            xlabel('calories (kcal/week)');
            ylabel(prm_nm);
            legend_size(pSize);
        end % CID field filter
    end % parameter loop
end % figure display