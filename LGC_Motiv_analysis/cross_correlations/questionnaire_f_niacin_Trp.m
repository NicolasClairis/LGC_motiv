
%% test correlation between some questionnaires and nutrition components
% extracted through the FFQ questionnaire

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
        gitPath = fullfile('C:','Users','Loco','Downloads');
    otherwise
        error('case not ready yet');
end
nutritionPath = [fullfile(gitPath,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'nutrition'),filesep];

%% load nutrition
[nutri.niacin, nutri.Trp, nutri.niacinPlusTrp,...
    nutri.calories,...
    nutri.niacin_div_totalCal,  nutri.Trp_div_totalCal,...
    nutri.niacinPlusTrp_div_totalCal] = deal(NaN(1,NS));
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

%% load questionnaire
potentialQuestionnaires = {'JPI_RScore', 'MADRS_SCorrected',...
    'MPSTEFSPhysicalTraitScore','MPSTEFSMentalTraitScore',...
    'PunishmentScore','RewardScore',...
    'IPAQ','IPAQInactivity'};
selectedCon = listdlg('PromptString','Which questionnaires?',...
    'ListString',potentialQuestionnaires);
nQuest = length(selectedCon);
questToCheck = potentialQuestionnaires(selectedCon);
[excelReadQuestionnairesFile, sub_CID_list] = load_questionnaires_data();
for iQuest = 1:nQuest
    quest_nm = questToCheck{iQuest};
    quest_data.(quest_nm) = NaN(1,NS);
end

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    %% load questionnaire
    sub_quest_idx = strcmp(sub_CID_list, sub_nm);
    for iQuest = 1:nQuest
        quest_nm = questToCheck{iQuest};
        quest_data.(quest_nm)(iS) = excelReadQuestionnairesFile.(quest_nm)(sub_quest_idx);
    end % questionnaire loop
    
    %% load nutrition score
    sub_niacin_idx = find(strcmp(niacin_table.CID, sub_nm));
    sub_Trp_idx = find(strcmp(Trp_table.CID, sub_nm));
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
for iQuest = 1:nQuest
    quest_nm = questToCheck{iQuest};
    for iNutri = 1:length(nutri_vars)
        nutri_nm = nutri_vars{iNutri};
        % filter bad subjects
        subs_ok = (~isnan(nutri.(nutri_nm)).*~isnan(quest_data.(quest_nm))) == 1;
        if sum(subs_ok) > 0
            % perform glm brain GSH = b0 + b1*nutrition variable
            [betas.(quest_nm).(['f_',nutri_nm]),~,stats_tmp] = glmfit(nutri.(nutri_nm)(subs_ok), quest_data.(quest_nm)(subs_ok),'normal');
            % extract p.value
            pval.(quest_nm).(['f_',nutri_nm]) = stats_tmp.p;
            % extract the fit
            nutri_bis.(quest_nm).(nutri_nm) = sort(nutri.(nutri_nm)(subs_ok));
            fittedData.(quest_nm).(['f_',nutri_nm]) = glmval(betas.(quest_nm).(['f_',nutri_nm]), nutri_bis.(quest_nm).(nutri_nm), 'identity');
        end
    end % loop through nutritional intake variables
end % loop through ROIs

%% figure
if figDisp == 1
    
    lWidth = 3;
    pSize = 25;
    rawDataCol = 'k';
    fitCol = [255 143 143]./255;
    %% show results
    for iQuest = 1:nQuest
        quest_nm = questToCheck{iQuest};
        % define y.label for questionnaires
        switch quest_nm
            case 'JPI_RScore'
                quest_nm_bis = 'JPI-R';
            case 'MADRS_SCorrected'
                quest_nm_bis = 'MADRS-S';
            case 'MPSTEFSPhysicalTraitScore'
                quest_nm_bis = 'MPSTEFS physical';
            case 'MPSTEFSMentalTraitScore'
                quest_nm_bis = 'MPSTEFS mental';
            case 'PunishmentScore'
                quest_nm_bis = 'PANAS P';
            case 'RewardScore'
                quest_nm_bis = 'PANAS R';
            case 'IPAQ'
                quest_nm_bis = 'IPAQ activity';
            case 'IPAQInactivity'
                quest_nm_bis = 'IPAQ inactivity';
        end
        
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
                    nutri_nm = 'niacinPlusTrp';
                    nutri_nm_bis = 'niacinPlusTrp_div_totalCal';
            end
            
            %% raw metabolite
            figure(fig1);
            subplot(1,3,iNutri);
            hold on;
            % show raw data
            curve_hdl = scatter(nutri.(nutri_nm), quest_data.(quest_nm));
            curve_hdl.LineWidth = lWidth;
            curve_hdl.MarkerEdgeColor = rawDataCol;
            % add the fit
            curve_fit_hdl = plot(nutri_bis.(quest_nm).(nutri_nm),...
                fittedData.(quest_nm).(['f_',(nutri_nm)]));
            curve_fit_hdl.LineStyle = '--';
            curve_fit_hdl.LineWidth = lWidth;
            curve_fit_hdl.Color = fitCol;
            switch nutri_nm
                case {'niacin','Trp'}
                    xlabel([nutri_nm,' (mg/week)']);
                case 'niacinPlusTrp'
                    xlabel('niacin + Trp (mg/week)');
            end
            ylabel(quest_nm_bis);
            legend_size(pSize);
            
            %% Nutrition/total calories
            figure(fig2);
            subplot(1,3,iNutri);
            hold on;
            % show data normalized
            curve_hdl = scatter(nutri.(nutri_nm_bis), quest_data.(quest_nm));
            curve_hdl.LineWidth = lWidth;
            curve_hdl.MarkerEdgeColor = rawDataCol;
            % add the fit
            curve_fit_hdl = plot(nutri_bis.(quest_nm).(nutri_nm_bis),...
                fittedData.(quest_nm).(['f_',(nutri_nm_bis)]));
            curve_fit_hdl.LineStyle = '--';
            curve_fit_hdl.LineWidth = lWidth;
            curve_fit_hdl.Color = fitCol;
            if ismember(nutri_nm,{'niacin','Trp'})
                xlabel([nutri_nm,'/calories']);
            elseif strcmp(nutri_nm,'niacinPlusTrp')
                xlabel('(niacin + Trp)/calories');
            end
            ylabel(quest_nm_bis);
            legend_size(pSize);
        end % nutrition metabolite
        
        %% add plot with calories only
        fig3 = fig;
        % show raw data
        curve_hdl = scatter(nutri.calories, quest_data.(quest_nm));
        curve_hdl.LineWidth = lWidth;
        curve_hdl.MarkerEdgeColor = rawDataCol;
        % add the fit
        curve_fit_hdl = plot(nutri_bis.(quest_nm).calories,...
            fittedData.(quest_nm).f_calories);
        curve_fit_hdl.LineStyle = '--';
        curve_fit_hdl.LineWidth = lWidth;
        curve_fit_hdl.Color = fitCol;
        xlabel('calories (kcal/week)');
        ylabel(quest_nm_bis);
        legend_size(pSize);
    end % questionnaire loop
end % figure display