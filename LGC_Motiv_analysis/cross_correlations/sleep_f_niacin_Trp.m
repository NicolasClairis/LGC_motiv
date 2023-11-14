function[betas, pval, r_corr, pval_signif] = sleep_f_niacin_Trp()
% [betas, pval, r_corr, pval_signif] = sleep_f_niacin_Trp()
% check correlation between sleep duration and nutrition intake of NAD
% precursors.
%
% INPUTS
%
% OUTPUTS
% betas, pval, r_corr: structure with betas, p.values and correlation 
% coefficients r_corr for GLM testing correlation between levels of 
% metabolites (as X variable) and average sleep/sleep the previous day of 
% the experiment/difference between sleep previous day of the experiment 
% and average amount of sleep (as Y variable).
%
% pval_signif: structure with significant p.values and correlation
% coefficients (for p<0.05)

%% define subjects to check
condition = subject_condition;
study_nm = 'study1';
[subject_id, NS] = LGCM_subject_selection(study_nm,condition,'all');
%% working directories
root = LGCM_root_paths;
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
%% initialize variables of interest
[nutri.niacin, nutri.Trp, nutri.niacinEquiv, nutri.niacinPlusTrp,...
    nutri.calories,...
    nutri.niacin_div_totalCal,  nutri.Trp_div_totalCal,...
    nutri.niacinEquiv_div_totalCal,...
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
%% load sleep
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
prevDaySleepTable = excelReadGeneralFile.HeuresDeSommeilLaVeilleDeL_exp_rience;
avgSleepTable = excelReadGeneralFile.HeuresDeSommeil_enMoyenne_;
sleepCID = excelReadGeneralFile.CID;
% replace 'h' of hours by ':'
prevDaySleepTable = strrep(prevDaySleepTable, 'h',':');
avgSleepTable = strrep(avgSleepTable, 'h',':');
% extract the data
[avgSleep, prevDaySleep] = deal(NaN(1,NS));

%% loop through subjects to extract nutrition + sleep
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
    
    %% load sleep
    % loop through subjects to find the correct one
    for iLine = 1:size(sleepCID,1)
        isItThisSubject = strcmp(sub_nm, sleepCID{iLine}(4:6));
        if isItThisSubject == true
            sub_sleep_idx = iLine;
        end
    end
    avgSleep(iS) = minutes(duration(avgSleepTable{sub_sleep_idx},'InputFormat','hh:mm'));
    prevDaySleep(iS) = minutes(duration(prevDaySleepTable{sub_sleep_idx},'InputFormat','hh:mm'));
end % subject loop

% look also at the difference between previous day and average
delta_PrevDay_AvgSleep = prevDaySleep - avgSleep;

%% perform correlations
for iNutri = 1:length(nutri_vars)
    nutri_nm = nutri_vars{iNutri};
    
    %% average sleep
    % filter bad subjects
    subs_avgS_ok = (~isnan(nutri.(nutri_nm)).*~isnan(avgSleep)) == 1;
    if sum(subs_avgS_ok) > 0
        % perform glm sleep = b0 + b1*nutrition variable
        [betas.avgSleep.(['f_',nutri_nm]),~,stats_tmp] = glmfit(nutri.(nutri_nm)(subs_avgS_ok), avgSleep(subs_avgS_ok),'normal');
        [r_corr.avgSleep.(['f_',nutri_nm])] = glmfit(zscore(nutri.(nutri_nm)(subs_avgS_ok)), zscore(avgSleep(subs_avgS_ok)),'normal');
        % extract p.value
        pval.avgSleep.(['f_',nutri_nm]) = stats_tmp.p;
        if stats_tmp.p(2) < 0.05
            pval_signif.avgSleep.(['f_',nutri_nm]).pval = stats_tmp.p(2);
            pval_signif.avgSleep.(['f_',nutri_nm]).r_corr = r_corr.avgSleep.(['f_',nutri_nm])(2);
        end
        % extract the fit
        nutri_bis.avgSleep.(nutri_nm) = sort(nutri.(nutri_nm)(subs_avgS_ok));
        fittedData.avgSleep.(['f_',nutri_nm]) = glmval(betas.avgSleep.(['f_',nutri_nm]), nutri_bis.avgSleep.(nutri_nm), 'identity');
    end
    
    %% previous day sleep
    subs_prevS_ok = (~isnan(nutri.(nutri_nm)).*~isnan(prevDaySleep)) == 1;
    if sum(subs_prevS_ok) > 0
        % perform glm sleep = b0 + b1*nutrition variable
        [betas.prevDaySleep.(['f_',nutri_nm]),~,stats_tmp] = glmfit(nutri.(nutri_nm)(subs_prevS_ok), prevDaySleep(subs_prevS_ok),'normal');
        [r_corr.prevDaySleep.(['f_',nutri_nm])] = glmfit(zscore(nutri.(nutri_nm)(subs_prevS_ok)), zscore(prevDaySleep(subs_prevS_ok)),'normal');
        % extract p.value
        pval.prevDaySleep.(['f_',nutri_nm]) = stats_tmp.p;
        if stats_tmp.p(2) < 0.05
            pval_signif.prevDaySleep.(['f_',nutri_nm]).pval = stats_tmp.p(2);
            pval_signif.prevDaySleep.(['f_',nutri_nm]).r_corr = r_corr.prevDaySleep.(['f_',nutri_nm])(2);
        end
        % extract the fit
        nutri_bis.prevDaySleep.(nutri_nm) = sort(nutri.(nutri_nm)(subs_prevS_ok));
        fittedData.prevDaySleep.(['f_',nutri_nm]) = glmval(betas.prevDaySleep.(['f_',nutri_nm]), nutri_bis.prevDaySleep.(nutri_nm), 'identity');
    end
    
    %% (previous day) - (average) sleep
    subs_deltaS_ok = (~isnan(nutri.(nutri_nm)).*~isnan(delta_PrevDay_AvgSleep)) == 1;
    if sum(subs_deltaS_ok) > 0
        % perform glm sleep = b0 + b1*nutrition variable
        [betas.delta_PrevDay_AvgSleep.(['f_',nutri_nm]),~,stats_tmp] = glmfit(nutri.(nutri_nm)(subs_deltaS_ok), delta_PrevDay_AvgSleep(subs_deltaS_ok),'normal');
        [r_corr.delta_PrevDay_AvgSleep.(['f_',nutri_nm])] = glmfit(zscore(nutri.(nutri_nm)(subs_deltaS_ok)), zscore(delta_PrevDay_AvgSleep(subs_deltaS_ok)),'normal');
        % extract p.value
        pval.delta_PrevDay_AvgSleep.(['f_',nutri_nm]) = stats_tmp.p;
        if stats_tmp.p(2) < 0.05
            pval_signif.delta_PrevDay_AvgSleep.(['f_',nutri_nm]).pval = stats_tmp.p(2);
            pval_signif.delta_PrevDay_AvgSleep.(['f_',nutri_nm]).r_corr = r_corr.delta_PrevDay_AvgSleep.(['f_',nutri_nm])(2);
        end
        % extract the fit
        nutri_bis.delta_PrevDay_AvgSleep.(nutri_nm) = sort(nutri.(nutri_nm)(subs_deltaS_ok));
        fittedData.delta_PrevDay_AvgSleep.(['f_',nutri_nm]) = glmval(betas.delta_PrevDay_AvgSleep.(['f_',nutri_nm]), nutri_bis.delta_PrevDay_AvgSleep.(nutri_nm), 'identity');
    end
end % loop through nutrition variables

%% show results
[pSize, lWidth, col, mSize] = general_fig_prm;
avgS_fig = fig;
prevS_fig = fig;
deltaS_fig = fig;
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
    
    %% average sleep figure
    figure(avgS_fig);
    % raw data
    subplot(2,3,iNutri);
    hold on;
    % show raw data
    avgSleep_hdl = scatter(nutri.(nutri_nm), avgSleep);
    avgSleep_hdl.LineWidth = lWidth;
    avgSleep_hdl.MarkerEdgeColor = 'k';
    % add the fit
    avgS_fit_hdl = plot(nutri_bis.avgSleep.(nutri_nm),...
        fittedData.avgSleep.(['f_',nutri_nm]));
    avgS_fit_hdl.LineStyle = '--';
    avgS_fit_hdl.LineWidth = lWidth;
    avgS_fit_hdl.Color = col.grey;
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
    ylabel('average Sleep (min)');
    legend_size(pSize);
    
    % Nutrition/total calories
    subplot(2,3,iNutri+3);
    hold on;
    % show data normalized
    avgSleep_norm_hdl = scatter(nutri.(nutri_nm_bis), avgSleep);
    avgSleep_norm_hdl.LineWidth = lWidth;
    avgSleep_norm_hdl.MarkerEdgeColor = 'k';
    % add the fit
    avgS_norm_fit_hdl = plot(nutri_bis.avgSleep.(nutri_nm_bis),...
        fittedData.avgSleep.(['f_',nutri_nm_bis]));
    avgS_norm_fit_hdl.LineStyle = '--';
    avgS_norm_fit_hdl.LineWidth = lWidth;
    avgS_norm_fit_hdl.Color = col.grey;
    
    switch nutri_nm
        case {'niacin','Trp'}
            xlabel([nutri_nm,'/calories']);
        case 'niacinPlusTrp'
            xlabel('(niacin + Trp)/calories');
        case 'niacinEquiv'
            xlabel('niacin equivalents/calories (mg/week)');
    end
    ylabel('average Sleep (min)');
    legend_size(pSize);
    
    %% previous day sleep figure
    figure(prevS_fig);
    % raw data
    subplot(2,3,iNutri);
    hold on;
    % show raw data
    prevDaySleep_hdl = scatter(nutri.(nutri_nm), prevDaySleep);
    prevDaySleep_hdl.LineWidth = lWidth;
    prevDaySleep_hdl.MarkerEdgeColor = 'k';
    % add the fit
    avgS_fit_hdl = plot(nutri_bis.prevDaySleep.(nutri_nm),...
        fittedData.prevDaySleep.(['f_',nutri_nm]));
    avgS_fit_hdl.LineStyle = '--';
    avgS_fit_hdl.LineWidth = lWidth;
    avgS_fit_hdl.Color = col.grey;
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
    ylabel('previous day Sleep (min)');
    legend_size(pSize);
    
    % Nutrition/total calories
    subplot(2,3,iNutri+3);
    hold on;
    % show data normalized
    prevDaySleep_norm_hdl = scatter(nutri.(nutri_nm_bis), prevDaySleep);
    prevDaySleep_norm_hdl.LineWidth = lWidth;
    prevDaySleep_norm_hdl.MarkerEdgeColor = 'k';
    % add the fit
    avgS_norm_fit_hdl = plot(nutri_bis.prevDaySleep.(nutri_nm_bis),...
        fittedData.prevDaySleep.(['f_',nutri_nm_bis]));
    avgS_norm_fit_hdl.LineStyle = '--';
    avgS_norm_fit_hdl.LineWidth = lWidth;
    avgS_norm_fit_hdl.Color = col.grey;
    
    switch nutri_nm
        case {'niacin','Trp'}
            xlabel([nutri_nm,'/calories']);
        case 'niacinPlusTrp'
            xlabel('(niacin + Trp)/calories');
        case 'niacinEquiv'
            xlabel('niacin equivalents/calories (mg/week)');
    end
    ylabel('previous day Sleep (min)');
    legend_size(pSize);
    
    %% delta previous day - average sleep figure
    figure(deltaS_fig);
    % raw data
    subplot(2,3,iNutri);
    hold on;
    % show raw data
    delta_PrevDay_AvgSleep_hdl = scatter(nutri.(nutri_nm), delta_PrevDay_AvgSleep);
    delta_PrevDay_AvgSleep_hdl.LineWidth = lWidth;
    delta_PrevDay_AvgSleep_hdl.MarkerEdgeColor = 'k';
    % add the fit
    avgS_fit_hdl = plot(nutri_bis.delta_PrevDay_AvgSleep.(nutri_nm),...
        fittedData.delta_PrevDay_AvgSleep.(['f_',nutri_nm]));
    avgS_fit_hdl.LineStyle = '--';
    avgS_fit_hdl.LineWidth = lWidth;
    avgS_fit_hdl.Color = col.grey;
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
    ylabel('previous day - average Sleep (min)');
    legend_size(pSize);
    
    % Nutrition/total calories
    subplot(2,3,iNutri+3);
    hold on;
    % show data normalized
    delta_PrevDay_AvgSleep_norm_hdl = scatter(nutri.(nutri_nm_bis), delta_PrevDay_AvgSleep);
    delta_PrevDay_AvgSleep_norm_hdl.LineWidth = lWidth;
    delta_PrevDay_AvgSleep_norm_hdl.MarkerEdgeColor = 'k';
    % add the fit
    avgS_norm_fit_hdl = plot(nutri_bis.delta_PrevDay_AvgSleep.(nutri_nm_bis),...
        fittedData.delta_PrevDay_AvgSleep.(['f_',nutri_nm_bis]));
    avgS_norm_fit_hdl.LineStyle = '--';
    avgS_norm_fit_hdl.LineWidth = lWidth;
    avgS_norm_fit_hdl.Color = col.grey;
    
    switch nutri_nm
        case {'niacin','Trp'}
            xlabel([nutri_nm,'/calories']);
        case 'niacinPlusTrp'
            xlabel('(niacin + Trp)/calories');
        case 'niacinEquiv'
            xlabel('niacin equivalents/calories (mg/week)');
    end
    ylabel('previous day - average Sleep (min)');
    legend_size(pSize);
end % nutrition metabolite


end % function