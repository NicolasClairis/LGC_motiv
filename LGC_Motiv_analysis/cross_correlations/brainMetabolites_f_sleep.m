function[betas, pval, r_corr] = metabolites_f_sleep(fig_disp, rmv_outliers_yn)
% [betas, pval, r_corr] = metabolites_f_sleep(fig_disp, rmv_outliers_yn)
% check correlation between level of metabolites in a given brain area and
% sleep previous day of the experiment.
%
% INPUTS
% fig_disp: display figure (1) or not (0)? By default will display it
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% betas, pval, r_corr: structure with betas, p.values and correlation 
% coefficients r_corr for GLM testing correlation between levels of 
% metabolites (as Y variable) and average sleep/sleep the previous day of 
% the experiment/difference between sleep previous day of the experiment 
% and average amount of sleep.
%
% See also metabolite_load, load_gal_data

%% default inputs
if ~exist('fig_disp','var') || isempty(fig_disp) || ~ismember(fig_disp,[0,1])
    fig_disp = 1;
end
% remove outliers by default
if ~exist('rmv_outliers_yn','var') || isempty(rmv_outliers_yn) || ~ismember(rmv_outliers_yn,[0,1])
    rmv_outliers_yn = 1;
end

%% define subjects to check
condition = subject_condition;
study_nm = 'study1';
[subject_id, NS] = LGCM_subject_selection(study_nm,condition,'all');

%% define metabolite of interest
[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);
[metabolite_nm_bis] = metab_div_rnm(metabolite_nm);
 
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
for iS = 1:NS
    sub_nm = subject_id{iS};

    %% load sleep
    % loop through subjects to find the correct one
    for iLine = 1:size(sleepCID,1)
        isItThisSubject = strcmp(sub_nm, sleepCID{iLine}(4:6));
        if isItThisSubject == true
            sub_idx = iLine;
        end
    end
    avgSleep(iS) = minutes(duration(avgSleepTable{sub_idx},'InputFormat','hh:mm'));
    prevDaySleep(iS) = minutes(duration(prevDaySleepTable{sub_idx},'InputFormat','hh:mm'));
end % subject loop

% look also at the difference between previous day and average
delta_PrevDay_AvgSleep = prevDaySleep - avgSleep;

%% remove outliers
if rmv_outliers_yn == 1
    [~,~,avgSleep] = rmv_outliers_3sd(avgSleep);
    [~,~,prevDaySleep] = rmv_outliers_3sd(prevDaySleep);
    [~,~,delta_PrevDay_AvgSleep] = rmv_outliers_3sd(delta_PrevDay_AvgSleep);
    [~,~,metabolite_allSubs] = rmv_outliers_3sd(metabolite_allSubs);
end % remove outliers

%% perform correlations
% metabolites = b0+b1*average sleep
goodSubjects_avgS = (~isnan(avgSleep).*~isnan(metabolite_allSubs)) == true;
[betas.met_f_avgSleep, ~, stats_met_f_avgSleep] =...
    glmfit(avgSleep(goodSubjects_avgS), metabolite_allSubs(goodSubjects_avgS), 'normal');
pval.met_f_avgSleep = stats_met_f_avgSleep.p(2);
metabolite_fit_avgSleep = glmval(betas.met_f_avgSleep,...
    avgSleep(goodSubjects_avgS), 'identity');
% extract also correlation coefficient
[r_corr.met_f_avgSleep] =...
    glmfit(nanzscore(avgSleep(goodSubjects_avgS)),...
    nanzscore(metabolite_allSubs(goodSubjects_avgS)), 'normal');

% metabolites = b0+b1*previous day sleep
goodSubjects_prevDayS = (~isnan(prevDaySleep).*~isnan(metabolite_allSubs)) == true;
[betas.met_f_prevDaySleep, ~, stats_met_f_prevDaySleep] =...
    glmfit(prevDaySleep(goodSubjects_prevDayS), metabolite_allSubs(goodSubjects_prevDayS), 'normal');
pval.met_f_prevDaySleep = stats_met_f_prevDaySleep.p(2);
metabolite_fit_prevDaySleep = glmval(betas.met_f_prevDaySleep,...
    prevDaySleep(goodSubjects_prevDayS), 'identity');
% extract also correlation coefficient
[r_corr.met_f_prevDaySleep] =...
    glmfit(nanzscore(prevDaySleep(goodSubjects_prevDayS)),...
    nanzscore(metabolite_allSubs(goodSubjects_prevDayS)), 'normal');

% metabolites = b0+b1*(previous day sleep - average sleep)
goodSubjects_deltaS = (~isnan(delta_PrevDay_AvgSleep).*~isnan(metabolite_allSubs)) == true;
[betas.met_f_delta_PrevDay_AvgSleep, ~, stats_met_f_delta_PrevDay_AvgSleep] =...
    glmfit(delta_PrevDay_AvgSleep(goodSubjects_deltaS), metabolite_allSubs(goodSubjects_deltaS), 'normal');
pval.met_f_delta_PrevDay_AvgSleep = stats_met_f_delta_PrevDay_AvgSleep.p(2);
metabolite_fit_delta_PrevDay_AvgSleep = glmval(betas.met_f_delta_PrevDay_AvgSleep,...
    delta_PrevDay_AvgSleep(goodSubjects_deltaS), 'identity');
% extract also correlation coefficient
[r_corr.met_f_delta_PrevDay_AvgSleep] =...
    glmfit(nanzscore(delta_PrevDay_AvgSleep(goodSubjects_deltaS)),...
    nanzscore(metabolite_allSubs(goodSubjects_deltaS)), 'normal');

%% display figures
if fig_disp == 1
    mSize = 100;
    lWidth = 3;
    lWidthScatter = 1.5;
    blackCol = [0 0 0];
    greyCol = [143 143 143]./255;
    yLabeling = [MRS_ROI_nm,' ',metabolite_nm_bis];
    % yLabeling = 'dmPFC/dACC lactate (mM)';
    % yLabeling = 'aIns lactate (mM)';
    
    % metabolite = f(avg sleep)
    fig;
    hold on;
    scatter(hours(minutes(avgSleep)), metabolite_allSubs,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',blackCol,'MarkerFaceColor',greyCol);
    plot(hours(minutes(avgSleep(goodSubjects_avgS))), metabolite_fit_avgSleep,...
        'LineWidth',lWidth,'LineStyle','-','Color',blackCol);
    ylabel(yLabeling);
    xlabel('average sleep (h)');
    place_r_and_pval(r_corr.met_f_avgSleep(2), pval.met_f_avgSleep);
    
    % metabolite = f(previous day sleep)
    fig;
    hold on;
    scatter(hours(minutes(prevDaySleep)), metabolite_allSubs,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',blackCol,'MarkerFaceColor',greyCol);
    plot(hours(minutes(prevDaySleep(goodSubjects_prevDayS))),...
        metabolite_fit_prevDaySleep,...
        'LineWidth',lWidth,'LineStyle','-','Color',blackCol);
    ylabel(yLabeling);
    xlabel('previous day sleep (h)');
    place_r_and_pval(r_corr.met_f_prevDaySleep(2), pval.met_f_prevDaySleep);
    
    % metabolite = f(avg sleep)
    fig;
    hold on;
    scatter(hours(minutes(delta_PrevDay_AvgSleep)), metabolite_allSubs,...
        'SizeData',mSize,'LineWidth',lWidthScatter,...
        'MarkerEdgeColor',blackCol,'MarkerFaceColor',greyCol);
    plot(hours(minutes(delta_PrevDay_AvgSleep(goodSubjects_deltaS))),...
        metabolite_fit_delta_PrevDay_AvgSleep,...
        'LineWidth',lWidth,'LineStyle','-','Color',blackCol);
    ylabel(yLabeling);
    xlabel('previous day - avg sleep (h)');
    place_r_and_pval(r_corr.met_f_delta_PrevDay_AvgSleep(2), pval.met_f_delta_PrevDay_AvgSleep);
end % figure display

end % function