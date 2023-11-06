function[betas, pval, r_corr] = metabolites_f_sleep()
% [betas, pval, r_corr] = metabolites_f_sleep()
% check correlation between level of metabolites in a given brain area and
% sleep previous day of the experiment.
%
% INPUTS
%
% OUTPUTS
% betas, pval, r_corr: structure with betas, p.values and correlation 
% coefficients r_corr for GLM testing correlation between levels of 
% metabolites (as Y variable) and average sleep/sleep the previous day of 
% the experiment/difference between sleep previous day of the experiment 
% and average amount of sleep.
%
% See also metabolite_load, load_gal_data

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

%% perform correlations
% metabolites = b0+b1*average sleep
goodSubjects = (~isnan(avgSleep).*~isnan(metabolite_allSubs)) == true;
[betas.met_f_avgSleep, ~, stats_met_f_avgSleep] =...
    glmfit(avgSleep(goodSubjects), metabolite_allSubs(goodSubjects), 'normal');
pval.met_f_avgSleep = stats_met_f_avgSleep.p(2);
metabolite_fit_avgSleep = glmval(betas.met_f_avgSleep,...
    avgSleep(goodSubjects), 'identity');
% extract also correlation coefficient
[r_corr.met_f_avgSleep] =...
    glmfit(nanzscore(avgSleep(goodSubjects)),...
    nanzscore(metabolite_allSubs(goodSubjects)), 'normal');

% metabolites = b0+b1*previous day sleep
goodSubjects = (~isnan(prevDaySleep).*~isnan(metabolite_allSubs)) == true;
[betas.met_f_prevDaySleep, ~, stats_met_f_prevDaySleep] =...
    glmfit(prevDaySleep(goodSubjects), metabolite_allSubs(goodSubjects), 'normal');
pval.met_f_prevDaySleep = stats_met_f_prevDaySleep.p(2);
metabolite_fit_prevDaySleep = glmval(betas.met_f_prevDaySleep,...
    prevDaySleep(goodSubjects), 'identity');
% extract also correlation coefficient
[r_corr.met_f_prevDaySleep] =...
    glmfit(nanzscore(prevDaySleep(goodSubjects)),...
    nanzscore(metabolite_allSubs(goodSubjects)), 'normal');

% metabolites = b0+b1*(previous day sleep - average sleep)
goodSubjects = (~isnan(delta_PrevDay_AvgSleep).*~isnan(metabolite_allSubs)) == true;
[betas.met_f_delta_PrevDay_AvgSleep, ~, stats_met_f_delta_PrevDay_AvgSleep] =...
    glmfit(delta_PrevDay_AvgSleep(goodSubjects), metabolite_allSubs(goodSubjects), 'normal');
pval.met_f_delta_PrevDay_AvgSleep = stats_met_f_delta_PrevDay_AvgSleep.p(2);
metabolite_fit_delta_PrevDay_AvgSleep = glmval(betas.met_f_delta_PrevDay_AvgSleep,...
    delta_PrevDay_AvgSleep(goodSubjects), 'identity');
% extract also correlation coefficient
[r_corr.met_f_delta_PrevDay_AvgSleep] =...
    glmfit(nanzscore(delta_PrevDay_AvgSleep(goodSubjects)),...
    nanzscore(metabolite_allSubs(goodSubjects)), 'normal');

%% display figures
mSize = 100;
lWidth = 3;
lWidthScatter = 1.5;
blackCol = [0 0 0];
greyCol = [143 143 143]./255;
pSize = 30;
yLabeling = [MRS_ROI_nm,' - ',metabolite_nm_bis];

% metabolite = f(avg sleep)
fig;
hold on;
hdl = scatter(hours(minutes(avgSleep)), metabolite_allSubs,...
    'SizeData',mSize,'LineWidth',lWidthScatter,...
    'MarkerEdgeColor',blackCol,'MarkerFaceColor',greyCol);
plot(hours(minutes(avgSleep(goodSubjects))), metabolite_fit_avgSleep,...
    'LineWidth',lWidth,'LineStyle','-','Color',blackCol);
ylabel(yLabeling);
xlabel('average sleep (h)');
legend_size(pSize);

% metabolite = f(previous day sleep)
fig;
hold on;
scatter(hours(minutes(prevDaySleep)), metabolite_allSubs,...
    'SizeData',mSize,'LineWidth',lWidthScatter,...
    'MarkerEdgeColor',blackCol,'MarkerFaceColor',greyCol);
plot(hours(minutes(prevDaySleep(goodSubjects))),...
    metabolite_fit_prevDaySleep,...
    'LineWidth',lWidth,'LineStyle','-','Color',blackCol);
ylabel(yLabeling);
xlabel('previous day sleep (h)');
legend_size(pSize);

% metabolite = f(avg sleep)
fig;
hold on;
scatter(hours(minutes(delta_PrevDay_AvgSleep)), metabolite_allSubs,...
    'SizeData',mSize,'LineWidth',lWidthScatter,...
    'MarkerEdgeColor',blackCol,'MarkerFaceColor',greyCol);
plot(hours(minutes(delta_PrevDay_AvgSleep(goodSubjects))),...
    metabolite_fit_delta_PrevDay_AvgSleep,...
    'LineWidth',lWidth,'LineStyle','-','Color',blackCol);
ylabel(yLabeling);
xlabel('previous day - avg sleep (h)');
legend_size(pSize);


end % function