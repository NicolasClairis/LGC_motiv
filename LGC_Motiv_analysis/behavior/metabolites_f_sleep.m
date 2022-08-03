function[betas, pval] = metabolites_f_sleep()
% check correlation between level of metabolites in a given brain area and
% sleep previous day of the experiment.
%
% INPUTS
%
% OUTPUTS
% betas, pval: structure with betas and p.values for GLM testing
% correlation between levels of metabolites (as Y variable) and average
% sleep/sleep the previous day of the experiment/difference between sleep
% previous day of the experiment and average amount of sleep
%
% See also metabolite_load, load_gal_data

%% define subjects to check
condition = subject_condition;
study_nm = 'study1';
[subject_id, NS] = LGCM_subject_selection(study_nm,condition,'all');

%% define metabolite of interest
switch study_nm
    case 'study1'
        %% which ROI?
        ROIs = {'dmPFC','aIns'};
        nROIs = length(ROIs);
        ROI_idx = spm_input('Metabolites in which brain area?',1,'m',...
            ROIs,1:nROIs,0);
        ROI_nm = ROIs{ROI_idx};
        %% select metabolite of interest
        metabolites = {'Mac','Ala','Asp','PCho','Cr','PCr','GABA',...
            'Gln','Glu','GSH','Gly','Ins','Lac','NAA','Scyllo','Tau',...
            'Asc','Glc','NAAG','GPC','PE','Ser',...
            'NAA_NAAG','Glu_Gln','GPC_PCho','Cr_PCr','Gly_Ins','Gln_div_Glu'};
    otherwise
        error(['not ready yet for ',study_nm]);
end
n_met = length(metabolites);
metabolite_idx = spm_input('Which metabolite to focus on?',1,'m',...
    metabolites,1:n_met,0);
metabolite_nm = metabolites{metabolite_idx};

%% load sleep
[excelReadGeneralFile] = load_gal_data();
prevDaySleepTable = excelReadGeneralFile.HeureDeSommeilLaVeilleDeL_exp_rience;
avgSleepTable = excelReadGeneralFile.HeuresDeSommeil_enMoyenne_;
sleepCID = excelReadGeneralFile.IDParticipant;
% replace 'h' of hours by ':'
prevDaySleepTable = strrep(prevDaySleepTable, 'h',':');
avgSleepTable = strrep(avgSleepTable, 'h',':');
% extract the data
[avgSleep, prevDaySleep, metabolite] = deal(NaN(1,NS));
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
    if ~strcmp(avgSleepTable{sub_idx},'NaN') && ~strcmp(avgSleepTable{sub_idx},'8:30 (visit 1) 8:30 (visit 2)')
        avgSleep(iS) = minutes(duration(avgSleepTable{sub_idx},'InputFormat','hh:mm'));
    elseif strcmp(avgSleepTable{sub_idx},'8:30 (visit 1) 8:30 (visit 2)') % both equal
        avgSleep(iS) = minutes(duration('8:30','InputFormat','hh:mm'));
    end
    if ~strcmp(prevDaySleepTable{sub_idx},'NaN') && ~strcmp(prevDaySleepTable{sub_idx},'9:00 (visit 1) 8:00 (visit 2)')
        prevDaySleep(iS) = minutes(duration(prevDaySleepTable{sub_idx},'InputFormat','hh:mm'));
    elseif strcmp(prevDaySleepTable{sub_idx},'9:00 (visit 1) 8:00 (visit 2)') % focus on 2nd visit
        prevDaySleep(iS) = minutes(duration('8:00','InputFormat','hh:mm'));
    end

    %% load metabolite
    [metabolites] = metabolite_load({sub_nm});
    % focus on metabolite and brain area selected
    metabolite(iS) = metabolites.(ROI_nm).(metabolite_nm);
end % subject loop

% look also at the difference between previous day and average
delta_PrevDay_AvgSleep = prevDaySleep - avgSleep;

%% perform correlations
% metabolites = b0+b1*average sleep
goodSubjects = (~isnan(avgSleep).*~isnan(metabolite)) == true;
[betas.met_f_avgSleep, ~, stats_met_f_avgSleep] =...
    glmfit(avgSleep(goodSubjects), metabolite(goodSubjects), 'normal');
pval.met_f_avgSleep = stats_met_f_avgSleep.p(2);
metabolite_fit_avgSleep = glmval(betas.met_f_avgSleep,...
    avgSleep(goodSubjects), 'identity');
% metabolites = b0+b1*previous day sleep
goodSubjects = (~isnan(prevDaySleep).*~isnan(metabolite)) == true;
[betas.met_f_prevDaySleep, ~, stats_met_f_prevDaySleep] =...
    glmfit(prevDaySleep(goodSubjects), metabolite(goodSubjects), 'normal');
pval.met_f_prevDaySleep = stats_met_f_prevDaySleep.p(2);
metabolite_fit_prevDaySleep = glmval(betas.met_f_prevDaySleep,...
    prevDaySleep(goodSubjects), 'identity');
% metabolites = b0+b1*(previous day sleep - average sleep)
goodSubjects = (~isnan(delta_PrevDay_AvgSleep).*~isnan(metabolite)) == true;
[betas.met_f_delta_PrevDay_AvgSleep, ~, stats_met_f_delta_PrevDay_AvgSleep] =...
    glmfit(delta_PrevDay_AvgSleep(goodSubjects), metabolite(goodSubjects), 'normal');
pval.met_f_delta_PrevDay_AvgSleep = stats_met_f_delta_PrevDay_AvgSleep.p(2);
metabolite_fit_delta_PrevDay_AvgSleep = glmval(betas.met_f_delta_PrevDay_AvgSleep,...
    delta_PrevDay_AvgSleep(goodSubjects), 'identity');


%% display figures
mSize = 100;
lWidth = 3;
lWidthScatter = 1.5;
blackCol = [0 0 0];
greyCol = [143 143 143]./255;
pSize = 30;

% metabolite = f(avg sleep)
fig;
hold on;
hdl = scatter(hours(minutes(avgSleep)), metabolite,...
    'SizeData',mSize,'LineWidth',lWidthScatter,...
    'MarkerEdgeColor',blackCol,'MarkerFaceColor',greyCol);
plot(hours(minutes(avgSleep(goodSubjects))), metabolite_fit_avgSleep,...
    'LineWidth',lWidth,'LineStyle','-','Color',blackCol);
ylabel(metabolite_nm);
xlabel('average sleep (h)');
legend_size(pSize);

% metabolite = f(previous day sleep)
fig;
hold on;
scatter(hours(minutes(prevDaySleep)), metabolite,...
    'SizeData',mSize,'LineWidth',lWidthScatter,...
    'MarkerEdgeColor',blackCol,'MarkerFaceColor',greyCol);
plot(hours(minutes(prevDaySleep(goodSubjects))),...
    metabolite_fit_prevDaySleep,...
    'LineWidth',lWidth,'LineStyle','-','Color',blackCol);
ylabel(metabolite_nm);
xlabel('previous day sleep (h)');
legend_size(pSize);

% metabolite = f(avg sleep)
fig;
hold on;
scatter(hours(minutes(delta_PrevDay_AvgSleep)), metabolite,...
    'SizeData',mSize,'LineWidth',lWidthScatter,...
    'MarkerEdgeColor',blackCol,'MarkerFaceColor',greyCol);
plot(hours(minutes(delta_PrevDay_AvgSleep(goodSubjects))),...
    metabolite_fit_delta_PrevDay_AvgSleep,...
    'LineWidth',lWidth,'LineStyle','-','Color',blackCol);
ylabel(metabolite_nm);
xlabel('previous day - avg sleep (h)');
legend_size(pSize);


end % function