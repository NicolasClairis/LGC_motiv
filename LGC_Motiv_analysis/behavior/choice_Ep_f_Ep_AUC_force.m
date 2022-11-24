% function[] = choice_Ep_f_Ep_AUC_force()
% choice_Ep_f_Ep_AUC_force aims at assessing whether the choices in the
% physical effort task depend on the area under the curve of the force 
% produced in preceding trials.
%
% INPUTS
% 
% OUTPUTS
%

%% subject identification
[study_nm, ~, ~, subject_id, NS] = sub_id;
%% working directories
% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end
%% main parameters
task_id = 'Ep';
task_fullName = 'physical';
nTrialsPerRun = 54;
n_Em_runs = 2;
nTotalEmTrials = nTrialsPerRun*n_Em_runs;
n_bins = 9;
[choice_f_prevAUC, prevAUC_f_prevAUC,...
    choice_f_sumPrevAUC, sumPrevAUC_f_sumPrevAUC,...
    choice_f_prevPeakF, prevPeakF_f_prevPeakF,...
    choice_f_trialN,...
    AUC_f_trialN, peakF_f_trialN,...
    trialN_f_trialN] = deal(NaN(n_bins, NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];
    %% extract average choice RT across all tasks
    % extract information of runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    nTotalRuns = length(runsStruct.tasks);
    runs = strcmp(runsStruct.tasks, task_id);
    [AUC_avg_run_tmp,...
        prevAUC_avg_run_tmp,...
        sumPrevAUC_avg_run_tmp,...
        peakF_avg_run_tmp,...
        prevPeakF_avg_run_tmp,...
        choice_hE_run_tmp,...
        trialN_run_tmp] = deal(NaN(nTotalEmTrials,1));

    jRun = 0;
    for iRun = 1:nTotalRuns
        runToInclude = 0;
        if runs(iRun) == 1
            jRun = jRun + 1;
            runToInclude = 1;
        end
        run_nm = num2str(iRun);

        if runToInclude == 1
            runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);

            %% load force for the current run
            [~, AUC_tmp, forcePeak_tmp, ~] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);

            %% load data
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                '_task.mat']);
            physicalPerf = behaviorStruct_tmp.physicalPerf;
            % choices
            [choice_hE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);

            for iTrial = 1:nTrialsPerRun
                jTrial = iTrial + nTrialsPerRun*(jRun - 1);
                % choice for current trial
                choice_hE_run_tmp(jTrial) = choice_hE(iTrial);
                % trial number
                trialN_run_tmp(jTrial) = iTrial;
                % performance
                AUC_avg_run_tmp(jTrial) = AUC_tmp.allTrials(iTrial);
                peakF_avg_run_tmp(jTrial) = forcePeak_tmp.allTrials(iTrial);
                % previous trial performance
                if iTrial == 1
                    prevAUC_avg_run_tmp(jTrial) = NaN;
                    prevPeakF_avg_run_tmp(jTrial) = NaN;
                    sumPrevAUC_avg_run_tmp(jTrial) = 0;
                else
                    prevAUC_avg_run_tmp(jTrial) = AUC_avg_run_tmp(jTrial - 1);
                    prevPeakF_avg_run_tmp(jTrial) = peakF_avg_run_tmp(jTrial - 1);
                    sumPrevAUC_avg_run_tmp(jTrial) = sumPrevAUC_avg_run_tmp(jTrial-1) + AUC_avg_run_tmp(jTrial - 1);
                end
            end % trial loop
        end % run to include?
    end % run loop
    %% average data per subject
    [choice_f_prevAUC(:,iS), prevAUC_f_prevAUC(:,iS)] = do_bin2(choice_hE_run_tmp, prevAUC_avg_run_tmp, n_bins,0);
    [choice_f_sumPrevAUC(:,iS), sumPrevAUC_f_sumPrevAUC(:,iS)] = do_bin2(choice_hE_run_tmp, sumPrevAUC_avg_run_tmp, n_bins,0);
    [choice_f_prevPeakF(:,iS), prevPeakF_f_prevPeakF(:,iS)] = do_bin2(choice_hE_run_tmp, prevPeakF_avg_run_tmp, n_bins,0);
    [choice_f_trialN(:,iS), trialN_f_trialN(:,iS)] = do_bin2(choice_hE_run_tmp, trialN_run_tmp, n_bins,0);
    [AUC_f_trialN(:,iS)] = do_bin2(AUC_avg_run_tmp, trialN_run_tmp, n_bins,0);
    [peakF_f_trialN(:,iS)] = do_bin2(peakF_avg_run_tmp, trialN_run_tmp, n_bins,0);
end % subject loop

%% average across subjects
% f force
[choice_f_prevAUC_avg,...
    choice_f_prevAUC_sem] = mean_sem_sd(choice_f_prevAUC,2);
[prevAUC_f_prevAUC_avg,...
    prevAUC_f_prevAUC_sem] = mean_sem_sd(prevAUC_f_prevAUC,2);
[choice_f_sumPrevAUC_avg,...
    choice_f_sumPrevAUC_sem] = mean_sem_sd(choice_f_sumPrevAUC,2);
[sumPrevAUC_f_sumPrevAUC_avg,...
    sumPrevAUC_f_sumPrevAUC_sem] = mean_sem_sd(sumPrevAUC_f_sumPrevAUC,2);
[choice_f_prevPeakF_avg,...
    choice_f_prevPeakF_sem] = mean_sem_sd(choice_f_prevPeakF,2);
[prevPeakF_f_prevPeakF_avg,...
    prevPeakF_f_prevPeakF_sem] = mean_sem_sd(prevPeakF_f_prevPeakF,2);
% f trial number
[choice_f_trialN_avg,...
    choice_f_trialN_sem] = mean_sem_sd(choice_f_trialN,2);
[trialN_f_trialN_avg,...
    trialN_f_trialN_sem] = mean_sem_sd(trialN_f_trialN,2);
[AUC_f_trialN_avg,...
    AUC_f_trialN_sem] = mean_sem_sd(AUC_f_trialN,2);
[peakF_f_trialN_avg,...
    peakF_f_trialN_sem] = mean_sem_sd(peakF_f_trialN,2);

%% figures
lWidth = 3;
pSize = 40;

%% choice = f(previous AUC)
fig;
hold on;
hdl = errorbar(prevAUC_f_prevAUC_avg,...
    choice_f_prevAUC_avg,...
    choice_f_prevAUC_sem);
hdl.LineWidth = lWidth;
xlabel('previous AUC');
ylabel('Choice hE (%)');
legend_size(pSize);

%% choice = f(sum previous AUC)
fig;
hold on;
hdl = errorbar(sumPrevAUC_f_sumPrevAUC_avg,...
    choice_f_sumPrevAUC_avg,...
    choice_f_sumPrevAUC_sem);
hdl.LineWidth = lWidth;
xlabel('sum previous AUC');
ylabel('Choice hE (%)');
legend_size(pSize);

%% choice = f(previous peak force)
fig;
hold on;
hdl = errorbar(prevPeakF_f_prevPeakF_avg,...
    choice_f_prevPeakF_avg,...
    choice_f_prevPeakF_sem);
hdl.LineWidth = lWidth;
xlabel('previous peak force (%)');
ylabel('Choice hE (%)');
legend_size(pSize);

%% AUC = f(trial number)
fig;
hold on;
hdl = errorbar(trialN_f_trialN_avg,...
    AUC_f_trialN_avg,...
    AUC_f_trialN_sem);
hdl.LineWidth = lWidth;
xticks(trialN_f_trialN_avg);
xLabelNames = cell(1,n_bins);
for iTick = 1:n_bins
    xLabelNames{iTick} = num2str(round(trialN_f_trialN_avg(iTick),0));
end
xticklabels(xLabelNames);
xlabel('Trial number');
ylabel('AUC');
legend_size(pSize);

%% peak force = f(trial number)
fig;
hold on;
hdl = errorbar(trialN_f_trialN_avg,...
    peakF_f_trialN_avg,...
    peakF_f_trialN_sem);
hdl.LineWidth = lWidth;
xticks(trialN_f_trialN_avg);
xLabelNames = cell(1,n_bins);
for iTick = 1:n_bins
    xLabelNames{iTick} = num2str(round(trialN_f_trialN_avg(iTick),0));
end
xticklabels(xLabelNames);
xlabel('Trial number');
ylabel('Peak force (%)');
legend_size(pSize);

% end % function