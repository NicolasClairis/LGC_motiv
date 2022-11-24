% function[] = choice_Ep_f_Ep_Trest()
% choice_Ep_f_Ep_Trest aims at assessing whether the choices in the
% physical effort task depend on the time resting in preceding trials.
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
n_bins = 6;
[choice_f_Trest1, Trest1_f_Trest1,...
    choice_f_Trest2, Trest2_f_Trest2,...
    choice_f_trialN,...
    Trest1_f_trialN, Trest2_f_trialN,...
    trialN_f_trialN] = deal(NaN(n_bins, NS));
n_sumPrevAUC_bins = 6;
for iAUCbin = 1:n_sumPrevAUC_bins
    bin_nm = ['bin',num2str(iAUCbin)];
    [choice_f_Trest1_AUCbin.(bin_nm),...
        Trest1_f_Trest1_AUCbin.(bin_nm)] = deal(NaN(n_bins, NS));
end % loop previous AUC bins

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
    [Trest1_avg_run_tmp,...
        Trest2_avg_run_tmp,...
        sumPrevAUC_avg_run_tmp,...
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

            %% load time resting for the current run
            [Trest1, Trest2] = extract_grip_Trest(subBehaviorFolder, sub_nm, run_nm);
            %% load force for the current run
            [~, AUC_tmp] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
            for iTrial = 1:nTrialsPerRun
                jTrial = iTrial + nTrialsPerRun*(jRun - 1);
                if iTrial == 1
                    sumPrevAUC_avg_run_tmp(jTrial) = 0;
                else
                    sumPrevAUC_avg_run_tmp(jTrial) = sumPrevAUC_avg_run_tmp(jTrial-1) +...
                        AUC_tmp.allTrials(iTrial-1);
                end
            end % trial loop

            %% load data
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                '_task.mat']);
            physicalPerf = behaviorStruct_tmp.physicalPerf;
            % choices
            [choice_hE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);

            % choice for current trial
            choice_hE_run_tmp(runTrials_idx) = choice_hE;
            % trial number
            trialN_run_tmp(runTrials_idx) = 1:nTrialsPerRun;
            % performance
            Trest1_avg_run_tmp(runTrials_idx) = Trest1;
            Trest2_avg_run_tmp(runTrials_idx) = Trest2;
        end % run to include?
    end % run loop
    %% average data per subject
    [choice_f_Trest1(:,iS), Trest1_f_Trest1(:,iS)] = do_bin2(choice_hE_run_tmp, Trest1_avg_run_tmp, n_bins,0);
    [choice_f_Trest2(:,iS), Trest2_f_Trest2(:,iS)] = do_bin2(choice_hE_run_tmp, Trest2_avg_run_tmp, n_bins,0);
    [choice_f_trialN(:,iS), trialN_f_trialN(:,iS)] = do_bin2(choice_hE_run_tmp, trialN_run_tmp, n_bins,0);
    [Trest1_f_trialN(:,iS)] = do_bin2(Trest1_avg_run_tmp, trialN_run_tmp, n_bins,0);
    [Trest2_f_trialN(:,iS)] = do_bin2(Trest2_avg_run_tmp, trialN_run_tmp, n_bins,0);


    %% split data according to level of fatigue
    [~,~,AUC_bin_idx] = do_bin2(choice_hE_run_tmp,sumPrevAUC_avg_run_tmp, n_sumPrevAUC_bins,0);
    for iBin = 1:n_sumPrevAUC_bins
        bin_nm = ['bin',num2str(iBin)];
        bin_idx = AUC_bin_idx == iBin;
        [choice_f_Trest1_AUCbin.(bin_nm)(:,iS),...
            Trest1_f_Trest1_AUCbin.(bin_nm)(:,iS)] = do_bin2(choice_hE_run_tmp(bin_idx), Trest1_avg_run_tmp(bin_idx), n_bins,0);
    end % bin loop
end % subject loop

%% average across subjects
% f time resting
[choice_f_Trest1_avg,...
    choice_f_Trest1_sem] = mean_sem_sd(choice_f_Trest1,2);
[Trest1_f_Trest1_avg,...
    Trest1_f_Trest1_sem] = mean_sem_sd(Trest1_f_Trest1,2);
[choice_f_Trest2_avg,...
    choice_f_Trest2_sem] = mean_sem_sd(choice_f_Trest2,2);
[Trest2_f_Trest2_avg,...
    Trest2_f_Trest2_sem] = mean_sem_sd(Trest2_f_Trest2,2);
% f trial number
[choice_f_trialN_avg,...
    choice_f_trialN_sem] = mean_sem_sd(choice_f_trialN,2);
[trialN_f_trialN_avg,...
    trialN_f_trialN_sem] = mean_sem_sd(trialN_f_trialN,2);
[Trest1_f_trialN_avg,...
    Trest1_f_trialN_sem] = mean_sem_sd(Trest1_f_trialN,2);
[Trest2_f_trialN_avg,...
    Trest2_f_trialN_sem] = mean_sem_sd(Trest2_f_trialN,2);

% f time resting split by fatigue level (sum AUC)
for iBin = 1:n_sumPrevAUC_bins
    bin_nm = ['bin',num2str(iBin)];
    [choice_f_Trest1_AUCbin_avg.(bin_nm),...
        choice_f_Trest1_AUCbin_sem.(bin_nm)] = mean_sem_sd(choice_f_Trest1_AUCbin.(bin_nm),2);
    [Trest1_f_Trest1_AUCbin_avg.(bin_nm),...
        Trest1_f_Trest1_AUCbin_sem.(bin_nm)] = mean_sem_sd(Trest1_f_Trest1_AUCbin.(bin_nm),2);
end

%% figures
lWidth = 3;
pSize = 40;

%% choice = f(Trest1)
fig;
hold on;
hdl = errorbar(Trest1_f_Trest1_avg,...
    choice_f_Trest1_avg,...
    choice_f_Trest1_sem);
hdl.LineWidth = lWidth;
xlabel('Time resting (s)');
ylabel('Choice hE (%)');
legend_size(pSize);

%% choice = f(Trest2)
fig;
hold on;
hdl = errorbar(Trest2_f_Trest2_avg,...
    choice_f_Trest2_avg,...
    choice_f_Trest2_sem);
hdl.LineWidth = lWidth;
xlabel('Time resting (s) ignoring easy efforts');
ylabel('Choice hE (%)');
legend_size(pSize);

%% Trest1 = f(trial number)
fig;
hold on;
hdl = errorbar(trialN_f_trialN_avg,...
    Trest1_f_trialN_avg,...
    Trest1_f_trialN_sem);
hdl.LineWidth = lWidth;
xticks(trialN_f_trialN_avg);
xLabelNames = cell(1,n_bins);
for iTick = 1:n_bins
    xLabelNames{iTick} = num2str(round(trialN_f_trialN_avg(iTick),0));
end
xticklabels(xLabelNames);
xlabel('Trial number');
ylabel('Time resting (s)');
legend_size(pSize);

%% Trest2 = f(trial number)
fig;
hold on;
hdl = errorbar(trialN_f_trialN_avg,...
    Trest2_f_trialN_avg,...
    Trest2_f_trialN_sem);
hdl.LineWidth = lWidth;
xticks(trialN_f_trialN_avg);
xLabelNames = cell(1,n_bins);
for iTick = 1:n_bins
    xLabelNames{iTick} = num2str(round(trialN_f_trialN_avg(iTick),0));
end
xticklabels(xLabelNames);
xlabel('Trial number');
ylabel('Time resting (s) ignoring easy efforts');
legend_size(pSize);

%% choice = f(Trest 1) depending on level of fatigue
fig;
hold on;
legend_hdl = NaN(1,n_sumPrevAUC_bins);
legend_nm_hdl = cell(1,n_sumPrevAUC_bins);
for iBin = 1:n_sumPrevAUC_bins
    bin_nm = ['bin',num2str(iBin)];
    legend_nm_hdl{iBin} = ['Fatigue ',num2str(iBin)];
    legend_hdl_tmp = errorbar(Trest1_f_Trest1_AUCbin_avg.(bin_nm),...
        choice_f_Trest1_AUCbin_avg.(bin_nm),...
        choice_f_Trest1_AUCbin_sem.(bin_nm));
    legend_hdl_tmp.LineWidth = lWidth;
    legend_hdl_tmp.Color = [1-iBin/(2*n_sumPrevAUC_bins),...
        1-iBin/(2*n_sumPrevAUC_bins),...
        iBin/n_sumPrevAUC_bins];
    legend_hdl(iBin) = legend_hdl_tmp;
end
legend(legend_hdl, legend_nm_hdl);
legend('boxoff');
xlabel('Time resting (s)');
ylabel('Choice hE (%)');
legend_size(pSize);

% end % function