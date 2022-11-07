function[successSpeed, n_errors, RT_avg] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm)
% [successSpeed, n_errors, RT_avg] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm)
%
%
% INPUTS
% subBehaviorFolder: folder where data is stored
%
% sub_nm: string with subject CID name
%
% run_nm: string with run name
%
% task_fullName: task full name 'mental'/'physical'
%
% OUTPUTS
% successSpeed: time taken to solve all the answers
%
% n_errors: number of errors made
%
% RT_avg: average reaction time for solving the N-back (ignoring the 2
% first presses which are only to initiate the sequence)

%% general parameters
nTrialsPerRun = 54;
n_hE_levels = 3;
n_Ech = 4;

[successSpeed.allTrials,...
    n_errors.allTrials,...
    RT_avg.allTrials,...
    trialSuccess.allTrials] = deal(NaN(1,nTrialsPerRun));
[successSpeed.per_Ech,...
    n_errors.per_Ech,...
    RT_avg.per_Ech] = deal(NaN(1,n_Ech));
[successSpeed.per_hE.choice_lowE,...
    n_errors.per_hE.choice_lowE,...
    RT_avg.per_hE.choice_lowE,...
    successSpeed.per_hE.choice_highE,...
    n_errors.per_hE.choice_highE,...
    RT_avg.per_hE.choice_highE] = deal(NaN(1,n_hE_levels));
%% load the data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_mental_task_behavioral_tmp.mat']);

%% extract grip force
for iTrial = 1:nTrialsPerRun
    % extract whether trial was a success
    trialSuccess.allTrials(iTrial) = behaviorStruct.summary.trial_was_successfull(iTrial);
    if trialSuccess.allTrials(iTrial) == 1 % only extract time when it was a success
        successSpeed.allTrials(iTrial) = behaviorStruct.summary.durations.effortPeriod(iTrial);
    end
    %% extract number of errors made
    n_errors.allTrials(iTrial) = behaviorStruct.summary.perfSummary{1,iTrial}.n_errorsMade;
    %% extract average RT for all answers after the 2 first aimed at initializing the trial
    RT_avg.allTrials(iTrial) = mean(behaviorStruct.summary.perfSummary{1,iTrial}.rt(1,3:end),2,'omitnan');
    
end % trial loop

%% split by effort level and effort chosen
% split by high effort level and choice
[hE_level] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, 'mental');
[choice_highE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, 'mental');
choice_lowE = choice_highE == 0;
choice_highE = choice_highE == 1;
for iEff = 1:n_hE_levels
    jEff_level = hE_level == iEff;
    j_hE_lowEChoice = (jEff_level.*choice_lowE) == 1;
    j_hE_hEChoice = (jEff_level.*choice_highE) == 1;
    % low effort chosen
    successSpeed.per_hE.choice_lowE(iEff) = mean(successSpeed.allTrials(j_hE_lowEChoice),2,'omitnan');
    n_errors.per_hE.choice_lowE(iEff) = mean(n_errors.allTrials(j_hE_lowEChoice),2,'omitnan');
    RT_avg.per_hE.choice_lowE(iEff) = mean(RT_avg.allTrials(j_hE_lowEChoice),2,'omitnan');
    % high effort chosen
    successSpeed.per_hE.choice_highE(iEff) = mean(successSpeed.allTrials(j_hE_hEChoice),2,'omitnan');
    n_errors.per_hE.choice_highE(iEff) = mean(n_errors.allTrials(j_hE_hEChoice),2,'omitnan');
    RT_avg.per_hE.choice_highE(iEff) = mean(RT_avg.allTrials(j_hE_hEChoice),2,'omitnan');
end

% split by effort chosen
E_chosen = behaviorStruct.summary.E_chosen;
for iEch = 1:n_Ech
    jEch_trials = E_chosen == (iEch - 1);
    successSpeed.per_Ech(iEch) = mean(successSpeed.allTrials(jEch_trials),2,'omitnan');
    n_errors.per_Ech(iEch) = mean(n_errors.allTrials(jEch_trials),2,'omitnan');
    RT_avg.per_Ech(iEch) = mean(RT_avg.allTrials(jEch_trials),2,'omitnan');
end % effort chosen
end % function