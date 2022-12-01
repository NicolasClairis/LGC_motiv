function[successSpeed, n_correct, n_errors, RT_avg,...
    RTtotal_with2firstUseless,...
    RTtotal_pureNback,...
    efficacy_with2first,...
    efficacy_pureNback] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm)
% [successSpeed, n_correct, n_errors, RT_avg,...
%     RTtotal_with2firstUseless,...
%     RTtotal_pureNback,...
%     efficacy_with2first,...
%     efficacy_pureNback] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm)
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
% n_correct: number of correct answers provided
%
% n_errors: number of errors made
%
% RT_avg: average reaction time for solving the N-back (ignoring the 2
% first presses which are only to initiate the sequence)
%
% RTtotal_with2firstUseless: time spent for effort exertion including the
% time spent solving the two first useless answers
%
% RTtotal_pureNback: time spent for effort exertion ignoring the time spent
% at the beginning for the two first useless answers
%
% efficacy_with2first: efficacy expressed as
% (n_correct-n_errors)/RT_total_with2firstUseless
%
% efficacy_pureNback: efficacy expressed as
% (n_correct-n_errors)/RTtotal_pureNback
%
% N. Clairis - november 2022

%% general parameters
nTrialsPerRun = 54;
n_hE_levels = 3;
n_Ech = 4;

[successSpeed.allTrials,...
    n_correct.allTrials,...
    n_errors.allTrials,...
    RT_avg.allTrials,...
    trialSuccess.allTrials,...
    RTtotal_with2firstUseless.allTrials,...
    RTtotal_pureNback.allTrials,...
    efficacy_with2first.allTrials,...
    efficacy_pureNback.allTrials] = deal(NaN(1,nTrialsPerRun));
[successSpeed.per_Ech,...
    n_correct.per_Ech,...
    n_errors.per_Ech,...
    RT_avg.per_Ech,...
    RTtotal_with2firstUseless.per_Ech,...
    RTtotal_pureNback.per_Ech,...
    efficacy_with2first.per_Ech,...
    efficacy_pureNback.per_Ech] = deal(NaN(1,n_Ech));
[successSpeed.per_hE.choice_lowE,...
    n_correct.per_hE.choice_lowE,...
    n_errors.per_hE.choice_lowE,...
    RT_avg.per_hE.choice_lowE,...,...
    RTtotal_with2firstUseless.per_hE.choice_lowE,...
    RTtotal_pureNback.per_hE.choice_lowE,...
    efficacy_with2first.per_hE.choice_lowE,...
    efficacy_pureNback.per_hE.choice_lowE,...
    successSpeed.per_hE.choice_highE,...
    n_correct.per_hE.choice_highE,...
    n_errors.per_hE.choice_highE,...
    RT_avg.per_hE.choice_highE,...
    RTtotal_with2firstUseless.per_hE.choice_highE,...
    RTtotal_pureNback.per_hE.choice_highE,...
    efficacy_with2first.per_hE.choice_highE,...
    efficacy_pureNback.per_hE.choice_highE] = deal(NaN(1,n_hE_levels));
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
    %% extract number of correct answers
    n_correct.allTrials(iTrial) = behaviorStruct.summary.perfSummary{1,iTrial}.n_correctAnswersProvided;
    %% extract number of errors made
    n_errors.allTrials(iTrial) = behaviorStruct.summary.perfSummary{1,iTrial}.n_errorsMade;
    %% extract average RT for all answers after the 2 first aimed at initializing the trial
    RT_avg.allTrials(iTrial) = mean(behaviorStruct.summary.perfSummary{1,iTrial}.rt(1,3:end),2,'omitnan');
    %% extract time spent to do the effort
    onset_names_tmp = fieldnames(behaviorStruct.summary.perfSummary{1,iTrial}.onsets);
    if ~isempty(onset_names_tmp) && sum(strcmp(onset_names_tmp,'nb_3')) == 1
        trialStart_tmp = behaviorStruct.summary.perfSummary{1,iTrial}.onsets.nb_1;
        trialNbackStart_tmp = behaviorStruct.summary.perfSummary{1,iTrial}.onsets.nb_3; % ignore 2 first useless answers
        trialEnd_tmp = behaviorStruct.summary.perfSummary{1,iTrial}.onsets.(['nb_',num2str(length(onset_names_tmp))]) +...
            behaviorStruct.summary.perfSummary{1,iTrial}.rt(end);
        RTtotal_with2firstUseless.allTrials(iTrial) = trialEnd_tmp - trialStart_tmp;
        RTtotal_pureNback.allTrials(iTrial) = trialEnd_tmp - trialNbackStart_tmp;
    end
    %% extract efficacy
    efficacy_with2first.allTrials(iTrial) = (n_correct.allTrials(iTrial) - n_errors.allTrials(iTrial))./RTtotal_with2firstUseless.allTrials(iTrial);
    efficacy_pureNback.allTrials(iTrial) = (n_correct.allTrials(iTrial) - n_errors.allTrials(iTrial))./RTtotal_pureNback.allTrials(iTrial);
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
    n_correct.per_hE.choice_lowE(iEff) = mean(n_correct.allTrials(j_hE_lowEChoice),2,'omitnan');
    n_errors.per_hE.choice_lowE(iEff) = mean(n_errors.allTrials(j_hE_lowEChoice),2,'omitnan');
    RT_avg.per_hE.choice_lowE(iEff) = mean(RT_avg.allTrials(j_hE_lowEChoice),2,'omitnan');
    RTtotal_with2firstUseless.per_hE.choice_lowE(iEff) = mean(RTtotal_with2firstUseless.allTrials(j_hE_lowEChoice),2,'omitnan');
    RTtotal_pureNback.per_hE.choice_lowE(iEff) = mean(RTtotal_pureNback.allTrials(j_hE_lowEChoice),2,'omitnan');
    efficacy_with2first.per_hE.choice_lowE(iEff) = mean(efficacy_with2first.allTrials(j_hE_lowEChoice),2,'omitnan');
    efficacy_pureNback.per_hE.choice_lowE(iEff) = mean(efficacy_pureNback.allTrials(j_hE_lowEChoice),2,'omitnan');
    % high effort chosen
    successSpeed.per_hE.choice_highE(iEff) = mean(successSpeed.allTrials(j_hE_hEChoice),2,'omitnan');
    n_correct.per_hE.choice_highE(iEff) = mean(n_correct.allTrials(j_hE_hEChoice),2,'omitnan');
    n_errors.per_hE.choice_highE(iEff) = mean(n_errors.allTrials(j_hE_hEChoice),2,'omitnan');
    RT_avg.per_hE.choice_highE(iEff) = mean(RT_avg.allTrials(j_hE_hEChoice),2,'omitnan');
    RTtotal_with2firstUseless.per_hE.choice_highE(iEff) = mean(RTtotal_with2firstUseless.allTrials(j_hE_hEChoice),2,'omitnan');
    RTtotal_pureNback.per_hE.choice_highE(iEff) = mean(RTtotal_pureNback.allTrials(j_hE_hEChoice),2,'omitnan');
    efficacy_with2first.per_hE.choice_highE(iEff) = mean(efficacy_with2first.allTrials(j_hE_hEChoice),2,'omitnan');
    efficacy_pureNback.per_hE.choice_highE(iEff) = mean(efficacy_pureNback.allTrials(j_hE_hEChoice),2,'omitnan');
end

% split by effort chosen
E_chosen = behaviorStruct.summary.E_chosen;
for iEch = 1:n_Ech
    jEch_trials = E_chosen == (iEch - 1);
    successSpeed.per_Ech(iEch) = mean(successSpeed.allTrials(jEch_trials),2,'omitnan');
    n_correct.per_Ech(iEch) = mean(n_correct.allTrials(jEch_trials),2,'omitnan');
    n_errors.per_Ech(iEch) = mean(n_errors.allTrials(jEch_trials),2,'omitnan');
    RT_avg.per_Ech(iEch) = mean(RT_avg.allTrials(jEch_trials),2,'omitnan');
    RTtotal_with2firstUseless.per_Ech(iEch) = mean(RTtotal_with2firstUseless.allTrials(jEch_trials),2,'omitnan');
    RTtotal_pureNback.per_Ech(iEch) = mean(RTtotal_pureNback.allTrials(jEch_trials),2,'omitnan');
    efficacy_with2first.per_Ech(iEch) = mean(efficacy_with2first.allTrials(jEch_trials),2,'omitnan');
    efficacy_pureNback.per_Ech(iEch) = mean(efficacy_pureNback.allTrials(jEch_trials),2,'omitnan');
end % effort chosen
end % function