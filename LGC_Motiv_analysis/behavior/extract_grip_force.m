function[latency, AUC, forcePeak] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm)
% [latency, AUC, forcePeak] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm)
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
% latency: estimated time it took for the participant to start squeezing
% the handgrip (careful as this might be very noisy)
%
% AUC: area under the curve expressing the total amount of force that was
% produced
%
% forcePeak: maximal force level for each trial

%% general parameters
nTrialsPerRun = 54;
n_hE_levels = 3;
n_Ech = 4;
diffFStartThreshold = 2; % threshold percentage to consider a force is being exerted

[latency.allTrials,...
    AUC.allTrials,...
    forcePeak.allTrials] = deal(NaN(1,nTrialsPerRun));
[latency.per_Ech,...
    AUC.per_Ech,...
    forcePeak.per_Ech] = deal(NaN(1,n_Ech));
[latency.per_hE.choice_lowE,...
    AUC.per_hE.choice_lowE,...
    forcePeak.per_hE.choice_lowE,...
    latency.per_hE.choice_highE,...
    AUC.per_hE.choice_highE,...
    forcePeak.per_hE.choice_highE] = deal(NaN(1,n_hE_levels));
%% load the data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_physical_task.mat']);

%% extract grip force
for iTrial = 1:nTrialsPerRun
    % extract force
    trialForceLevels = behaviorStruct.physicalPerf.perfSummary{1,iTrial}.displayInformations.forceLevel;
    timeForce = behaviorStruct.physicalPerf.perfSummary{1,iTrial}.displayInformations.time +...
        -behaviorStruct.physicalPerf.onsets.effortPeriod{1,iTrial}.effort_phase;
    % remove baseline which seems to vary weirdly
    trialForceLevels = trialForceLevels - trialForceLevels(1,1);
    
    %% extract peak force
    forcePeak.allTrials(iTrial) = max(trialForceLevels,[],2,'omitnan');
    %% now get the latency
    latency.allTrials(iTrial) = find(diff(trialForceLevels) > diffFStartThreshold, 1, 'first');
    
    %% extract AUC
    AUC.allTrials(iTrial) = trapz(timeForce, trialForceLevels);
end % trial loop

%% split by effort level and effort chosen
% split by high effort level and choice
[hE_level] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, 'physical');
[choice_highE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, 'physical');
choice_lowE = choice_highE == 0;
choice_highE = choice_highE == 1;
for iEff = 1:n_hE_levels
    jEff_level = hE_level == iEff;
    j_hE_lowEChoice = (jEff_level.*choice_lowE) == 1;
    j_hE_hEChoice = (jEff_level.*choice_highE) == 1;
    % low effort chosen
    latency.per_hE.choice_lowE(iEff) = mean(latency.allTrials(j_hE_lowEChoice),2,'omitnan');
    AUC.per_hE.choice_lowE(iEff) = mean(AUC.allTrials(j_hE_lowEChoice),2,'omitnan');
    forcePeak.per_hE.choice_lowE(iEff) = mean(forcePeak.allTrials(j_hE_lowEChoice),2,'omitnan');
    % high effort chosen
    latency.per_hE.choice_highE(iEff) = mean(latency.allTrials(j_hE_hEChoice),2,'omitnan');
    AUC.per_hE.choice_highE(iEff) = mean(AUC.allTrials(j_hE_hEChoice),2,'omitnan');
    forcePeak.per_hE.choice_highE(iEff) = mean(forcePeak.allTrials(j_hE_hEChoice),2,'omitnan');
end

% split by effort chosen
E_chosen = behaviorStruct.physicalPerf.E_chosen;
for iEch = 1:n_Ech
    jEch_trials = E_chosen == (iEch - 1);
    latency.per_Ech(iEch) = mean(latency.allTrials(jEch_trials),2,'omitnan');
    AUC.per_Ech(iEch) = mean(AUC.allTrials(jEch_trials),2,'omitnan');
    forcePeak.per_Ech(iEch) = mean(forcePeak.allTrials(jEch_trials),2,'omitnan');
end % effort chosen
end % function