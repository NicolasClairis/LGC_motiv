function[latency, AUC, forcePeak, AUC_overshoot] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm)
% [latency, AUC, forcePeak, AUC_overshoot] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm)
% extract_grip_force will extract several parameters related to performance
% in the phyical effort task, such as latency to squeeze, the area under
% the curve (AUC) of effort performed, the peak of force reached during
% the performance and the area under the curve for the overshoot AUC_overshoot) above the red zone.
%
% INPUTS
% subBehaviorFolder: folder where data is stored
%
% sub_nm: string with subject CID name
%
% run_nm: string with run name (only the number 'X')
%
% OUTPUTS
% latency: estimated time it took for the participant to start squeezing
% the handgrip (careful as this might be very noisy)
%
% AUC: area under the curve expressing the total amount of force that was
% produced
%
% forcePeak: maximal force level for each trial
%
% AUC_overshoot: area under the curve expressing the total amount of
% "useless" force that was produced during the trial (i.e. force above the
% required threshold indicated by the red zone on the screen)
%
% N. Clairis - november 2022

%% general parameters
nTrialsPerRun = 54;
n_hE_levels = 3;
n_Ech = 4;
diffFStartThreshold = 2; % threshold percentage to consider a force is being exerted
Fthreshold = 55;
F_tolerance = 2.5;
F_upper_threshold = Fthreshold + F_tolerance;

[latency.allTrials,...
    AUC.allTrials,...
    forcePeak.allTrials,...
    AUC_overshoot.allTrials] = deal(NaN(1,nTrialsPerRun));
[latency.per_Ech,...
    AUC.per_Ech,...
    forcePeak.per_Ech,...
    AUC_overshoot.per_Ech] = deal(NaN(1,n_Ech));
[latency.per_hE.choice_lowE,...
    AUC.per_hE.choice_lowE,...
    forcePeak.per_hE.choice_lowE,...
    AUC_overshoot.per_hE.choice_lowE,...
    latency.per_hE.choice_highE,...
    AUC.per_hE.choice_highE,...
    forcePeak.per_hE.choice_highE,...
    AUC_overshoot.per_hE.choice_highE] = deal(NaN(1,n_hE_levels));
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
    trialForceLevels_corrected = trialForceLevels - trialForceLevels(1,1);
    
    %% extract peak force
    forcePeak.allTrials(iTrial) = max(trialForceLevels_corrected,[],2,'omitnan');
    %% now get the latency
    latency_squeeze_idx = find(diff(trialForceLevels_corrected) > diffFStartThreshold, 1, 'first');
    if ~isempty(latency_squeeze_idx)
        latency.allTrials(iTrial) = timeForce(latency_squeeze_idx);
    end
    
    %% extract AUC for the whole trial
    AUC.allTrials(iTrial) = trapz(timeForce, trialForceLevels_corrected);

    %% extract AUC overshoot for the whole trial (need to work with actual values, not corrected for this measure)
    timePointsOverShoot = trialForceLevels >= F_upper_threshold;
    drv_clusterIdentif = diff(timePointsOverShoot);
    n_clusters = sum(drv_clusterIdentif == 1);
    if n_clusters == 0
        AUC_overshoot.allTrials(iTrial) = 0;
    elseif n_clusters >= 1
        clusterStarts = find(drv_clusterIdentif == 1);
        clusterEnds = find(drv_clusterIdentif == -1) + 1; % need to add (+1) for the derivative
        % if subject did the task well, last sample should also be a
        % cluster end:
        if (size(clusterStarts,2) > size(clusterEnds,2)) && trialForceLevels(end) >= F_upper_threshold
            clusterEnds = [clusterEnds, length(trialForceLevels)];
        end
        AUC_overshoot.allTrials(iTrial) = 0;
        for iCluster = 1:n_clusters
            if clusterEnds(iCluster) > clusterStarts(iCluster) % no sense to compute AUC if almost no samples
                samples_idx = clusterStarts(iCluster):clusterEnds(iCluster);
                AUC_overshoot.allTrials(iTrial) = AUC_overshoot.allTrials(iTrial) +...
                    trapz(timeForce(samples_idx), trialForceLevels(samples_idx)-F_upper_threshold);
            end
        end % cluster loop
    end
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
    AUC_overshoot.per_hE.choice_lowE(iEff) = mean(AUC_overshoot.allTrials(j_hE_lowEChoice),2,'omitnan');
    % high effort chosen
    latency.per_hE.choice_highE(iEff) = mean(latency.allTrials(j_hE_hEChoice),2,'omitnan');
    AUC.per_hE.choice_highE(iEff) = mean(AUC.allTrials(j_hE_hEChoice),2,'omitnan');
    forcePeak.per_hE.choice_highE(iEff) = mean(forcePeak.allTrials(j_hE_hEChoice),2,'omitnan');
    AUC_overshoot.per_hE.choice_highE(iEff) = mean(AUC_overshoot.allTrials(j_hE_hEChoice),2,'omitnan');
end

% split by effort chosen
E_chosen = behaviorStruct.physicalPerf.E_chosen;
for iEch = 1:n_Ech
    jEch_trials = E_chosen == (iEch - 1);
    latency.per_Ech(iEch) = mean(latency.allTrials(jEch_trials),2,'omitnan');
    AUC.per_Ech(iEch) = mean(AUC.allTrials(jEch_trials),2,'omitnan');
    forcePeak.per_Ech(iEch) = mean(forcePeak.allTrials(jEch_trials),2,'omitnan');
    AUC_overshoot.per_Ech(iEch) = mean(AUC_overshoot.allTrials(jEch_trials),2,'omitnan');
end % effort chosen
end % function