function[latency, AUC, forcePeak, AUC_overshoot,...
    AUC_N, forcePeak_N, AUC_overshoot_N,...
    perf_duration, percTime_above_threshold, percTime_out_of_forceBox] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm)
% [latency, AUC, forcePeak, AUC_overshoot,...
%   AUC_N, forcePeak_N, AUC_overshoot_N,...
%   perf_duration, percTime_above_threshold, percTime_out_of_forceBox] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm)
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
% produced (in volts, normalized by MVC)
%
% forcePeak: maximal force level for each trial (in volts, normalized by MVC)
%
% AUC_overshoot: area under the curve expressing the total amount of
% "useless" force that was produced during the trial (i.e. force above the
% required threshold indicated by the red zone on the screen) (in volts, normalized by MVC)
%
% AUC: area under the curve expressing the total amount of force that was
% produced (in newtons)
%
% forcePeak: maximal force level for each trial (in newtons)
%
% AUC_overshoot: area under the curve expressing the total amount of
% "useless" force that was produced during the trial (i.e. force above the
% required threshold indicated by the red zone on the screen) (in newtons)
%
% perf_duration: duration between start of squeeze above threshold (ie
% latency) and end of the trial (ie when no more effort needs to be
% exerted)
%
% percTime_above_threshold: percentage of time the force was above the 55%
% threshold to reach
%
% percTime_out_of_forceBox: percentage of time the force was below or above
% the red square of required force (ie spending either too much or too low
% force compared to what is necessary)
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
F_lower_threshold = Fthreshold - F_tolerance;

[latency.allTrials, perf_duration.allTrials] = deal(NaN(1,nTrialsPerRun));
[AUC.allTrials,...
    forcePeak.allTrials,...
    AUC_overshoot.allTrials,...
    AUC_N.allTrials,...
    forcePeak_N.allTrials,...
    AUC_overshoot_N.allTrials,...
    percTime_above_threshold.allTrials,...
    percTime_out_of_forceBox.allTrials] = deal(zeros(1,nTrialsPerRun));
[latency.per_Ech, perf_duration.per_Ech] = deal(NaN(1,n_Ech));
[AUC.per_Ech,...
    forcePeak.per_Ech,...
    AUC_overshoot.per_Ech,...
    AUC_N.per_Ech,...
    forcePeak_N.per_Ech,...
    AUC_overshoot_N.per_Ech,...
    percTime_above_threshold.per_Ech,...
    percTime_out_of_forceBox.per_Ech] = deal(zeros(1,n_Ech));
[latency.per_hE.choice_lowE,...
    latency.per_hE.choice_highE,...
    perf_duration.per_hE.choice_lowE,...
    perf_duration.per_hE.choice_highE] = deal(NaN(1,n_hE_levels));
[AUC.per_hE.choice_lowE,...
    forcePeak.per_hE.choice_lowE,...
    AUC_overshoot.per_hE.choice_lowE,...
    AUC_N.per_hE.choice_lowE,...
    forcePeak_N.per_hE.choice_lowE,...
    AUC_overshoot_N.per_hE.choice_lowE,...
    percTime_above_threshold.per_hE.choice_lowE,...
    percTime_out_of_forceBox.per_hE.choice_lowE,...
    AUC.per_hE.choice_highE,...
    forcePeak.per_hE.choice_highE,...
    AUC_overshoot.per_hE.choice_highE,...
    AUC_N.per_hE.choice_highE,...
    forcePeak_N.per_hE.choice_highE,...
    AUC_overshoot_N.per_hE.choice_highE,...
    percTime_above_threshold.per_hE.choice_highE,...
    percTime_out_of_forceBox.per_hE.choice_highE] = deal(zeros(1,n_hE_levels));
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
    
    % extract force in newtons (instead of (F in volts)/MVC)
    timeForce_N = behaviorStruct.physicalPerf.perfSummary{1,iTrial}.force_levels(:,2) +...
        -behaviorStruct.physicalPerf.onsets.effortPeriod{1,iTrial}.effort_phase;
    % extract voltage
    F_voltage = behaviorStruct.physicalPerf.perfSummary{1,iTrial}.force_levels(:,3);
    nSamples = length(F_voltage);
    % fix situation where no variable recorded
    for iSample = 1:nSamples
        if isnan(F_voltage(iSample))
            if iSample == 1
                F_voltage = 0;
            elseif iSample < nSamples % interpolate between previous and next sample
                F_voltage(iSample) = mean([F_voltage(iSample - 1), F_voltage(iSample+1)],'omitnan');
            elseif iSample == nSamples % no next sample => keep previous value
                F_voltage(iSample) = F_voltage(iSample - 1);
            end
        end
    end
    % convert into newtons
    trialForceLevels_N = grip_biopac_volts_to_newtons_conversion(F_voltage);
    trialForceLevels_N = trialForceLevels_N';
    trialForceLevels_corrected_N = trialForceLevels_N - trialForceLevels_N(1,1);
    MVC_volts = behaviorStruct.physicalPerf.MVC;
    threshold_newtons = grip_biopac_volts_to_newtons_conversion(MVC_volts)*(F_upper_threshold/100);
    
    %% extract peak force
    forcePeak.allTrials(iTrial) = max(trialForceLevels_corrected,[],2,'omitnan');
    forcePeak_N.allTrials(iTrial) = max(trialForceLevels_corrected_N,[],2,'omitnan');
    %% now get the latency
    latency_squeeze_idx = find(diff(trialForceLevels_corrected) > diffFStartThreshold, 1, 'first');
    if ~isempty(latency_squeeze_idx)
        latency.allTrials(iTrial) = timeForce(latency_squeeze_idx);
    end
    
    %% extract performance duration
    if ~isempty(latency_squeeze_idx) % filter trials where some force was exerted
        E_period_dur_tmp = behaviorStruct.physicalPerf.durations.effortPeriod(iTrial);
        perf_duration.allTrials(iTrial) = E_period_dur_tmp - timeForce(latency_squeeze_idx);
    end
    %% extract AUC for the whole trial
    rmv_AUC_NaNs = ~isnan(trialForceLevels_corrected);
    AUC.allTrials(iTrial) = trapz(timeForce(rmv_AUC_NaNs), trialForceLevels_corrected(rmv_AUC_NaNs));
    rmv_AUC_N_NaNs = ~isnan(trialForceLevels_corrected_N);
    AUC_N.allTrials(iTrial) = trapz(timeForce_N(rmv_AUC_N_NaNs), trialForceLevels_corrected_N(rmv_AUC_N_NaNs));
    
    %% extract AUC overshoot for the whole trial (need to work with actual values, not corrected for this measure)
    timePointsOverShoot = trialForceLevels >= F_upper_threshold;
    drv_clusterIdentif = diff(timePointsOverShoot);
    n_clusters = sum(drv_clusterIdentif == 1);
    if n_clusters == 0
        AUC_overshoot.allTrials(iTrial) = 0;
        AUC_overshoot_N.allTrials(iTrial) = 0;
    elseif n_clusters >= 1
        clusterStarts = find(drv_clusterIdentif == 1);
        clusterEnds = find(drv_clusterIdentif == -1) + 1; % need to add (+1) for the derivative
        % if subject did the task well, last sample should also be a
        % cluster end:
        if (size(clusterStarts,2) > size(clusterEnds,2)) && trialForceLevels(end) >= F_upper_threshold
            clusterEnds = [clusterEnds, length(trialForceLevels)];
        end
        AUC_overshoot.allTrials(iTrial) = 0;
        AUC_overshoot_N.allTrials(iTrial) = 0;
        for iCluster = 1:n_clusters
            if clusterEnds(iCluster) > clusterStarts(iCluster) % no sense to compute AUC if almost no samples
                samples_idx = clusterStarts(iCluster):clusterEnds(iCluster);
                % based on volts values
                AUC_overshoot.allTrials(iTrial) = AUC_overshoot.allTrials(iTrial) +...
                    trapz(timeForce(samples_idx), trialForceLevels(samples_idx)-F_upper_threshold);
                % based on newton values
                AUC_overshoot_N.allTrials(iTrial) = AUC_overshoot_N.allTrials(iTrial) +...
                    trapz(timeForce_N(samples_idx), trialForceLevels_N(samples_idx)-threshold_newtons);
            end
        end % cluster loop
    end
    
    %% extract percentage of the effort period in which force was below the threshold
    timeInfo_tmp = behaviorStruct.physicalPerf.perfSummary{1,iTrial}.force_levels(:,2) - behaviorStruct.physicalPerf.perfSummary{1,iTrial}.force_levels(1,2);
    forceLevel_tmp = behaviorStruct.physicalPerf.perfSummary{1,iTrial}.force_levels(:,1);
    trialDur_tmp = timeInfo_tmp(end) - timeInfo_tmp(1);
    jTime_above_threshold = 0;
    jTime_out_of_threshold_box = 0;
    for iT = 2:length(timeInfo_tmp) % ignore the unusual case where subject would start above threshold (should not happen: script prevents this)
        % time force above tolerance threshold
        if forceLevel_tmp(iT) >= F_lower_threshold % force is above tolerance threshold
            jTime_above_threshold = jTime_above_threshold + (timeInfo_tmp(iT) - timeInfo_tmp(iT-1));
        end % lower force threshold
        
        % time force not within tolerance box
        if (forceLevel_tmp(iT) < F_lower_threshold) || (forceLevel_tmp(iT) > F_upper_threshold)
            jTime_out_of_threshold_box = jTime_out_of_threshold_box + (timeInfo_tmp(iT) - timeInfo_tmp(iT-1));
        end % out of red box
    end % time loop
    percTime_above_threshold.allTrials(iTrial) = jTime_above_threshold/trialDur_tmp;
    percTime_out_of_forceBox.allTrials(iTrial) = jTime_out_of_threshold_box/trialDur_tmp;
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
    perf_duration.per_hE.choice_lowE(iEff) = mean(perf_duration.allTrials(j_hE_lowEChoice),2,'omitnan');
    AUC.per_hE.choice_lowE(iEff) = mean(AUC.allTrials(j_hE_lowEChoice),2,'omitnan');
    forcePeak.per_hE.choice_lowE(iEff) = mean(forcePeak.allTrials(j_hE_lowEChoice),2,'omitnan');
    AUC_overshoot.per_hE.choice_lowE(iEff) = mean(AUC_overshoot.allTrials(j_hE_lowEChoice),2,'omitnan');
    AUC_N.per_hE.choice_lowE(iEff) = mean(AUC_N.allTrials(j_hE_lowEChoice),2,'omitnan');
    forcePeak_N.per_hE.choice_lowE(iEff) = mean(forcePeak_N.allTrials(j_hE_lowEChoice),2,'omitnan');
    AUC_overshoot_N.per_hE.choice_lowE(iEff) = mean(AUC_overshoot_N.allTrials(j_hE_lowEChoice),2,'omitnan');
    percTime_above_threshold.per_hE.choice_lowE(iEff) = mean(percTime_above_threshold.allTrials(j_hE_lowEChoice),2,'omitnan');
    percTime_out_of_forceBox.per_hE.choice_lowE(iEff) = mean(percTime_out_of_forceBox.allTrials(j_hE_lowEChoice),2,'omitnan');
    % high effort chosen
    latency.per_hE.choice_highE(iEff) = mean(latency.allTrials(j_hE_hEChoice),2,'omitnan');
    perf_duration.per_hE.choice_highE(iEff) = mean(perf_duration.allTrials(j_hE_hEChoice),2,'omitnan');
    AUC.per_hE.choice_highE(iEff) = mean(AUC.allTrials(j_hE_hEChoice),2,'omitnan');
    forcePeak.per_hE.choice_highE(iEff) = mean(forcePeak.allTrials(j_hE_hEChoice),2,'omitnan');
    AUC_overshoot.per_hE.choice_highE(iEff) = mean(AUC_overshoot.allTrials(j_hE_hEChoice),2,'omitnan');
    AUC_N.per_hE.choice_highE(iEff) = mean(AUC_N.allTrials(j_hE_hEChoice),2,'omitnan');
    forcePeak_N.per_hE.choice_highE(iEff) = mean(forcePeak_N.allTrials(j_hE_hEChoice),2,'omitnan');
    AUC_overshoot_N.per_hE.choice_highE(iEff) = mean(AUC_overshoot_N.allTrials(j_hE_hEChoice),2,'omitnan');
    percTime_above_threshold.per_hE.choice_highE(iEff) = mean(percTime_above_threshold.allTrials(j_hE_hEChoice),2,'omitnan');
    percTime_out_of_forceBox.per_hE.choice_highE(iEff) = mean(percTime_out_of_forceBox.allTrials(j_hE_hEChoice),2,'omitnan');
end

% split by effort chosen
E_chosen = behaviorStruct.physicalPerf.E_chosen;
for iEch = 1:n_Ech
    jEch_trials = E_chosen == (iEch - 1);
    latency.per_Ech(iEch) = mean(latency.allTrials(jEch_trials),2,'omitnan');
    perf_duration.per_Ech(iEch) = mean(perf_duration.allTrials(jEch_trials),2,'omitnan');
    AUC.per_Ech(iEch) = mean(AUC.allTrials(jEch_trials),2,'omitnan');
    forcePeak.per_Ech(iEch) = mean(forcePeak.allTrials(jEch_trials),2,'omitnan');
    AUC_overshoot.per_Ech(iEch) = mean(AUC_overshoot.allTrials(jEch_trials),2,'omitnan');
    AUC_N.per_Ech(iEch) = mean(AUC_N.allTrials(jEch_trials),2,'omitnan');
    forcePeak_N.per_Ech(iEch) = mean(forcePeak_N.allTrials(jEch_trials),2,'omitnan');
    AUC_overshoot_N.per_Ech(iEch) = mean(AUC_overshoot_N.allTrials(jEch_trials),2,'omitnan');
    percTime_above_threshold.per_Ech(iEch) = mean(percTime_above_threshold.allTrials(jEch_trials),2,'omitnan');
    percTime_out_of_forceBox.per_Ech(iEch) = mean(percTime_out_of_forceBox.allTrials(jEch_trials),2,'omitnan');
end % effort chosen
end % function