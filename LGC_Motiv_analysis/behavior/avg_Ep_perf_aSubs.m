function[latency, AUC, forcePeak, AUC_overshoot,...
    AUC_N, forcePeak_N, AUC_overshoot_N,...
    percTime_above_threshold, percTime_out_of_forceBox] = avg_Ep_perf_aSubs(study_nm, subject_id, condition, NS)
% [latency, AUC, forcePeak, AUC_overshoot,...
%     AUC_N, forcePeak_N, AUC_overshoot_N,...
%     percTime_above_threshold, percTime_out_of_forceBox] = avg_Ep_perf_aSubs(study_nm, subject_id, condition, NS)
% avg_Ep_perf_aSubs
%
% INPUTS
% study_nm: study name (study1 by default)
%
% subject_id: list of subjects
%
% condition: will define which runs should be considered
%
% NS: number of subjects
%
% OUTPUTS
% physical effort (Ep) performance variables in structures splitting and
% averaging the data across chosen effort levels, chosen incentive levels
% and reward/punishment
%
% latency: latency to start squeezing
%
% AUC: area under the curve for force performed
%
% forcePeak: peak force (in % of Fmax)
%
% AUC_overshoot: area under the curve overshoot above the threshold of 
% force to produce to do the trial
%
% AUC_N: area under the curve (based on force in Newtons)
%
% forcePeak_N: peak force in Newtons
%
% AUC_overshoot_N: area under the curve overshoot above the threshold of 
% force to produce to do the trial (based on force in Newtons)
%
% percTime_above_threshold: percentage of time the force was above the 55%
% threshold to reach
%
% percTime_out_of_forceBox: percentage of time the force was below or above
% the red square of required force (ie spending either too much or too low
% force compared to what is necessary)

%% subject selection
% study 1 by default if left empty
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id) 
    [subject_id, NS] = LGCM_subject_selection(study_nm,condition);
end
%% initialize variables of interest
% Effort chosen
Ech_levels = {'Ech0','Ech1','Ech2','Ech3'}; n_Ech = length(Ech_levels);
for iEch = 1:n_Ech
    Ech_nm = Ech_levels{iEch};
    [latency.Ech.(Ech_nm).perSub, AUC.Ech.(Ech_nm).perSub, forcePeak.Ech.(Ech_nm).perSub, AUC_overshoot.Ech.(Ech_nm).perSub,...
        AUC_N.Ech.(Ech_nm).perSub, forcePeak_N.Ech.(Ech_nm).perSub, AUC_overshoot_N.Ech.(Ech_nm).perSub,...
        percTime_above_threshold.Ech.(Ech_nm).perSub, percTime_out_of_forceBox.Ech.(Ech_nm).perSub] = deal(NaN(1,NS));
end % Ech loop

% incentives
Ich_levels = {'inc0','inc1','inc2','inc3'}; n_Ich = length(Ich_levels);
for iInc = 1:n_Ich
    Ich_nm = Ich_levels{iInc};
    [latency.Ich.(Ich_nm).perSub, AUC.Ich.(Ich_nm).perSub, forcePeak.Ich.(Ich_nm).perSub, AUC_overshoot.Ich.(Ich_nm).perSub,...
        AUC_N.Ich.(Ich_nm).perSub, forcePeak_N.Ich.(Ich_nm).perSub, AUC_overshoot_N.Ich.(Ich_nm).perSub,...
        percTime_above_threshold.Ich.(Ich_nm).perSub, percTime_out_of_forceBox.Ich.(Ich_nm).perSub] = deal(NaN(1,NS));
end

% RP
RP_levels = {'R','P'}; n_RP = length(RP_levels);
for iRP = 1:n_RP
    RP_nm = RP_levels{iRP};
    [latency.RP.(RP_nm).perSub, AUC.RP.(RP_nm).perSub, forcePeak.RP.(RP_nm).perSub, AUC_overshoot.RP.(RP_nm).perSub,...
        AUC_N.RP.(RP_nm).perSub, forcePeak_N.RP.(RP_nm).perSub, AUC_overshoot_N.RP.(RP_nm).perSub,...
        percTime_above_threshold.RP.(RP_nm).perSub, percTime_out_of_forceBox.RP.(RP_nm).perSub] = deal(NaN(1,NS));
end

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    % load data for current subject
    [latency_tmp, AUC_tmp, forcePeak_tmp, AUC_overshoot_tmp,...
        AUC_N_tmp, forcePeak_N_tmp, AUC_overshoot_N_tmp,...
        percTime_above_threshold_tmp, percTime_out_of_forceBox_tmp] = avg_Ep_perf_perSub(study_nm, sub_nm, condition);
    
    % extract data per effort level
    for iEch = 1:n_Ech
        Ech_nm = Ech_levels{iEch};
        latency.Ech.(Ech_nm).perSub(iS) = latency_tmp.allRuns.(Ech_nm);
        AUC.Ech.(Ech_nm).perSub(iS) = AUC_tmp.allRuns.(Ech_nm);
        forcePeak.Ech.(Ech_nm).perSub(iS) = forcePeak_tmp.allRuns.(Ech_nm);
        AUC_overshoot.Ech.(Ech_nm).perSub(iS) = AUC_overshoot_tmp.allRuns.(Ech_nm);
        AUC_N.Ech.(Ech_nm).perSub(iS) = AUC_N_tmp.allRuns.(Ech_nm);
        forcePeak_N.Ech.(Ech_nm).perSub(iS) = forcePeak_N_tmp.allRuns.(Ech_nm);
        AUC_overshoot_N.Ech.(Ech_nm).perSub(iS) = AUC_overshoot_N_tmp.allRuns.(Ech_nm);
        percTime_above_threshold.Ech.(Ech_nm).perSub(iS) = percTime_above_threshold_tmp.allRuns.(Ech_nm);
        percTime_out_of_forceBox.Ech.(Ech_nm).perSub(iS) = percTime_out_of_forceBox_tmp.allRuns.(Ech_nm);
    end % effort chosen
    
    % extract data per incentive level
    for iInc = 1:n_Ich
        Ich_nm = Ich_levels{iInc};
        latency.Ich.(Ich_nm).perSub(iS) = latency_tmp.allRuns.(Ich_nm);
        AUC.Ich.(Ich_nm).perSub(iS) = AUC_tmp.allRuns.(Ich_nm);
        forcePeak.Ich.(Ich_nm).perSub(iS) = forcePeak_tmp.allRuns.(Ich_nm);
        AUC_overshoot.Ich.(Ich_nm).perSub(iS) = AUC_overshoot_tmp.allRuns.(Ich_nm);
        AUC_N.Ich.(Ich_nm).perSub(iS) = AUC_N_tmp.allRuns.(Ich_nm);
        forcePeak_N.Ich.(Ich_nm).perSub(iS) = forcePeak_N_tmp.allRuns.(Ich_nm);
        AUC_overshoot_N.Ich.(Ich_nm).perSub(iS) = AUC_overshoot_N_tmp.allRuns.(Ich_nm);
        percTime_above_threshold.Ich.(Ich_nm).perSub(iS) = percTime_above_threshold_tmp.allRuns.(Ich_nm);
        percTime_out_of_forceBox.Ich.(Ich_nm).perSub(iS) = percTime_out_of_forceBox_tmp.allRuns.(Ich_nm);
    end % incentive level chosen
    
    % extract data reward/punishment
    for iRP = 1:n_RP
        RP_nm = RP_levels{iRP};
        latency.RP.(RP_nm).perSub(iS) = latency_tmp.allRuns.(RP_nm);
        AUC.RP.(RP_nm).perSub(iS) = AUC_tmp.allRuns.(RP_nm);
        forcePeak.RP.(RP_nm).perSub(iS) = forcePeak_tmp.allRuns.(RP_nm);
        AUC_overshoot.RP.(RP_nm).perSub(iS) = AUC_overshoot_tmp.allRuns.(RP_nm);
        AUC_N.RP.(RP_nm).perSub(iS) = AUC_N_tmp.allRuns.(RP_nm);
        forcePeak_N.RP.(RP_nm).perSub(iS) = forcePeak_N_tmp.allRuns.(RP_nm);
        AUC_overshoot_N.RP.(RP_nm).perSub(iS) = AUC_overshoot_N_tmp.allRuns.(RP_nm);
        percTime_above_threshold.RP.(RP_nm).perSub(iS) = percTime_above_threshold_tmp.allRuns.(RP_nm);
        percTime_out_of_forceBox.RP.(RP_nm).perSub(iS) = percTime_out_of_forceBox_tmp.allRuns.(RP_nm);
    end % R/P
end % subject loop

%% average across subjects
% effort chosen
for iEch = 1:n_Ech
    Ech_nm = Ech_levels{iEch};
    latency.Ech.(Ech_nm).aSubs = mean(latency.Ech.(Ech_nm).perSub,2,'omitnan');
    AUC.Ech.(Ech_nm).aSubs = mean(AUC.Ech.(Ech_nm).perSub,2,'omitnan');
    forcePeak.Ech.(Ech_nm).aSubs = mean(forcePeak.Ech.(Ech_nm).perSub,2,'omitnan');
    AUC_overshoot.Ech.(Ech_nm).aSubs = mean(AUC_overshoot.Ech.(Ech_nm).perSub,2,'omitnan');
    AUC_N.Ech.(Ech_nm).aSubs = mean(AUC_N.Ech.(Ech_nm).perSub,2,'omitnan');
    forcePeak_N.Ech.(Ech_nm).aSubs = mean(forcePeak_N.Ech.(Ech_nm).perSub,2,'omitnan');
    AUC_overshoot_N.Ech.(Ech_nm).aSubs = mean(AUC_overshoot_N.Ech.(Ech_nm).perSub,2,'omitnan');
    percTime_above_threshold.Ech.(Ech_nm).aSubs = mean(percTime_above_threshold.Ech.(Ech_nm).perSub,2,'omitnan');
    percTime_out_of_forceBox.Ech.(Ech_nm).aSubs = mean(percTime_out_of_forceBox.Ech.(Ech_nm).perSub,2,'omitnan');
end

% incentive chosen
for iInc = 1:n_Ich
    Ich_nm = Ich_levels{iInc};
    latency.Ich.(Ich_nm).aSubs = mean(latency.Ich.(Ich_nm).perSub,2,'omitnan');
    AUC.Ich.(Ich_nm).aSubs = mean(AUC.Ich.(Ich_nm).perSub,2,'omitnan');
    forcePeak.Ich.(Ich_nm).aSubs = mean(forcePeak.Ich.(Ich_nm).perSub,2,'omitnan');
    AUC_overshoot.Ich.(Ich_nm).aSubs = mean(AUC_overshoot.Ich.(Ich_nm).perSub,2,'omitnan');
    AUC_N.Ich.(Ich_nm).aSubs = mean(AUC_N.Ich.(Ich_nm).perSub,2,'omitnan');
    forcePeak_N.Ich.(Ich_nm).aSubs = mean(forcePeak_N.Ich.(Ich_nm).perSub,2,'omitnan');
    AUC_overshoot_N.Ich.(Ich_nm).aSubs = mean(AUC_overshoot_N.Ich.(Ich_nm).perSub,2,'omitnan');
    percTime_above_threshold.Ich.(Ich_nm).aSubs = mean(percTime_above_threshold.Ich.(Ich_nm).perSub,2,'omitnan');
    percTime_out_of_forceBox.Ich.(Ich_nm).aSubs = mean(percTime_out_of_forceBox.Ich.(Ich_nm).perSub,2,'omitnan');
end % incentive level chosen

% R/P
for iRP = 1:n_RP
    RP_nm = RP_levels{iRP};
    latency.RP.(RP_nm).aSubs = mean(latency.RP.(RP_nm).perSub,2,'omitnan');
    AUC.RP.(RP_nm).aSubs = mean(AUC.RP.(RP_nm).perSub,2,'omitnan');
    forcePeak.RP.(RP_nm).aSubs = mean(forcePeak.RP.(RP_nm).perSub,2,'omitnan');
    AUC_overshoot.RP.(RP_nm).aSubs = mean(AUC_overshoot.RP.(RP_nm).perSub,2,'omitnan');
    AUC_N.RP.(RP_nm).aSubs = mean(AUC_N.RP.(RP_nm).perSub,2,'omitnan');
    forcePeak_N.RP.(RP_nm).aSubs = mean(forcePeak_N.RP.(RP_nm).perSub,2,'omitnan');
    AUC_overshoot_N.RP.(RP_nm).aSubs = mean(AUC_overshoot_N.RP.(RP_nm).perSub,2,'omitnan');
    percTime_above_threshold.RP.(RP_nm).aSubs = mean(percTime_above_threshold.RP.(RP_nm).perSub,2,'omitnan');
    percTime_out_of_forceBox.RP.(RP_nm).aSubs = mean(percTime_out_of_forceBox.RP.(RP_nm).perSub,2,'omitnan');
end % R/P

end % function