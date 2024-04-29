function[forcePeak, percTime_above_threshold, percTime_out_of_forceBox] = extract_Ep_force_accuracy(study_nm, subject_id, condition)
% [forcePeak, percTime_above_threshold, percTime_out_of_forceBox] = extract_Ep_force_accuracy(study_nm, subject_id, condition)
% extract_Ep_force_accuracy will extract variables related to force
% accuracy in the physical effort task
%
% INPUTS
% study_nm: study name
%
% subject_id: cell with list of subjects
%
% condition: condition for filtering subjects
%
% OUTPUTS
% forcePeak: average peak force across all trials in % of MVC
%
% percTime_above_threshold: average percentage of effort period time where
% force was above the required threshold (gives an indication of the
% average time spent doing the required force vs not doing enough force)
%
% percTime_out_of_forceBox: average percentage of effort period time where
% force was either below the required threshold or above the required
% threshold (ie too low force exerted = yellow circle does not reduce; or
% too high force exerted = yellow circle reduces but too much energy is
% being spent)

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
% condition selection
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
% subject selection
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm,condition);
else
    NS = length(subject_id);
end

%% initialize variables of interest
[forcePeak, percTime_above_threshold, percTime_out_of_forceBox] = deal(NaN(1,NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    [~, ~, forcePeak_tmp, ~,...
        ~, ~, ~,...
        percTime_above_threshold_tmp, percTime_out_of_forceBox_tmp] = avg_Ep_perf_perSub(study_nm, sub_nm, condition);
    forcePeak(iS) = forcePeak_tmp.allRuns.allTrials;
    percTime_above_threshold(iS) = percTime_above_threshold_tmp.allRuns.allTrials;
    percTime_out_of_forceBox(iS) = percTime_out_of_forceBox_tmp.allRuns.allTrials;
end % subject loop

end % function