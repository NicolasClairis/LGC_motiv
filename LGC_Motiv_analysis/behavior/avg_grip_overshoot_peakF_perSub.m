function[AUC_overshoot, peakF, condition] = avg_grip_overshoot_peakF_perSub(subject_id, condition)

%% subject selection
study_nm = 'study1';
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
end

end % function