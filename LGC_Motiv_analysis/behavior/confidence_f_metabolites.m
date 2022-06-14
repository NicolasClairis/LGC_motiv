%% script to check whether there is a global difference in confidence 
% ratings and in inferred confidence based on the level of brain
% metabolites.


%% working directory
[computerRoot] = LGCM_root_paths();
study_nm = 'study1';
studyFolder = fullfile(computerRoot, study_nm);

%% define all subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection('study1', condition);

%% define metabolite and ROI you want to focus on and extract subjects accordingly
[low_met_subs, high_met_subs, metabolite_nm] = medSplit_metabolites(subject_id);