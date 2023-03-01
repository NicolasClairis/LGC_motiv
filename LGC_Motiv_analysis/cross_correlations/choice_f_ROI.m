function[] = choice_f_ROI()
% choice_f_ROI will look at percentage of choices depending on ROI level of
% activity. Then it will also split the data depending on the effort level
% proposed and try to look whether the results vary accordingly.


%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% subject selection
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% general parameters
if ~exist('nBins','var') || isempty(nBins)
    nBins = 9;
end

%% extract ROI activity for all subjects
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);
% define which ROI, and which time period is of interest to you
% define ROI
[fMRI_ROI_nm, ROI_short_nm,...
    ROI_task_to_look,...
    timePeriod_nm] = extract_ROI_betas_onsets_only_questInfos(ROI_trial_b_trial);

%% loop through subjects
for iS = 1:NS
    
end % subject loop

%% average 

%% figure

end % function