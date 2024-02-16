function[ROI_trial_b_trial, ROI_subList,...
    ROI_nm, ROI_short_nm,...
    timePeriod_nm] = extract_ROI_betas_onsets_only_bis(computerRoot,...
    study_nm, subject_id, condition)
% [ROI_trial_b_trial, ROI_subList,...
%     ROI_nm, ROI_short_nm,...
%     timePeriod_nm] = extract_ROI_betas_onsets_only_bis(computerRoot,...
%     study_nm, subject_id, condition)
% extract_ROI_betas_onsets_only_bis calls extract_ROI_betas_onsets_only
% to extract data trial by trial and then will ask you more specifically
% what you want to focus on within this data to facilitate further use by
% other scripts.
%
% INPUTS
% computerRoot: path to computer root
%
% study_nm: string with study name
%
% subject_id: list of subjects to extract
%
% condition: which subjects and runs should be considered
%
% OUTPUTS
% ROI_trial_b_trial: ROI extracted trial by trial for each subject and each
% task and time period
%
% ROI_subList: list of the subjects that were extracted
%
% ROI_nm: ROI where data was extracted full name
%
% ROI_short_nm: short name asked to facilitate figures and stuff
%
% timePeriod_nm: time period where to look

%% extract the data (slower part)
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);

%% define which ROI, and which time period is of interest to you
%% define ROI of interest
ROI_names = fieldnames(ROI_trial_b_trial);
ROI_subList = ROI_trial_b_trial.subject_id;
ROI_names(strcmp(ROI_names,'subject_id')) = [];
if length(ROI_names) > 1
    which_ROI = listdlg('PromptString','Which ROI?','ListString',ROI_names);
    ROI_nm = ROI_names(which_ROI);
else
    ROI_nm = ROI_names;
end
ROI_short_nm = inputdlg('ROI short name?');
ROI_short_nm = ROI_short_nm{1};
% in case spaces were entered, replace them by '_'
ROI_short_nm = strrep(ROI_short_nm,' ','_');

%% define time period (depends on GLM selected)
timePeriods = fieldnames(ROI_trial_b_trial.(ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};
end