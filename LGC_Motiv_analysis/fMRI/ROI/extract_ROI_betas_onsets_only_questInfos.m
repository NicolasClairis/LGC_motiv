function[fMRI_ROI_nm, fMRI_ROI_short_nm, ROI_task_to_look, timePeriod_nm] = extract_ROI_betas_onsets_only_questInfos(ROI_trial_b_trial)
% [fMRI_ROI_nm, ROI_short_nm, ROI_task_to_look, timePeriod_nm] = extract_ROI_betas_onsets_only_questInfos(ROI_trial_b_trial)
% extract_ROI_betas_onsets_only_questInfos will get ROI_trial_b_trial
% structure as an input (extracted with extract_ROI_betas_onsets_only.m)
% and it will then ask you which ROI you want to focus on, which name you
% want to give it, which task you want to look and also which time period
% you want to look for.
%
% INPUTS
% ROI_trial_b_trial: structure with information of ROI trial by trial
%
% OUTPUTS
% fMRI_ROI_nm: full name of the ROI selected
%
% fMRI_ROI_short_nm: short name of the ROI selected, entered by the user in
% order to give it a more friendly and short name
%
% ROI_task_to_look: focus on physical (Ep), mental (Em) or both (EpEmPool) tasks?
%
% timePeriod_nm: select which time period of the trial you want to focus on
% based on what was modeled in the GLM selected.
%
% See also extract_ROI_betas_onsets_only.m

% define which ROI, and which time period is of interest to you
% define ROI

%% define ROI name
ROI_names = fieldnames(ROI_trial_b_trial);
ROI_names(strcmp(ROI_names,'subject_id')) = [];
if length(ROI_names) > 1
    error(['There should be only 1 ROI selected, not ',num2str(length(ROI_names))])
else
    fMRI_ROI_nm = ROI_names;
end
fMRI_ROI_short_nm = inputdlg('ROI short name?');
fMRI_ROI_short_nm = fMRI_ROI_short_nm{1};

%% define task
task_names = {'Ep','Em','EpEmPool'};
which_ROI_task = listdlg('PromptString','Which task for ROI?','ListString',task_names);
ROI_task_to_look = task_names{which_ROI_task};

%% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(fMRI_ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};

end % function