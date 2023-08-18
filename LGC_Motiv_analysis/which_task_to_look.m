function[task_to_look] = which_task_to_look()
% [task_to_look] = which_task_to_look()
% which_task_to_look will ask you which task to look for between physical
% ('Ep'), mental ('Em') or the pool of both ('EpEmPool'). Especially useful
% for fMRI extraction trial by trial with
% extract_ROI_betas_onsets_only.m output and functions using it.
%
% OUTPUTS
% task_to_look: indicating which task to look for
% 'Ep': physical effort task
% 'Em': mental effort task
% 'EpEmPool': physical + mental effort

% define task (physical/mental/both)
task_names = {'Ep','Em','EpEmPool'};
which_task = listdlg('PromptString','Which task to look for in fMRI?','ListString',task_names);
task_to_look = task_names{which_task};

end % function