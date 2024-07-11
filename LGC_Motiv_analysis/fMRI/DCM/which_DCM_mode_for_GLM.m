function[DCM_mode] = which_DCM_mode_for_GLM
%[DCM_mode] = which_DCM_mode_for_GLM
% which_DCM_mode_for_GLM asks you how you want to concatenate the onsets
% and regressors in your first level GLM
%
% OUTPUTS
% DCM_mode:
% (1) all sessions modeled independently like in a classic univariate GLM
% => hard to manipulate for DCM but could be useful for testing
% session-specific effects or comparing sessions
% (2) sessions pooled within each task (ex: session 1 and 3 of physical
% effort will be concatenated into one single regressor) but each task will
% be modeled separately
% (3) all sessions pooled together
% (4) all trial periods are pooled together across sessions except for
% choice and effort which are modeled independently for each task (but
% pooled across sessions of the same task)
% (5) all trial periods are pooled together across sessions except for
% the effort period which is modeled independently for each task (but
% pooled across sessions of the same task)

% define list of possible configurations for DCM mode
DCM_modes = {'1: indpdt sessions';... % each session entered independently
    '2: indpdt Ep/Em';... % sessions pooled for each task independently, but Ep and Em modeled separately
    '3: all concatenated';... % pool all sessions together
    '4; all concatenated but choice and effort indpdt for Ep/Em';... % pool all sessions together for all events, except choice + effort
    '5: all concatenated but effort indpdt for Ep/Em'}; % pool all sessions together for all events, except effort

% ask which DCM mode to use?
DCM_mode_idx = listdlg('PromptString','Select which DCM mode you want',...
    'ListString',DCM_modes,'SelectionMode','single');

% extract corresponding number
switch DCM_modes{DCM_mode_idx}
    case '1: indpdt sessions'
        DCM_mode = 1;
    case '2: indpdt Ep/Em'
        DCM_mode = 2;
    case '3: all concatenated'
        DCM_mode = 3;
    case '4; all concatenated but choice and effort indpdt for Ep/Em'
        DCM_mode = 4;
    case '5: all concatenated but effort indpdt for Ep/Em'
        DCM_mode = 5;
end % which DCM mode was selected
end % function