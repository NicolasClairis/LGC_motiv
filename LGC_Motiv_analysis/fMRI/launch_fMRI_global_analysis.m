
%% define main parameters
if ~exist('GLM','var') || isempty(GLM)
    GLM_str = inputdlg('GLM number?');
    GLM = str2double(GLM_str); % convert from cell to numeric
end
checking = 0;
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
study_nm = 'study1';

%% launch first level, contrasts and 2nd level subsequently
First_level_batch(GLM,checking,condition,study_nm);
LGCM_contrasts_spm(GLM,checking,condition,study_nm);
Second_level_batch(GLM,condition);