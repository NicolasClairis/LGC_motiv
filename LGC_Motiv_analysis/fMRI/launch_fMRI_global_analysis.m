
%% define main parameters
GLM_str = inputdlg('GLM number?');
GLM = str2double(GLM_str); % convert from cell to numeric
checking = 0;
condition = subject_condition;
study_nm = 'study1';

%% launch first level, contrasts and 2nd level subsequently
First_level_batch(GLM,checking,condition,study_nm);
LGCM_contrasts_spm(GLM,checking,condition,study_nm);
Second_level_batch(GLM,condition);