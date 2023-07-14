
%% define main parameters
GLM = inputdlg('GLM number?');
checking = 0;
condition = subject_condition;
study_nm = 'study1';

%% launch first level, contrasts and 2nd level subsequently
First_level_batch(GLM,checking,condition,study_nm);
LGCM_contrasts_spm(GLM,checking,condition,study_nm);
Second_level_batch(GLM,condition);