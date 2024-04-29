
%% define main parameters
% which GLM to use?
if ~exist('GLM','var') || isempty(GLM)
    GLM_str = inputdlg('GLM number?');
    GLM = str2double(GLM_str); % convert from cell to numeric
end

% no checking (launch it all)
checking = 0;

% condition
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end

% study name
study_nm = 'study1';

% use bias-field corrected  images or raw images?
biasFieldCorr = 0;

%% launch first level, contrasts and 2nd level subsequently
First_level_batch(GLM,checking,condition,study_nm, [],[], biasFieldCorr);
LGCM_contrasts_spm(GLM,checking,condition,study_nm, [],[], biasFieldCorr);
Second_level_batch(GLM,condition, [],biasFieldCorr);