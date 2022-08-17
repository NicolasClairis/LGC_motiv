function[betas, pval] = RT_GLM(figDisp)
% [betas, pval] = RT_GLM()
%RT_GLM will perform a GLM on the reaction times
% 
% INPUTS
% figDisp: display figure
%
% OUTPUTS
% betas: structure with betas corresponding to each regressor
%
% pval: structure with corresponding p.value for each beta
%
% See also which_RT_GLM

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% study names
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% by default, display group figures
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
    disp(['figDispGroup was not defined in the inputs so that by default ',...
        'group figures are displayed.']);
end
%% working directories
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];
resultFolder_a = [studyBehaviorFolder,'results',filesep];
if ~exist(resultFolder_a,'dir')
    mkdir(resultFolder_a);
end
resultFolder = [resultFolder_a,'figures',filesep];
if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end

%% subject selection
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% initialize variables of interest
GLM_str = inputdlg('Which GLM?');
GLM = str2double(GLM_str);
[GLMprm] = which_RT_GLM(GLM);
potentialRegressors = fieldnames(GLMprm.regs);
for iReg = 1:length(potentialRegressors)
    curr_reg_nm = potentialRegressors{iReg};
    isRegON = strcmp(GLMprm.regs.(curr_reg_nm), 'on');
    if isRegOn == true
        beta.(curr_reg_nm) = NaN(1,NS);
    end
end % regressor loop

%% perform the correlation
for iS = 1:NS
    %% load the data
    
    %% 
    
end % subject loop

%% perform the test

%% display figure
if figDisp == 1
    
end % figure display

end % function