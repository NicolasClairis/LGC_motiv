function [conNamesPerTaskRun, con_simpleVectors, n_conPerTaskRun] = list_contrasts_perTaskRun(GLM, task_nm)
%[conNamesPerTaskRun, con_simpleVectors, n_conPerTaskRun] = list_contrasts_perTaskRun(GLM, task_nm)
%
% INPUTS
% GLM: GLM number
%
% task_nm: task name ('physical'/'mental')
%
% OUTPUTS
% conNamesPerTaskRun: cell with name of each contrast
% 
% con_simpleVectors: matrix with contrasts of interest
%
% n_conPerTaskRun: number of contrasts per run

%% extract GLM parameters
GLMprm = which_GLM(GLM);

%% task id
switch task_nm
    case 'physical'
        task_id = 'Ep';
    case 'mental'
        task_id = 'Em';
end

%% initialize variables of interest
conNamesPerTaskRun = {};
conVec = [];
n_conPerTaskRun = 0;

%% fixation cross
if GLMprm.model_onset.(task_id).cross == 1
    conNamesPerTaskRun{1} = 'cross';
    n_conPerTaskRun = 
    conVec
end

end % function