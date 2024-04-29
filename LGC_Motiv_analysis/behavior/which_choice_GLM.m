function[GLMprm] = which_choice_GLM(GLM)
% [GLMprm] = which_choice_GLM(GLM)
% which_choice_GLM will specify some general parameters for the GLM to perform
% with choices.
%
% INPUTS
% GLM: GLM number
%
% OUTPUTS
% GLMprm: structure with GLM parameters
%   .main: main parameters
%       .orth_vars: orthogonalize the regressors or not? ('on'/'off'), by
%       default no orthogonalization will be applied
%
%   .regs: which regressors to include ('on'/'off')
%       .run_cstt: one constant per run
%       .task_cstt: one constant per task (pooling runs of each task)
%       .deltaMoney: difference in terms of money between highE and default
%       option
%       .deltaEffort: level of effort proposed
%       .trialN: trial number

%% by default all vars are off
% general parameters
GLMprm.main.orth_vars = 'off';
potentialRegressors = {'run_cstt','task_cstt','choice_bias',...
    'RP','deltaMoney','deltaR','deltaP',...
    'deltaEffort','deltaEp','deltaEm',...
    'trialN','physical_Fatigue','mental_Facilitation'};
for iReg = 1:length(potentialRegressors)
    GLMprm.regs.(potentialRegressors{iReg}) = 'off';
end % regressor loop

%% adapt variables depending on the selected GLM
switch GLM
    case 1
        % general
        GLMprm.main.orth_vars = 'off';
        % which regressors to include
        GLMprm.regs.run_cstt = 'off';
        GLMprm.regs.task_cstt = 'off';
        GLMprm.regs.choice_bias = 'on';
        GLMprm.regs.deltaR = 'on';
        GLMprm.regs.deltaP = 'on';
        GLMprm.regs.deltaEp = 'on';
        GLMprm.regs.deltaEm = 'on';
        GLMprm.regs.physical_Fatigue = 'on';
        GLMprm.regs.mental_Facilitation = 'on';
    otherwise
        error(['choice GLM number ',num2str(GLM),' does not exist yet.']);
end % GLM

%% check that all regressors exist by default
actualRegs = fieldnames(GLMprm.regs);
for iRegBis = 1:length(actualRegs)
    curr_reg_nm = actualRegs{iRegBis};
    isRegIncludedByDefault = sum(strcmp(curr_reg_nm,potentialRegressors));
    if isRegIncludedByDefault == 0
        error(['Please add ',curr_reg_nm,' to the list of potential ',...
            'regressors as it is currently missing and may create problems in the future']);
    end
end % regressor loop

end % function