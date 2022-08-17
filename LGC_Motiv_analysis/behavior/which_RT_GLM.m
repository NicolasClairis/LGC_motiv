function[GLMprm] = which_RT_GLM(GLM)
% [GLMprm] = which_RT_GLM(GLM)
% which_RT_GLM will specify some general parameters for the GLM to perform
% with reaction times.
%
% INPUTS
% GLM: GLM number
%
% OUTPUTS
% GLMprm: structure with GLM parameters
%   .main: main parameters
%       .RT_format: 'raw'/'zscore'/'log'
%   .regs: which regressors to include
%       .run_cstt: one constant per run
%       .task_cstt: one constant per task (pooling runs of each task)
%       .uncertainty: uncertainty derived from the model
%       .uncertaintyRtg: binary uncertainty rating based on button pressed
%       .deltaMoney: difference in terms of money between highE and default
%       option
%       .deltaEffort: level of effort proposed
%       .trialN: trial number
%       .RT_prevTrial: RT at previous trial

%% by default all vars are off
% general parameters
GLMprm.main.RT_format = 'raw';
potentialRegressors = {'run_cstt','task_cstt','uncertainty','uncertaintyRtg',...
    'deltaMoney','deltaEffort','trialN','RT_prevTrial'};
for iReg = 1:length(potentialRegressors)
    GLMprm.regs.(potentialRegressors{iReg}) = 'off';
end % regressor loop

%% adapt variables depending on the selected GLM
switch GLM
    case 1
        % general
        GLMprm.main.RT_format = 'raw';
        % which regressors to include
        GLMprm.regs.run_cstt = 'off';
        GLMprm.regs.task_cstt = 'on';
        GLMprm.regs.uncertainty = 'on';
        GLMprm.regs.uncertaintyRtg = 'off';
        GLMprm.regs.deltaMoney = 'on';
        GLMprm.regs.deltaEffort = 'on';
        GLMprm.regs.trialN = 'on';
        GLMprm.regs.RT_prevTrial = 'on';
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