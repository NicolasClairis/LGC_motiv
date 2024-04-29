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
%       .RT_format: 'raw'/'zscorePerRun'/'zscoreARuns'/'zscorePerRunARuns'/'log'
%       .orth_vars: orthogonalize the regressors or not? ('on'/'off'), by
%       default no orthogonalization will be applied
%       .confMdlType: use 'bayesian' or 'simple' model for extracting
%       confidence
%       .confMdlN: string indicating which model number to consider
%
%   .regs: which regressors to include ('on'/'off')
%       .run_cstt: one constant per run
%       .task_cstt: one constant per task (pooling runs of each task)
%       .conf: confidence derived from the model
%       .confRtg: binary confidence rating based on button pressed
%       .deltaMoney: difference in terms of money between highE and default
%       option
%       .deltaEffort: level of effort proposed
%       .trialN: trial number
%       .RT_raw_prevTrial: RT (raw in seconds) at previous trial

%% by default all vars are off
% general parameters
GLMprm.main.RT_format = 'raw';
GLMprm.main.orth_vars = 'off';
GLMprm.main.confMdlType = 'simple';
GLMprm.main.confMdlN = '';
potentialRegressors = {'run_cstt','task_cstt',...
    'conf','confRtg',...
    'RP','deltaMoney','deltaR','deltaP',...
    'deltaEffort','deltaEp','deltaEm',...
    'trialN','physical_Fatigue','mental_Facilitation','RT_raw_prevTrial'};
for iReg = 1:length(potentialRegressors)
    GLMprm.regs.(potentialRegressors{iReg}) = 'off';
end % regressor loop

%% adapt variables depending on the selected GLM
switch GLM
    case 1
        % general
        GLMprm.main.RT_format = 'raw';
        GLMprm.main.orth_vars = 'off';
        GLMprm.main.confMdlType = 'simple';
        GLMprm.main.confMdlN = '3';
        % which regressors to include
        GLMprm.regs.run_cstt = 'off';
        GLMprm.regs.task_cstt = 'on';
        GLMprm.regs.conf = 'on';
        GLMprm.regs.confRtg = 'off';
        GLMprm.regs.deltaMoney = 'on';
        GLMprm.regs.deltaEffort = 'on';
        GLMprm.regs.trialN = 'on';
        GLMprm.regs.RT_raw_prevTrial = 'on';
    case 2
        % general
        GLMprm.main.RT_format = 'raw';
        GLMprm.main.orth_vars = 'off';
        GLMprm.main.confMdlType = 'bayesian';
        GLMprm.main.confMdlN = '3';
        % which regressors to include
        GLMprm.regs.run_cstt = 'on';
        GLMprm.regs.task_cstt = 'off';
        GLMprm.regs.conf = 'on';
        GLMprm.regs.RP = 'off';
        GLMprm.regs.deltaR = 'on';
        GLMprm.regs.deltaP = 'on';
        GLMprm.regs.deltaEp = 'on';
        GLMprm.regs.deltaEm = 'on';
        GLMprm.regs.physical_Fatigue = 'on';
        GLMprm.regs.mental_Facilitation = 'on';
    case 3 % same as GLM2 but adding RP constant
        % general
        GLMprm.main.RT_format = 'raw';
        GLMprm.main.orth_vars = 'off';
        GLMprm.main.confMdlType = 'bayesian';
        GLMprm.main.confMdlN = '3';
        % which regressors to include
        GLMprm.regs.run_cstt = 'on';
        GLMprm.regs.task_cstt = 'off';
        GLMprm.regs.conf = 'on';
        GLMprm.regs.RP = 'on';
        GLMprm.regs.deltaR = 'on';
        GLMprm.regs.deltaP = 'on';
        GLMprm.regs.deltaEp = 'on';
        GLMprm.regs.deltaEm = 'on';
        GLMprm.regs.physical_Fatigue = 'on';
        GLMprm.regs.mental_Facilitation = 'on';
    case 4
        % general
        GLMprm.main.RT_format = 'raw';
        GLMprm.main.orth_vars = 'off';
        GLMprm.main.confMdlType = 'bayesian';
        GLMprm.main.confMdlN = '3';
        % which regressors to include
        GLMprm.regs.run_cstt = 'on';
        GLMprm.regs.task_cstt = 'off';
        GLMprm.regs.conf = 'on';
        GLMprm.regs.RP = 'on';
        GLMprm.regs.deltaR = 'on';
        GLMprm.regs.deltaP = 'on';
        GLMprm.regs.deltaEp = 'on';
        GLMprm.regs.deltaEm = 'on';
        GLMprm.regs.trialN = 'on';
    otherwise
        error(['RT GLM number ',num2str(GLM),' does not exist yet.']);
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