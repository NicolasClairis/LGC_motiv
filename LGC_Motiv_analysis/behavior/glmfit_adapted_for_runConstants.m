function[betas_struct, Y_fit, betas] = glmfit_adapted_for_runConstants(r1_cstt, r2_cstt,...
    X_regs_bf_r_cstt, Y_var, x_reg_bf_r_cstt_names, modelType)
% [betas_struct, Y_fit, betas] = glmfit_adapted_for_runConstants(r1_cstt, r2_cstt,...
%   X_regs_bf_r_cstt, Y_var, x_reg_bf_r_cstt_names, modelType)
% glmfit_adapted_for_runConstants will perform a glmfit between variable X
% and Y for the current task (physical or mental) but after taking into 
% account that for the current subject maybe not all runs are included
% => The script will automatically remove the runs not included to avoid
% making glmfit not work properly.
%
% INPUTS
% r1_cstt: 1*n constant equal to 1 for the run 1, 0 on session 2 and NaN for
% any trial not to be included.
%
% r2_cstt: 1*n constant equal to 0 on session 1, 1 for the run 2, and NaN for
% any trial not to be included.
%
% X_regs_bf_r_cstt: nTrials*nRegressors matrix containing all the
% regressors to be included in the GLM apart from run constants
%
% Y_var: variable on which you want to do the fit
%
% x_reg_bf_r_cstt_names: 1*nRegressors cell including name of each of the
% regressors to be included (apart for run constants)
%
% modelType: string indicating which model to use
% ('normal'/'binomial'/etc.)
%
% OUTPUTS
% betas_struct: structure with one subfield for each beta
%
% Y_fit: fit for Y_var with glmfit
%
% betas: n*1 vector with the betas output from glmfit (n will depend on the
% number of regressors included, ie the size will vary depending on if all
% runs are ok or only one is ok).
%
% Note: glmfit and glmval require the Matlab Statistics and Machine Learning Toolbox


%% check if runs are ok
r1_ok = sum(r1_cstt,'omitnan') > 0;
r2_ok = sum(r2_cstt,'omitnan') > 0;

%% build regressors based on that + extract index for the fit
if r1_ok == 1 && r2_ok == 1
    x_regs = [r1_cstt, r2_cstt, X_regs_bf_r_cstt];
elseif r1_ok == 1 && r2_ok == 0
    x_regs = [r1_cstt, X_regs_bf_r_cstt];
elseif r1_ok == 0 && r2_ok == 1
    x_regs = [r2_cstt, X_regs_bf_r_cstt];
else
    error('wtf happened? No runs seems to be ok');
end

%% define which model to use
if ~exist('modelType','var') || isempty(modelType)
    modelType = 'normal';
    modelType_for_glmval = 'identity';
else
    switch modelType
        case 'normal'
            modelType_for_glmval = 'identity';
        case 'binomial'
            modelType_for_glmval = 'logit';
        otherwise
            error(['model ',modelType,' not ready yet, please update script.']);
    end
end

%% perform glmfit + corresponding fit
betas = glmfit(x_regs, Y_var, modelType,'Constant','off');
Y_fit = glmval(betas, x_regs, modelType_for_glmval,'Constant','off');

%% distribute betas according to how many runs were used
% attribute values for each run depending on if run included or not
jReg = 0;
if r1_ok == 1 && r2_ok == 1
    jReg = jReg + 1;
    betas_struct.r1_cstt = betas(jReg);
    jReg = jReg + 1;
    betas_struct.r2_cstt = betas(jReg);
elseif r1_ok == 1 && r2_ok == 0
    jReg = jReg + 1;
    betas_struct.r1_cstt = betas(jReg);
    betas_struct.r2_cstt = NaN;
elseif r1_ok == 0 && r2_ok == 1
    betas_struct.r1_cstt = NaN;
    jReg = jReg + 1;
    betas_struct.r2_cstt = betas(jReg);
else
    error('wtf happened? No runs seems to be ok');
end

% attribute other regressors
for iReg = 1:size(X_regs_bf_r_cstt,2)
    jReg = jReg + 1;
    betas_struct.(x_reg_bf_r_cstt_names{iReg}) = betas(jReg);
end

end % function