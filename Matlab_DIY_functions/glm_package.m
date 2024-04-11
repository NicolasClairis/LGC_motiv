function[r_corr, betas, pval, y_fit, x_sorted, y_fit_xSorted, z_pval_corr] = glm_package(x_var, y_var, distrib, constant)
% [r_corr, betas, pval, y_fit, x_sorted, y_fit_xSorted, z_pval_corr] = glm_package(x_var, y_var, distrib, constant)
% glm_package will perform all the package required by a classic GLM to
% have all the data you may want in output. In other words, it will provide
% you with correlation coefficients, betas, p.values and the corresponding
% fit.
%
% INPUTS
% x_var: nTrials*nVariables vector with predictor variables
%
% y_var: nTrials*1 response variable (dependent variable)
%
% distrib: distribution to test 'normal' by default
%
% constant: indicate if you want to include ('on') or not ('off'). 'on' by
% default if not specified
%
% OUTPUTS
% r_corr: vector with correlation coefficients of size (nVariables*1) if
% constant 'off' or ((1+nVariables)*1) if constant 'on' and first variable
% = constant in that case
%
% betas: vector with betas of size (nVariables*1) if constant 'off' or 
% ((1+nVariables)*1) if constant 'on' and first variable = constant in that
% case
%
% pval: vector with p.values of size (nVariables*1) if constant 'off' or 
% ((1+nVariables)*1) if constant 'on' and first variable = constant in that
% case
%
% y_fit: nTrials*1 vector with fit correspond to x_var values
%
% x_sorted: x_var sorted in ascending order (only if there is only 1
% variable in X)
%
% y_fit_xSorted: corresponding fit for x_sorted (to make figure displays
% easier) (only if there is only 1 variable in X)
%
% z_pval_corr: p.value for the correlation coefficients in r_corr
%
% See also glmfit.m and glmval.m
% This script requires the Statistics and Machine Learning Toolbox
% https://www.mathworks.com/products/statistics.html

%% control for correct orientation of variables of interest
% x_var should be also in the following format: nTrials*nVariables
if isempty(x_var)
    error('x_var is empty');
elseif size(x_var,2) > size(x_var,1)
    if  size(x_var,1) == 1 % flip X if entered as 1*nTrials instead of nTrials*1
        x_var = x_var';
    elseif  size(x_var,1) > 1 % ask the experimenter to modify X in that more complex case just to be sure all is good
        error(['Careful X should be entered as nTrials*nVars ',...
            'but it seems that there are more trials than variables ',...
            'in the entered input.']);
    end
end

% y_var has to be in the following format: nTrials*1
if isempty(y_var)
    error('y_var is empty');
elseif size(y_var,1) == 1 && size(y_var,2) > size(y_var,1) % flip y if entered as 1*nTrials instead
    y_var = y_var';
end

%% initialize inputs
% distribution
list_possible_distribs = {'normal','binomial',...
    'poisson','gamma','inverse gaussian'};
if ~exist('distrib','var') || ~ismember(distrib, list_possible_distribs)
    distrib = 'normal';
    warning(['distrib empty in inputs so set to ',distrib,' by default']);
end

% constant
if ~exist('constant','var') || ~ismember(constant, {'off','on'})
    constant = 'on';
end


%% perform GLM
[betas,~,stats] = glmfit(x_var, y_var, distrib, 'Constant',constant);
pval = stats.p;
%% extract correlation coefficients
if size(x_var,2) == 1
    [z_betas,~,stats_z_vars] = glmfit(nanzscore(x_var), nanzscore(y_var),distrib,'Constant',constant);
elseif size(x_var,2) > 1
    % if more then 1 x.variables => need to zscore each independently
    z_x_var = NaN(size(x_var));
    n_x_vars = size(x_var,2);
    for iX_var = 1:n_x_vars
        z_x_var(:,iX_var) = nanzscore(x_var(:,iX_var));
    end % loop through x_var variables
    
    [z_betas,~,stats_z_vars] = glmfit(z_x_var, nanzscore(y_var),distrib,'Constant',constant);
end
switch constant
    case 'off'
        r_corr = z_betas;
        z_pval_corr = stats_z_vars.p;
    case 'on' % ignore constant
        r_corr = z_betas(2:end);
        z_pval_corr = stats_z_vars.p(2:end);
end

%% extract corresponding fit
switch distrib
    case 'normal'
        distrib_fit = 'identity';
    case 'binomial'
        distrib_fit = 'logit';
    case 'poisson'
        distrib_fit = 'log';
    otherwise
        error(['case where distrib = ',distrib,' not ready yet']);
end
y_fit = glmval(betas, x_var, distrib_fit, 'Constant',constant);

%% if only 1 x.variable => sort to make fit nicer, otherwise just keep all
% the values
if size(x_var,2) == 1
    x_sorted = sort(x_var);
    y_fit_xSorted = glmval(betas, x_sorted, distrib_fit, 'Constant',constant);
else
    x_sorted = [];
    y_fit_xSorted = [];
end
end % function