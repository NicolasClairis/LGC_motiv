function[betas, pval, xSort, yFit] = linear_fit(x, y, cstt)
% [betas, pval, xSort, yFit] = linear_fit(x, y, cstt)
% linear_fit will extract the linear fit beta and p.value for the
% significativity of the betas using the glmfit.m function. It will also
% provide you with the fit in case you want to display the output
%
% INPUT
% x: n*nT matrix with n = the number of variables for the linear fit and nT
% = the number of trials/observations to make the fit onto
%
% y: 1*nT vector with output to predict
%
% cstt: include a constant ('on' by default) or no constant ('off')
%
% OUTPUT
% betas: 1*n vector with beta corresponding to each regressor
%
% pval: 1*n vector indicating the p.value for each beta
%
% xSort: n*nT matrix
%
% yFit: fitted y values with the linear fit
%
% See also glmfit.m

%% check data is in the correct format
% add constant by default
if ~exist('cstt','var') || isempty(cstt) || ~ismember(cstt,{'on','off'})
    cstt = 'on';
    disp('Linear model will include a constant by default');
end

% flip y if necessary
if size(y,1) > 1 && size(y,2) == 1
    y = y';
end
nT = size(y,2);

% flip x if necessary
if size(x,1) > 1 && (size(x,2) == 1 || size(x,2) ~= nT)
    x = x';
end
% recheck size x is ok
if size(x,2) ~= nT
    error(['problem with X vs Y dimensions. ',...
        'Please enter y in a 1*nTrials and x in (nVars)*(nTrials) format']);
end

%% perform the fit
[betas, ~, stats] = glmfit(x, y, 'normal','Constant',cstt);
pval = stats.p;

%% extract the fit
if size(x,1) == 1
    xSort = sort(x);
else
    xSort = x; % no sorting if more than 1 predicting variable
end
yFit = glmval(betas, xSort, 'identity');

end % function