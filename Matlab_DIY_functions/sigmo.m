function[y] = sigmo(x, lambda)
% [y]=sigmo(x, lambda)
% sigmo will compute the sigmoid of x eventually multiplied by a parameter
% lambda (inverse temperature)
% y = 1/(1 + exp(-lambda*x)
%
% If lambda left empty, will use lambda = 1 by default
%
% INPUTS
% x: the variable passed through the sigmoid
%
% lambda: the multiplying variable of x (=inverse temperature)
%
% OUTPUTS
% y: the output, sigmoid of x

%% define lambda
if ~exist('lambda','var') || isempty(lambda)
    lambda = 1;
end

%% compute the sigmoid
y = 1./( 1 + exp( -lambda.*x) );

end % function end