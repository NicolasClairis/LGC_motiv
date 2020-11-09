function[param] = fn_for_prior(X, prm_priors)
%[param] = fn_for_prior(prior, prm_priors)
% fn_for_prior transforms X into exp(X) or keeps it just as X depending on
% prm_priors ('free'/'pos')
%
% INPUTS
% X: prior to transform
%
% prm_priors = parameter for the prior to know how to transform X
%   'free': no constraint on the prior => keeps param = X
%   'pos': prior is constrained to be positive => param = exp(X)
%
% OUTPUTS
% param: X after being transformed

switch prm_priors
    case 'pos' % positive constraint on prior => use of exponential
        param = exp(X);
    case 'free' % no constraint => use directly the variable
        param = X;
end

end % function