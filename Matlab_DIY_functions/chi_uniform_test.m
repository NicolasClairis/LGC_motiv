function [ proba ] = chi_uniform_test( X, n_bins )
% [ proba ] = chi_uniform_test( X, n_bins )
% chi_uniform_test tests whether the data in X is uniform or not.
%
% INPUTS
% X: data 1*n_bins with the data to test. Each column should contain the
% number of observations for the corresponding bin.
%
% n_bins number of bins in X
%
% OUTPUTS
% proba: probability that X is uniform
%
% Function initially written by Nicolas Clairis and
% inspired by the great Jules Brochard - june 2019.

%% check data entered correctly
if size(X,1) > 1 || size(X,2) == 1
    error('problem with X, should be entered as 1*n_bin vector. Please fix it and come back');
end

%% get number of bins
if ~exist('n_bins','var') || isempty(n_bins)
    n_bins = length(X);
end

%% extract theoretical values for uniform distribution to get chi-square
n_total_samples = sum(X);
pX = X./n_total_samples; % pass from frequency for each bin to %
ref_value_uniform = 1/n_bins;%n_total_samples/n_bins;

ref_uni_X = ones(1,n_bins)*ref_value_uniform;

bins_ok = (ref_uni_X ~= 0); % cannot divide by zero
n_bins_ok = sum(bins_ok);

chi_square = n_bins_ok*sum( ((pX(bins_ok) - ref_uni_X(bins_ok)).^2)./ref_uni_X(bins_ok) );

%% extract p.value
proba = 1 - chi2cdf(chi_square, n_bins_ok);
disp(['Probability that the distribution is uniform (null hypothesis): ',num2str(proba)])

end % function