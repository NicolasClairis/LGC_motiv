function[z_diff, pvalue] = compare_indpdt_correlations(r1, r2, n1, n2)
% [z_diff, pvalue] = compare_indpdt_correlations(r1, r2, n1, n2)
% compare_indpdt_correlations aims at comparing two correlation coefficients using
% a Fisher z-transformation test. Importantly, the correlations have to be
% coming from two independent groups:
% ex: A<=>B (rAB) vs C<=>D (rCD) and compare_indpdt_correlations will
% compare whether rAB is significantly different from rCD.
% However, if some groups are dependent or overlapping (ex1: if you want to
% compare A1<=>B1 (rAB1) to A2<=>B2 (rAB2) (same groups but before/after treatment)
% or ex2: if you have one overlapping variable and want to compare A<=>B vs
% A<=>C (rAB vs rAC) then you should use a different test.
%
% INPUTS
% r1: first correlation coefficient
%
% r2: second correlation coefficient
%
% n1: number of subjects for obtaining r1
%
% n2: number of subjects for obtaining r2
%
% OUTPUTS
% z_diff: zscore for the difference
%
% pvalue: pvalue for the difference based on a two-tailed test

%% 1) transform the correlation coefficients into zscores
z1 = (1/2).*log((1 + r1)./(1-r1));
z2 = (1/2).*log((1 + r2)./(1-r2));

%% 2) compute the Standard Errors
SE1 = 1./sqrt(n1 - 3);
SE2 = 1./sqrt(n2 - 3);

%% 3) compute the combined Standard Error
SEdiff = sqrt((SE1.^2) + (SE2.^2));

%% 4) compute the zscore for the difference
z_diff = (z1 - z2)./SEdiff;

%% 5) calculate the p.value for a two-tailed test
p_value = 2 * (1 - normcdf(abs(z_diff)));

%% Output results
fprintf('z-score for difference: %.3f\n', z_diff);
fprintf('p-value for difference: %.3f\n', p_value);

end % function