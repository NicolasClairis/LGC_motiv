function [cleaned_X, idx_badSubs, cleaned_X_bis] = rmv_outliers_3sd(variable_X)
%[cleaned_X, idx_badSubs, cleaned_X_bis] = rmv_outliers_3sd(variable_X)
% rmv_outliers_3sd will remove outliers (>mean+3SD or <mean-3SD) from
% variable_X and give you the corresponding index. It will also remove any
% NaN and include the index as well.
%
% INPUTS
% variable_X: 1*n1 vector with the variable to clean
%
% OUTPUTS
% cleaned_X: 1*n2 vector with the cleaned variable (n2=n1-any bad sample)
%
% idx_badSubs: 1*n1 vector indicating the index, in variable_X, of the bad
% subjects
%
% cleaned_X_bis: 1*n1 vector where any outlier is replaced by a NaN

%% check if there are NaNs
idx_badSubs = isnan(variable_X);

%% extract mean and SD
mu_X = median(variable_X, 2, 'omitnan');
std_X = std(variable_X, [],2, 'omitnan');

%% extract outliers index
too_high = variable_X > mu_X + 3.*std_X;
too_low = variable_X < mu_X - 3.*std_X;
idx_badSubs(too_high | too_low) = true;

%% remove outliers
cleaned_X = variable_X(idx_badSubs == false);
cleaned_X_bis = variable_X;
cleaned_X_bis(idx_badSubs == true) = NaN;

end % function