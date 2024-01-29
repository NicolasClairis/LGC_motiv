function [cleaned_X, idx_badSubs, cleaned_X_bis, idx_goodSubs] = rmv_outliers_3sd(variable_X)
%[cleaned_X, idx_badSubs, cleaned_X_bis, idx_goodSubs] = rmv_outliers_3sd(variable_X)
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
%
% idx_goodSubs: 1*n1 vector indicating the index, in variable_X, of the
% good subjects

%% check that input is ok in terms of lines/columns and adapt it if not
if size(variable_X,1) > 1 && size(variable_X,2) == 1
    % need to flip X to be in the correct range
    variable_X = variable_X';
elseif size(variable_X,1) > 1 && size(variable_X,2) > 1
    error('variable_X is not a vector but a matrix');
elseif size(variable_X,1) == 0 && size(variable_X,2) == 0
    error('variable_X is empty');
end

%% check if there are NaNs
idx_badSubs = isnan(variable_X);
idx_goodSubs = ~isnan(variable_X);

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