function[meanVar, semVar, sdVar, medianVar]=mean_sem_sd(var, dim)
% [meanVar, semVar, sdVar, medianVar]=mean_sem_sd(var, dim)
% mean_sem_sd will get you the mean, the standard error of the mean, the
% standard deviation and the median of the variable var along the dimension
% specified in dim.
%
% INPUTS
% var: vector/matrix with the variable of interest
%
% bin: dimension along which the mean, sem and sd will be extracted
%
% OUTPUTS
% meanVar: mean of var along dim
%
% semVar: standard error of the mean along dim
%
% sd: standard error deviation along dim
%
% medianVar: median of var along dim

if ~exist('dim','var') || isempty(dim)
    error('Please provide a value for ''dim''.');
end

%% extract mean, SEM and SD
meanVar = mean(var, dim,'omitnan');
semVar = sem(var, dim);
sdVar = std(var, 0, dim,'omitnan');
medianVar = median(var, dim, 'omitnan');

end % function