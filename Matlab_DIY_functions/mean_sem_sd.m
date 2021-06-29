function[meanVar, semVar, sdVar]=mean_sem_sd(var, dim)
% [meanVar, semVar, sdVar]=mean_sem_sd(var, dim)
% mean_sem_sd will get you the mean, the standard error of the mean and the
% standard deviation of the variable var along the dimension specified in
% dim.
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

meanVar = mean(var, dim);
semVar = sem(var, dim);
sdVar = std(var, 0, dim);

end % function