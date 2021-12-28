function [standError] = sem(array, dimension)
%[standError] = sem(array, dimension)
% sem computes the standard error of the mean of a given vector
% or array (works with NaNs too)
%
% INPUT
% - array: (self-explanatory) n*p array of data
% - dimension: (optional, default = 1) dimension along which you want to compute sem
%
% OUTPUT :
% sem value(s)
%
% Originally developed by Emmanuelle Bioud - 24/01/2018



try
    %  The line that does it all:
    standError = std(array, 0, dimension,'omitnan')./sqrt(sum(~isnan(array), dimension));
    
catch ME
    
    %%  Take care of possible issues with the 'dimension' argument
    if nargin < 2
        warning('''dimension'' argument not specified -> Set to 1 or relevant dimension if 1D vector.')
    elseif isempty(dimension)
        warning('''dimension'' argument not specified -> Set to 1 or relevant dimension if 1D vector.')
    elseif ~isnumeric(dimension)
        warning('''dimension'' argument should be numeric. Can only be [1; 2 ; 3]. -> Set to 1 or relevant dimension if 1D vector.')
    elseif floor(dimension) ~= dimension % not an integer
        warning('''dimension'' argument should be an integer. Can only be [1; 2 ; 3]. -> Set to 1 or relevant dimension if 1D vector.')
    elseif ~ismember(dimension, [1 2 3])
        warning('Requested dimension is invalid. Can only be [1; 2 ; 3]. -> Set to 1 or relevant dimension if 1D vector.')
    else
        getReport(ME)
    end
    
    
    %% Define default dimension
    if isvector(array)
        [~, dimension] = max(size(array));
    else
        dimension = 1;
    end
    
    %% Try again the computation
    standError = std(array, 0, dimension,'omitnan')./sqrt(sum(~isnan(array), dimension));
    
    
end

end