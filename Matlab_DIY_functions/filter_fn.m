function[goodS_idx] = filter_fn(raw_or_corr_nm, x, y)
% [goodS_idx] = filter_fn(raw_or_corr_nm, x, y)
% filter_fn will filter the subjects who are outliers. If 'raw' is entered
% in inputs, we remove only those who have no values (NaN). If 'filtered'
% is entered, the script will remove any outlier based on mean +/- 3*SD.
%
% INPUTS
% raw_or_corr_nm: variable indicating how to filter
% 'raw': remove only those who have no values (NaN)
% 'filtered': remove any outlier based on mean +/- 3*SD
%
% x, y: variables to filter, as 1*nObservations vector
%
% OUTPUTS
% goodS_idx: 1*N vector with indication towards which subjects to keep or
% not as a binary variable

switch raw_or_corr_nm
    case 'raw'
        goodS_idx = ~isnan(x.*y);
    case 'filtered'
        [~,~,filtered_x] = rmv_outliers_3sd(x);
        [~,~,filtered_y] = rmv_outliers_3sd(y);
        goodS_idx = ~isnan(filtered_x.*filtered_y);
end
            
end % sub-function