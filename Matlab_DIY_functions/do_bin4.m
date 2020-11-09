function [ var_to_bin_binned, bin_var1_binned, bin_var2_binned ] = do_bin4( var_to_bin, bin_var1, bin_var2, n_bins1, n_bins2 )
%[ var_to_bin1_binned, bin_var1_binned, bin_var2_binned ] = do_bin4( var_to_bin1, bin_var1, bin_var2, n_bins1, n_bins2 )
% do_bin4 splits the data first according to bin_var1 in n_bins1 bins. Then
% for each set of data, it creates bins based on bin_var2. In the output,
% you then get 3D-data.
%
% INPUT
% var_to_bin: variable to bin
%
% bin_var1: variable to create the first split
%
% bin_var2: variable to make bins after splitting the data based on
% bin_var1
%
% OUTPUT
% var_to_binned: var_to_bin binned according to bin_var1 and bin_var2
% Each line is a different bin of bin_var1 (in an ascending order).
% Each column is a different bin of bin_var2 (in an ascending order).
%
% bin_var1_bin_binned, bin_var2_binned: bin_var1 and bin_var2 re-ordered
% similarly as var_to_bin_binned
%
% Created by Nicolas Clairis -
% 20/11/2018
%
% See also do_bin.m, do_bin2.m and do_bin3.m

%% initialize the variables to fill
[var_to_bin_binned, bin_var1_binned, bin_var2_binned] = deal(nan(n_bins1,n_bins2));

%% by default, the scripts does not work if the number of bins entered in the input is not compatible with having exactly the same number of samples for each bin
% if you set 'same_n_samples_per_bin' to zero, it will ignore this and just
% use the best approximation possible

if ~exist('same_n_samples_per_bin','var') || isempty(same_n_samples_per_bin) || same_n_samples_per_bin == 1
    n_samples_pBin1 = length(bin_var1)/n_bins1; % get number of samples per bin
    n_samples_pBin2 = length(bin_var2)/(n_bins2*n_bins1); % get number of samples per bin
    % check that n_sample_pBin is an integer
    if floor(n_samples_pBin1) ~= n_samples_pBin1
        error(['Please enter another value for n_bins1 ',...
            'so that length(bin_var)/nb_bins is an integer.']);
    end
    % check that n_sample_pBin is an integer
    if floor(n_samples_pBin2) ~= n_samples_pBin2
        error(['Please enter another value for n_bins2 ',...
            'so that length(bin_var)/nb_bins is an integer.']);
    end
else
    n_samples_pBin1 = floor(length(bin_var1)/n_bins1); % get number of samples per bin
end

%% 1) extract the ascending order of bin_var1
[~, bin_var1_ascOrder_idx] = sort(bin_var1);


%% 2) bin the data based on bin_var2 for each bin of bin_var1
% loop through bin_var1 bins
for iBin1 = 1:n_bins1
    current_interval1    = (1:n_samples_pBin1) + n_samples_pBin1*(iBin1 - 1);
    
    % average bin_var1
    bin_var1_tmp        = nanmean( bin_var1( bin_var1_ascOrder_idx(current_interval1) ) );
    
    % extract corresponding values for var_to_bin and bin_var2
    var_to_bin_tmp      = var_to_bin( bin_var1_ascOrder_idx(current_interval1) );
    bin_var2_tmp        = bin_var2( bin_var1_ascOrder_idx(current_interval1) );
    
    % bin according to ascending order bin_var2_tmp values
    [bin_var2_tmp_ascOrder, bin_var2_tmp_ascOrder_idx] = sort(bin_var2_tmp);
    var_to_bin_tmp_ascOrder = var_to_bin_tmp(bin_var2_tmp_ascOrder_idx);
    
    % bin according to bin_var2
    for iBin2 = 1:n_bins2
        current_interval2    = (1:n_samples_pBin2) + n_samples_pBin2*(iBin2 - 1);
        bin_var1_binned(iBin1, iBin2)       = bin_var1_tmp; % (first bin type, should be constant for all
    % bin_var2)
        bin_var2_binned(iBin1, iBin2)       = nanmean(bin_var2_tmp_ascOrder(current_interval2));
        var_to_bin_binned(iBin1, iBin2)    = nanmean(var_to_bin_tmp_ascOrder(current_interval2));
    end
    
end


end

