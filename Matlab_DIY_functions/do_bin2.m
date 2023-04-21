function [var_binned, Bin_val, bin_idx] = do_bin2(var_to_bin, bin_var, nb_bins, same_n_samples_per_bin)
%[var_binned, Bin_val, bin_idx] = do_bin2(var_to_bin, bin_var, nb_bins, same_n_samples_per_bin)
% function to create 'nb_bins' bins of equal size (same number of samples in each bin) according to 'bin_var' values
% used to bin 'var_to_bin' and 'bin_var'.
%
% INPUTS
% var_to_bin: vector variable that you want to bin
%
% bin_var: vector variable according to which you want to make your 'nb_bins' bins
%
% nb_bins: number of bins that you want
%
% same_n_samples_per_bin: (optional) by default, the script does not work 
% if the number of bins entered in the input is not compatible with having 
% exactly the same number of samples for each bin if you set 
% 'same_n_samples_per_bin' to zero, it will ignore this and just
% use the best approximation possible
%
% OUTPUTS
% var_binned: 'var_to_bin' binned in 'nb_bins' bins according to 'bin_var'
%
% Bin_val: 'bin_var' binned in 'nb_bins' bins in ascending order
%
% bin_idx: vector with the index corresponding to each bin for all trials
%
% Created by Nicolas Clairis (based on Alizee do_bin.m function) -
% 04/01/2018
%
% See also do_bin.m and do_bin3.m

%% Remove NaN samples from the variable before binning
NaN_trials  = isnan(var_to_bin) | isnan(bin_var);
var_to_bin_noNaN = var_to_bin;
bin_var_noNaN = bin_var;
if sum(NaN_trials) > 0
    var_to_bin_noNaN(NaN_trials)  = [];
    bin_var_noNaN(NaN_trials)     = [];
end

%% initialize the variables to fill
[Bin_val, var_binned] = deal(nan(1,nb_bins));

%% by default, the scripts does not work if the number of bins entered in the input is not compatible with having exactly the same number of samples for each bin
% if you set 'same_n_samples_per_bin' to zero, it will ignore this and just
% use the best approximation possible

if ~exist('same_n_samples_per_bin','var') || isempty(same_n_samples_per_bin) || same_n_samples_per_bin == 1
    n_samples_pBin = length(bin_var_noNaN)/nb_bins; % get number of samples per bin
    % check that n_sample_pBin is an integer
    if floor(n_samples_pBin) ~= n_samples_pBin
        error(['Please enter another value for nb_bins ',...
            'so that length(bin_var after removing NaN trials)/nb_bins is an integer.']);
    end
    
else
    n_samples_pBin = floor(length(bin_var_noNaN)/nb_bins); % get number of samples per bin
end

%% bin according to ascending order based on bin_var_noNaN values
[bin_var_ascOrder, bin_var_ascOrder_noNaN_idx] = sort(bin_var_noNaN);
var_to_bin_ascOrder = var_to_bin_noNaN(bin_var_ascOrder_noNaN_idx);

%% extract index corresponding to each bin (but careful to keep track of NaN values)
bin_idx = NaN(1,length(bin_var_noNaN));
for iBin = 1:nb_bins
    current_interval = (1:n_samples_pBin) + n_samples_pBin*(iBin - 1);
    curr_bin_idx_tmp = bin_var_ascOrder_noNaN_idx(current_interval);
    bin_idx(curr_bin_idx_tmp) = iBin;
end % bin loop
% if there were NaN values, you need to correct the index
if sum(NaN_trials) > 0
    for iTrial = 1:length(NaN_trials)
        if NaN_trials(iTrial) == true
            bin_idx = [bin_idx(1:(iTrial - 1)), NaN, bin_idx(iTrial:end)];
        end
    end % trial loop
end

%% extract values for each bin
for iBin = 1:nb_bins
    current_interval    = (1:n_samples_pBin) + n_samples_pBin*(iBin - 1);
    Bin_val(iBin)       = mean(bin_var_ascOrder(current_interval),'omitnan');
    var_binned(iBin)    = mean(var_to_bin_ascOrder(current_interval),'omitnan');
end % bin loop

end