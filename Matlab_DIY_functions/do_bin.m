function [var_binned, Bin_val] = do_bin(var_to_bin, scale, bin_var)
%[var_binned, Bin_val] = do_bin(var_to_bin, scale, bin_var)
% function to create bins according to a specific range of 'bin_var' values
% entered in 'scale' used to bin 'var_to_bin' and 'bin_var'.
%
% INPUTS
% var_to_bin: vector variable that you want to bin
%
% scale: vector containg the range of values on which you want to make
% your bins
% 'scale = -N1:u:N2' will create one datapoint for each step
% ([-N1,-N1+u], [-N1+u,-N1+2*u], etc. until [N2-u,N2]) values of bin_var
% and extract the corresponding mean of var_to_bin
%
% bin_var: vector variable according to which you want to make your bins
% based on scale
%
% OUTPUTS
% var_binned: 'var_to_bin' binned according to 'bin_var' and 'scale'
%
% Bin_val: 'bin_var' binned according to 'scale'
%
% Created by Alizee Lopez (2016) - adapted by Nicolas Clairis 04/01/2018

step = scale(2) - scale(1);

[Bin_val, var_binned] = deal(nan(1,length(scale)));

v = 0;
for vr_ind = scale
    v = v + 1;
    lower_bound = bin_var >= vr_ind;
    higher_bound = bin_var < vr_ind + step;
    current_interval = lower_bound & higher_bound;
    Bin_val(v)      = nanmean(bin_var(current_interval));
    var_binned(v)   = nanmean(var_to_bin(current_interval));
end

end