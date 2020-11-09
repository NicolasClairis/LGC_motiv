function [var_binned1_low, var_binned1_high,...
    Bin_val_low, Bin_val_high,...
    var_binned2_low, var_binned2_high,...
    var_to_bin1_ascOrder_low, var_to_bin1_ascOrder_high,...
    bin_var_ascOrder_low, bin_var_ascOrder_high,...
    var_to_bin2_ascOrder_low, var_to_bin2_ascOrder_high] = do_bin3(var_to_bin1, bin_var, nb_bins, var_to_bin2, medSpl_perBin_or_Global)
% [var_binned1_low, var_binned1_high,...
%     Bin_val_low, Bin_val_high,...
%     var_binned2_low, var_binned2_high,...
%     var_to_bin1_ascOrder_low, var_to_bin1_ascOrder_high,...
%     bin_var_ascOrder_low, bin_var_ascOrder_high,...
%     var_to_bin2_ascOrder_low, var_to_bin2_ascOrder_high] = do_bin3(var_to_bin1, bin_var, nb_bins, var_to_bin2)
%
% do_bin3 creates 'nb_bins' bins of equal size (same number of samples in each bin) according to 'bin_var' values
% used to bin 'var_to_bin1', 'var_to_bin2' and 'bin_var' as in do_bin2.m. Here, in addition, for each
% bin, a median split is made on the 'var_to_bin1' variable and data (var_to_bin1, var_to_bin2 and bin_var) are
% spread in a '_low' group or in a '_high' group depending on this.
%
% INPUTS
% var_to_bin1: vector variable that you want to bin (and used for the median
% split)
%
% bin_var: vector variable according to which you want to make your 'nb_bins' bins
%
% nb_bins: number of bins that you want
%
% var_to_bin2: second vector variable that you want to bin
%
% medSpl_perBin_or_Global
% (0) make a median split of 'var_to_bin1' for each 'bin_var' bin
% (1) make a global median split of 'var_to_bin1' and then make the bins on
% the 2 groups of data separated based on the 'var_to_bin1' median split
%
% OUTPUTS
% var_binned1_low & var_binned1_high: 'var_to_bin1' binned in 'nb_bins' bins according to
% 'bin_var' values and for each bin, the lowest (respectively the highest) var_to_bin1 values are
% extracted
%
% Bin_val_low & Bin_val_high: 'bin_var' binned in 'nb_bins' bins in ascending order
% separated according to 'var_to_bin' values
%
% var_binned2_low & var_binned2_high: 'var_to_bin2' binned in 'nb_bins' bins according to
% 'bin_var' values and for each bin, the lowest (respectively the highest) var_to_bin1 values are
% used to split the data into two groups
%
% var_to_bin1_ascOrder_low, var_to_bin1_ascOrder_high,...
%     bin_var_ascOrder_low, bin_var_ascOrder_high,...
%     var_to_bin2_ascOrder_low, var_to_bin2_ascOrder_high: these variables
%     contain the data before averaging per bin
%
% Created by Nicolas Clairis (based on Alizee do_bin.m function) -
% 04/01/2018
%
% See also do_bin2.m and do_bin.m

%% initialize output variables
[var_binned1_low, var_binned1_high,...
    Bin_val_low, Bin_val_high,...
    var_binned2_low, var_binned2_high] = deal(nan(1,nb_bins));

%% extract the number of samples per bin and check it's ok
if medSpl_perBin_or_Global == 0 %% case median split per bin
    
    n_samples_pBin = length(bin_var)/nb_bins; % get number of samples per bin
    
    % check that n_sample_pBin is an integer and that each bin contains enough
    % samples to be split into two groups
    if floor(n_samples_pBin) ~= n_samples_pBin
        error(['Please enter another value for nb_bins ',...
            'so that length(bin_var)/nb_bins is an integer.']);
    elseif n_samples_pBin == 1
        error(['Only one sample per bin doesn''t allow to split them into low/high values according to var_to_bin.',...
            ' Please choose a lower value for nb_bins so as to allow at least 2 samples/bin']);
    elseif mod(n_samples_pBin,2) ~= 0
        warning(['length(bin_var)/nb_bins is impair',...
            ' which implies that the ''_low'' variables will contain one more sample per bin',...
            ' than the ''_high'' variables. If you don''t mind, just ignore this message,',...
            ' otherwise, try another value for nb_bins to get an equal number of samples for each gtroup.'])
    end
    
elseif medSpl_perBin_or_Global == 1 % case median split global
    
     n_samples_pBin = (length(bin_var)/nb_bins)/2; % get number of samples per bin
     
     % check that n_sample_pBin is an integer
     if floor(n_samples_pBin) ~= n_samples_pBin
         error(['Please enter another value for nb_bins ',...
             'so that (length(bin_var)/nb_bins)/2 is an integer.']);
     elseif n_samples_pBin == 1
         warning('Only one sample per bin.');
     end
end

%% sort all variables according to ascending order 'bin_var' values
[bin_var_ascOrder, bin_var_ascOrder_idx] = sort(bin_var);
var_to_bin1_ascOrder = var_to_bin1(bin_var_ascOrder_idx);
var_to_bin2_ascOrder = var_to_bin2(bin_var_ascOrder_idx);

%% bin split: initiate variables without bins but just median split per bin separation
if medSpl_perBin_or_Global == 0
    [var_to_bin1_ascOrder_low, var_to_bin1_ascOrder_high,...
        bin_var_ascOrder_low, bin_var_ascOrder_high,...
        var_to_bin2_ascOrder_low, var_to_bin2_ascOrder_high] = deal(NaN(length(bin_var),1));

%% global median split: spread into two groups
elseif medSpl_perBin_or_Global == 1
    % extract index of lower/higher than median
    
    median_var_to_bin1  = median(var_to_bin1(~isnan(var_to_bin1))); % ignore NaN values
    low_var_to_bin      = var_to_bin1_ascOrder <= median_var_to_bin1;
    high_var_to_bin     = var_to_bin1_ascOrder > median_var_to_bin1;
    
    % make the median split of all the variables based on var_to_bin1
    % median
    var_to_bin1_ascOrder_low    = var_to_bin1_ascOrder(low_var_to_bin);
    var_to_bin1_ascOrder_high   = var_to_bin1_ascOrder(high_var_to_bin);
    bin_var_ascOrder_low        = bin_var_ascOrder(low_var_to_bin);
    bin_var_ascOrder_high       = bin_var_ascOrder(high_var_to_bin);
    var_to_bin2_ascOrder_low    = var_to_bin2_ascOrder(low_var_to_bin);
    var_to_bin2_ascOrder_high   = var_to_bin2_ascOrder(high_var_to_bin);
    
    % if necessary, complete with NaN for NaN trials in order that the
    % script works (nanmean will ignore them however)
    n_NaN_trials = sum(isnan(var_to_bin1));
    NaN_trials = sum(low_var_to_bin) - sum(high_var_to_bin);
    if n_NaN_trials > 0
        
        if NaN_trials > 0 % less samples in high samples
            % compensate for difference between high and low samples
            var_to_bin1_ascOrder_high   = [var_to_bin1_ascOrder_high;   NaN(NaN_trials,1)];
            bin_var_ascOrder_high       = [bin_var_ascOrder_high;       NaN(NaN_trials,1)];
            var_to_bin2_ascOrder_high   = [var_to_bin2_ascOrder_high;   NaN(NaN_trials,1)];
            
            if NaN_trials ~= n_NaN_trials % more NaN trials in high sample, but also some in low sample
                bonus_NaN_trials = (n_NaN_trials - NaN_trials)/2;
                var_to_bin1_ascOrder_high   = [var_to_bin1_ascOrder_high;   NaN(bonus_NaN_trials,1)];
                bin_var_ascOrder_high       = [bin_var_ascOrder_high;       NaN(bonus_NaN_trials,1)];
                var_to_bin2_ascOrder_high   = [var_to_bin2_ascOrder_high;   NaN(bonus_NaN_trials,1)];
                var_to_bin1_ascOrder_low    = [var_to_bin1_ascOrder_low;    NaN(bonus_NaN_trials,1)];
                bin_var_ascOrder_low        = [bin_var_ascOrder_low;        NaN(bonus_NaN_trials,1)];
                var_to_bin2_ascOrder_low    = [var_to_bin2_ascOrder_low;    NaN(bonus_NaN_trials,1)];
            end
            
        elseif NaN_trials < 0 % less samples in low samples
            % compensate for difference between low and high samples
            var_to_bin1_ascOrder_low    = [var_to_bin1_ascOrder_low;    NaN(NaN_trials,1)];
            bin_var_ascOrder_low        = [bin_var_ascOrder_low;        NaN(NaN_trials,1)];
            var_to_bin2_ascOrder_low    = [var_to_bin2_ascOrder_low;    NaN(NaN_trials,1)];
            
            if abs(NaN_trials) ~= n_NaN_trials % more NaN trials in low sample, but also some in high samples
                bonus_NaN_trials = (n_NaN_trials - NaN_trials)/2;
                var_to_bin1_ascOrder_low    = [var_to_bin1_ascOrder_low;    NaN(bonus_NaN_trials,1)];
                bin_var_ascOrder_low        = [bin_var_ascOrder_low;        NaN(bonus_NaN_trials,1)];
                var_to_bin2_ascOrder_low    = [var_to_bin2_ascOrder_low;    NaN(bonus_NaN_trials,1)];
                var_to_bin1_ascOrder_high   = [var_to_bin1_ascOrder_high;   NaN(bonus_NaN_trials,1)];
                bin_var_ascOrder_high       = [bin_var_ascOrder_high;       NaN(bonus_NaN_trials,1)];
                var_to_bin2_ascOrder_high   = [var_to_bin2_ascOrder_high;   NaN(bonus_NaN_trials,1)];
            end
            
        elseif NaN_trials == 0 % as many data missing in both low and high
            var_to_bin1_ascOrder_high   = [var_to_bin1_ascOrder_high;   NaN(n_NaN_trials/2,1)];
            bin_var_ascOrder_high       = [bin_var_ascOrder_high;       NaN(n_NaN_trials/2,1)];
            var_to_bin2_ascOrder_high   = [var_to_bin2_ascOrder_high;   NaN(n_NaN_trials/2,1)];
            var_to_bin1_ascOrder_low    = [var_to_bin1_ascOrder_low;    NaN(n_NaN_trials/2,1)];
            bin_var_ascOrder_low        = [bin_var_ascOrder_low;        NaN(n_NaN_trials/2,1)];
            var_to_bin2_ascOrder_low    = [var_to_bin2_ascOrder_low;    NaN(n_NaN_trials/2,1)];
        end
        
    end
end


%% extract values for each bin
for iBin = 1:nb_bins
    current_interval    = (1:n_samples_pBin) + n_samples_pBin*(iBin - 1);
    
    %% median split for each bin
    if medSpl_perBin_or_Global == 0
        % extract data on which to focus for this bin
        curr_bin_var     = bin_var_ascOrder(current_interval);
        curr_var_to_bin1 = var_to_bin1_ascOrder(current_interval);
        curr_var_to_bin2 = var_to_bin2_ascOrder(current_interval);
        
        % extract median
        med_var_to_bin = median(curr_var_to_bin1(~isnan(curr_var_to_bin1))); % ignore NaN values
        
        % split according to var_to_bin
        % extract index of low/high var_to_bin
        low_var_to_bin          = curr_var_to_bin1 <= med_var_to_bin;
        high_var_to_bin         = curr_var_to_bin1 > med_var_to_bin;
        
        % data averaged for each bin
        % extract var_to_bin1
        var_binned1_low(iBin)   = nanmean(curr_var_to_bin1(low_var_to_bin));
        var_binned1_high(iBin)  = nanmean(curr_var_to_bin1(high_var_to_bin));
        
        % extract bin_var
        Bin_val_low(iBin)       = nanmean(curr_bin_var(low_var_to_bin));
        Bin_val_high(iBin)      = nanmean(curr_bin_var(high_var_to_bin));
        
        % extract var_to_bin2
        var_binned2_low(iBin)   = nanmean(curr_var_to_bin2(low_var_to_bin));
        var_binned2_high(iBin)  = nanmean(curr_var_to_bin2(high_var_to_bin));
        
        % data separated but not averaged
        % extract var_to_bin1
        var_to_bin1_ascOrder_low(current_interval)    = curr_var_to_bin1.*low_var_to_bin;
        var_to_bin1_ascOrder_high(current_interval)   = curr_var_to_bin1.*high_var_to_bin;
        
        % extract bin_var
        bin_var_ascOrder_low(current_interval)        = curr_bin_var.*low_var_to_bin;
        bin_var_ascOrder_high(current_interval)       = curr_bin_var.*high_var_to_bin;
        
        % extract var_to_bin2
        var_to_bin2_ascOrder_low(current_interval)    = curr_var_to_bin2.*low_var_to_bin;
        var_to_bin2_ascOrder_high(current_interval)   = curr_var_to_bin2.*high_var_to_bin;
        
    %% global median split
    elseif medSpl_perBin_or_Global == 1
        
        % data averaged for each bin
        % extract var_to_bin1
        var_binned1_low(iBin)   = nanmean(var_to_bin1_ascOrder_low(current_interval));
        var_binned1_high(iBin)  = nanmean(var_to_bin1_ascOrder_high(current_interval));
        
        % extract bin_var
        Bin_val_low(iBin)       = nanmean(bin_var_ascOrder_low(current_interval));
        Bin_val_high(iBin)      = nanmean(bin_var_ascOrder_high(current_interval));
        
        % extract var_to_bin2
        var_binned2_low(iBin)   = nanmean(var_to_bin2_ascOrder_low(current_interval));
        var_binned2_high(iBin)  = nanmean(var_to_bin2_ascOrder_high(current_interval));
    end
end

%% for bin split, delete NaN values to keep only low/high trials
if medSpl_perBin_or_Global == 0
    % var_to_bin1
    var_to_bin1_ascOrder_low = var_to_bin1_ascOrder_low(~isnan(var_to_bin1_ascOrder_low));
    var_to_bin1_ascOrder_high = var_to_bin1_ascOrder_high(~isnan(var_to_bin1_ascOrder_high));
    
    % extract bin_var
    bin_var_ascOrder_low = bin_var_ascOrder_low(~isnan(bin_var_ascOrder_low));
    bin_var_ascOrder_high = bin_var_ascOrder_high(~isnan(bin_var_ascOrder_high));
    
    % extract var_to_bin2
    var_to_bin2_ascOrder_low = var_to_bin2_ascOrder_low(~isnan(var_to_bin2_ascOrder_low));
    var_to_bin2_ascOrder_high = var_to_bin2_ascOrder_high(~isnan(var_to_bin2_ascOrder_high));

end

end