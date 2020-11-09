function [ X_smoothed ] = smooth_NC( X_data, smooth_kernel )
%[ X_smoothed ] = smooth_NC( X_data, smooth_kernel ) smooths the data
%contained in X_data with a smooth_kernel size.
%
% INPUTS
% X_data: (1xN) vector of the data to smooth (data will be flipped if
% entered as a (Nx1) vector)
%
% smooth_kernel: kernel size for the smoothing
%
% OUTPUT
% X_smoothed: (1xN) vector containing the X_data smoothed
%
% See also conv
% Written by Nicolas Clairis - february 2019

%% check X_data format
if size(X_data,1) > size(X_data,2)
    X_data = X_data';
    warning('Be careful, X_data had to be flipped, X_smoothed will be formatted differently than X_data');
end

%% X_data initial size
X_length = size(X_data,2);

%% increase X_data size to avoid weird border effects
X_data_bis = [X_data(end:-1:1), X_data, X_data(end:-1:1)];

%% convolution parameter for smoothing
convprm = 1/smooth_kernel*ones(smooth_kernel,1);

%% perform the smoothing
X_smoothed_big = conv(X_data_bis, convprm, 'same');

%% keep only the relevant part
X_idx = (X_length + 1):(X_length*2);
X_smoothed = X_smoothed_big( X_idx );

end

