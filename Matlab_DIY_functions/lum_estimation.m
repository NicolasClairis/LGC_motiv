function [luminance] = lum_estimation(window)
%[luminance] = lum_estimation(window)
% lum_estimation extracts the luminance of the stimulus currently displayed
% on screen by Psychtoolbox (PTB).
%
% INPUTS
% window: PTB window where stimulus is displayed
%
% OUTPUTS
% luminance: estimated based on the red/green/blue extraction formula
%
% Initially written by Nicolas Clairis - 12/03/18 based on a formula found
% on internet

rgb_img = Screen(window,'GetImage',[]); % extract the stimulus as a specific RGB composition

% extract each colour layer
% (normalize to have values ranking from 0 to 1)
R = rgb_img(:,:,1)./255 ; % red layer
G = rgb_img(:,:,2)./255 ; % green layer
B = rgb_img(:,:,3)./255 ; % blue layer

% make an estimate of the luminance based on this RGB composition
luminance = mean(mean(0.299*R + 0.587*G + 0.114*B)); % 0 to 1 (0=black; 1= white)

end

