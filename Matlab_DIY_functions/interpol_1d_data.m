function [dataInterpolated, interpolIndex] = interpol_1d_data(rawData, saturationThresholds, dispInterp)
%[dataInterpolated, interpolIndex] = interpol_1d_data(rawData, saturationThresholds)
% interpol_1d_data will interpolate the 1D data inside the rawData vector
% based on the lower and upper thresholds defined in saturationThresholds.
% The function will also apply interpolation on missing data if those data
% are indicated by NaN values.
%
% INPUTS
% rawData: 1D data containing samples to interpolate (either because empty
% data or because saturated data)
%
% saturationThresholds: [x1 x2] = [lower and higher] thresholds defining
% saturation in rawData values
%
% dispInterp: display interpolated data and original data (1) or no figure
% at all (0)
%
% OUTPUTS
% dataInterpolated: 1D data same size as rawData entered in inputs but with
% bad samples replaced by interpolated data
%
% interpolIndex: index of the samples that were replaced by interpolated
% data
%

%% extract index for all samples
nSamples = length(rawData);
x = 1:nSamples;

%% extract thresholds for saturation
lowThreshold = saturationThresholds(1);
highThreshold = saturationThresholds(2);

%% define samples where saturation happens
badSamples = false(1,nSamples);
iCluster = 0;
for iSample = x
    %% interpolate all missing data
    if isnan(rawData(iSample))
        badSamples(iSample) = true;
    end
    
    %% interpolate values for data reaching saturation (at lower or at higher threshold)
    if (rawData(iSample) == highThreshold) ||...
            (rawData(iSample) == lowThreshold)
        iCluster = iCluster + 1;
        if iCluster > 1
            badSamples(iSample) = true;
        end
    else
        if iCluster > 1
            badSamples(iSample - 1) = false;
        end
        iCluster = 0;
    end
end % loop through all samples
interpolIndex = find(badSamples);
% define good samples (ie where no interpolation is needed)
goodSamples = x;
% remove all bad samples
goodSamples(badSamples) = [];
% perform the interpolation
interpDataSamples = interp1(goodSamples,rawData(goodSamples),interpolIndex,'spline');
% replace bad samples by interpolated data
dataInterpolated = NaN(1,nSamples);
dataInterpolated(goodSamples) = rawData(goodSamples);
dataInterpolated(badSamples) = interpDataSamples;

%% display original vs interpolated  data
if dispInterp == 1
    figure;
    % display original data
    plot(x, rawData,'-','Color','k');
    % display interpolated data
    hold on;
    plot(x, dataInterpolated,'-','Color','r');
    
end % display data

end % function