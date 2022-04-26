% function to simulate data and interpolate it with interp1 function from
% matlab.
%
% Developed by Nicolas Clairis with Catherine Bratschi - March 2022


% define threshold for saturation
highThreshold = 10;
lowThreshold = 1;

% create oscillating variable
d = repmat([1:10,10:-1:1],1,10);
% add saturation on top
d(11:30) = 10;
% add saturation at bottom
d(41:60) = 1;

% extract x index
x = 1:length(d);

% define samples where saturation happens
badSamples = false(1,length(x));
iCluster = 0;
for iSample = x
    if (d(iSample) == highThreshold) || (d(iSample) == lowThreshold)
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
badSamplesIdx = find(badSamples);
% define good samples (ie where no interpolation is needed)
goodSamples = x;
% remove all bad samples
goodSamples(badSamples) = [];
% perform the interpolation
d3 = interp1(goodSamples,d(goodSamples),badSamplesIdx,'spline');
% look at the data
plot(goodSamples,d(goodSamples),'o',badSamplesIdx,d3,'*')
% replace bad samples by interpolated data
d_corrected = NaN(size(d));
d_corrected(goodSamples) = d(goodSamples);
d_corrected(badSamples) = d3;
figure;
figure;plot(x,d,'-','Color','k');hold on;
plot(x,d_corrected,'-','Color','r');