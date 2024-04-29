% script to look at SPM design matrix to understand where the bug comes
% from

%% load data
load('SPM.mat');

%% display data
% figure;imagesc(SPM.xX.X);

% this will display the data loaded into SPM => try to see if one session
% is missing (in that case it means there must be a bug somewhere) or
% whether two regressors are (almost) completely identical, meaning that
% they are too highly correlated and it makes SPM crash.
figure;imagesc(SPM.xX.nKX);