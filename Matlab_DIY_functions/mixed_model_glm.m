% Script to perform a mixed-model GLM including both fixed and random
% (ie, subjects) effects based on chatgpt help.
% This script requires the Statistics and Machine Learning Matlab Toolbox
% which contains the fitlme.m function.
%
% Designed by N. Clairis - november 2023

% pool data inside a table
dataTable = table([subject_id';subject_id';subject_id';subject_id'],...
    [ones(NS,1); ones(NS,1)*2; ones(NS,1)*3; ones(NS,1)*4],...
    [Y1'; Y2'; Y3'; Y4'],... 
    [X';X';X';X';],...
    'VariableNames', {'SubjectID', 'TimePoint', 'Y', 'X'});

% Define the mixed-effects model formula
formula = 'Y ~ 1 + X + TimePoint + (1|SubjectID)';

% Fit the mixed-effects model using fitlme
mdl = fitlme(dataTable, formula);

% Display the global results
disp(mdl);

%% main results
% mdl.Coefficients.Estimate will provide the corresponding betas
% mdl.Coefficients.pValue will provide the corresponding p.values
% mdl.Coefficients.Name will provide the name for each variable in the
% model