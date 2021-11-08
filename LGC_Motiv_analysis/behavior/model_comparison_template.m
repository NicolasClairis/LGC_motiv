
nModels = 2; % nber of models to test
NS = 50; % nber of subjects

free_energy_mtrx = NaN(nModels, NS);

%% run the models you want to run
for iModel = 1:nModels
    %% run the model for each individual
    for iS = 1:NS
        % inversion
        [posterior,out] = VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options);
        free_energy_mtrx(iMdl, iS) = out.F;
    end % loop over subjects
end % loop over models

%% compare models
options.DisplayWin = 1;
options.verbose = 1;
options.modelNames = {'model1_linearF','model2_FwithRecovery'};

% run the model comparison
[posterior, out] = VBA_groupBMC(free_energy_mtrx, options);

% VBA_groupBMC will display the exceedance model probability for each model
% and show you which one is the best