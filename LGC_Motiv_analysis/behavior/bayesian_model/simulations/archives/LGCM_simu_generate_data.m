function[choice_left] = LGCM_simu_generate_data(theta, phi, u_t, n_trials, mdlPrm)
%[choice_left] = LGCM_simu_generate_data(theta, phi, u_t, n_trials, mdlPrm)
%LGCM_simu_generate_data will generate data for the behavioral tasks provided the parameters you will define.
%
% INPUTS
% theta: information about evolution parameters
%
% phi: information about observation parameters
%
% u_t: matrix with task parameters ordered to fit with LGCM_g_observation
% and LGCM_f_evolution functions.
%
% n_trials: total number of trials to simulate
%
% mdlPrm: structure with information about the model to use
%
% OUTPUTS
% choice_left: binary variable equal to (1) when left option selected and
% to (0) when right option selected

%% extract model parameters
fatigueEvolOrObs    = mdlPrm.fatigueEvolOrObs;
fatigueType         = mdlPrm.fatigueType;
kRP_tasksSplitOrPool = mdlPrm.kRP_tasksSplitOrPool;

%% define functions used for the model

% evolution function
switch fatigueEvolOrObs
    case 'hidden state'
        f_fname = @LGCM_f_evolution;
        
        % initial fatigue state = null a priori
        switch kRP_tasksSplitOrPool
            case 'split' % fatigue
                x0 = 0;
                n_F_prm = 1;
                n_hiddenStates = 1;
            case 'pool' % physical and mental effort fatigue
                x0 = [0 0];
                n_F_prm = 2;
                n_hiddenStates = 2;
        end
    case 'observation parameter' % fatigue in G function
        f_fname = []; % no evolution function if no hidden state considered
        x0 = [];
        n_hiddenStates = 0;
        n_F_prm = 0;
end

% observation function
g_fname = @LGCM_g_observation;
switch fatigueEvolOrObs
    case 'hidden state'
        switch kRP_tasksSplitOrPool
            case 'pool'
                G_Phi_prm_nm = {'kR','kP','kEm','kEp'};
            case 'split'
                G_Phi_prm_nm = {'kR','kP','kE'};
        end
    case 'observation parameter'
        switch kRP_tasksSplitOrPool
            case 'pool'
                G_Phi_prm_nm = {'kR','kP','kEm','kEp','kFEm','kFEp'};
            case 'split'
                G_Phi_prm_nm = {'kR','kP','kE','kF'};
        end
end
n_G_prm = length(G_Phi_prm_nm);


% define number of parameters to estimate
dim = struct('n', n_hiddenStates,... % number of hidden states
    'n_t', n_trials,... % number of trials across runs
    'n_theta',n_F_prm,... % number of evolution parameters
    'n_phi',n_G_prm); % number of observation parameters

%% define parameters for simulations


options = struct;
options.sources.type = 1; % binary data
% options.GnFigs      = 0;
options.verbose     = 1; % display text during inversion (1) or not (0)

% store model parameters in the inputs of each function
if strcmp(fatigueEvolOrObs,'hidden state')
    options.inF.mdlPrm = mdlPrm;
end
options.inG.mdlPrm = mdlPrm;

options.dim = dim;

%% simulate the data
% [y,x,x0,eta,e,u] = VBA_simulate(n_t, f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0,fb);

% fixed the two following parameters according to Antonius but no idea why
% they have this value (Inf), would be nice to check their roles
alpha = Inf; % precision of the stochastic innovations
sigma = Inf; % precision of the measurement error
choice_left = VBA_simulate(n_trials, f_fname, g_fname,...
    theta, phi, u_t, alpha, sigma, options, x0);

end % function