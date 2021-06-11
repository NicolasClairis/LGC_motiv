% script for simulating first model. where f = [] and everything is in observation

%% clean workspace before startingclearvars;
clearvars;
close all;
clc;

%% variables initialization

% place to store results
results_folder = 'D:\Matlab codes\Data simulation\Simulated Data\';

% define main parameters: number of trials, number of conditions and
% parameters to simulate
n_sessions = 4;
n_conditionsPerSession = 2; % reward AND punishment
n_trialsPerCondition = 22;
n_trialsPerSession = n_trialsPerCondition*n_conditionsPerSession;
n_trials = n_trialsPerSession*n_sessions;
% number of conditions
n_R_levels = 3;
n_E_levels = 3;
punishment_yn = 'yes';
[R_money] = R_amounts(n_R_levels, punishment_yn);
% define parameters to simulate
kFp = (1/4)*(1/43); % define fatigue as 25% effect at the last trial
kFm = (1/10)*(1/43); % define fatigue mental as 10% effect at the last trial
kR = 0.6; % predict a high sensitivity
kP = 0.6; % predict a high sensitivity
kEp = 0.3; % compensate for the levels of effort which is 1-2-3 while kR/kP is 0.3
kEm = 0.3; % compensate for the levels of effort which is 1-2-3 while kR/kP is 0.3

% define functions to use. in this case we do not use evaluation, only observation
f_fname = [];
g_fname = @g_observation;

% define parameters to simulate
theta = [];
n_F_prm = length(theta);
phi = [kR kP kEp kEm kFp kFm];
% phi = [kR kP kEp kEm];
n_G_prm = length(phi);
alpha = Inf;
sigma = Inf;
% in case we have f_evaluation, it would be our fatigue param but not here.
% (set initial physical and mental fatigue at zero (but could be higher depending on the baseline))
x0 = [];
n_hiddenStates = length(x0);

%% compute choice_matrices, prepare the variables 

% enter design matrix, compute all 4 blocks and then merge them into one combined matrix
var = NaN(5, n_trials/4,n_sessions);
for iSession = 1:n_sessions
    [choice_opt] = choice_option_design(n_R_levels, n_E_levels, punishment_yn, n_trials/4, R_money);
    toBeSaved.choice_opt.(['Session',num2str(iSession)]) = choice_opt;
    deltaR = choice_opt.R.left - choice_opt.R.right;
    deltaE = choice_opt.E.left - choice_opt.E.right;
    R_or_P = choice_opt.R_or_P;
    RP_trials = NaN(1,n_trials/4);
    for iTrial = 1:n_trials/4
        switch R_or_P{iTrial}
            case 'P'
                RP_trials(iTrial) = 0;
            case 'R'
                RP_trials(iTrial) = 1;
        end
    end % trial loop
    
    %  alternate between physical and mental effort
    if ismember(iSession,[1,3])
        Ep_or_Em_trials = ones(1,n_trialsPerSession);
    elseif ismember(iSession,[2,4])
        Ep_or_Em_trials = zeros(1,n_trialsPerSession);
    end
    
    % prepare the matrix of our variables
    var(:,:,iSession) = [deltaR; deltaE; RP_trials; Ep_or_Em_trials;(1:44)-1];
end % loop through blocks

% reshape it as one huge trial
var = reshape(permute(var,[2 3 1]),size(var,2)*4,[])';

options = struct;
options.sources.type = 1; % binary data
options.DisplayWin  = 1; % display figure during inversion
options.verbose     = 1; % display text during inversion (1) or not (0)
% define number of parameters to estimate
options.dim = struct('n', n_hiddenStates,... % number of hidden states
    'n_t', n_trials,... % number of trials across runs
    'n_theta',n_F_prm,... % number of evolution parameters
    'n_phi',n_G_prm); % number of observation parameters


% % use multisession to pool 3 runs together
% options.multisession.split = repmat(n_trialsPerSession, 1, n_sessions);
% % fix the parameters (only the hidden states should vary in each
% % session, in theory)
% options.multisession.fixed.theta = 1:n_F_prm; % fix the evolution parameters across sessions (learning rate, R/P weights)
% options.multisession.fixed.phi = 1:n_G_prm; % fix the observation parameters across sessions (temperature and side bias)

%% Launch VBA simulations
[choiceLeftOptionSimu, x, x0, eta, e, variable] = VBA_simulate(n_trials, f_fname, g_fname,...
    theta, phi, var, alpha, sigma, options, x0);
% figure();
% subplot(1,2,1)
% var(1),choiceLeftOptionSimu(var(1,:) == -2)
% subplot(1,2,2)
% bar(var(2,:), choiceLeftOptionSimu)
lol(1) = sum(choiceLeftOptionSimu(var(1,:) == -2))
lol(2) = sum(choiceLeftOptionSimu(var(1,:) == -1))
lol(3) = sum(choiceLeftOptionSimu(var(1,:) == 1))
lol(4)= sum(choiceLeftOptionSimu(var(1,:) == 2))
histogram(lol)
%% Launch posterior computation
% for posterior calculation, need priors on our k sensitivities. define them with a mean of 0 and sd of 1
options.priors.muPhi = zeros(length(phi), 1)*0.5; %(moyenne des priors sur les paramètres)
options.priors.SigmaPhi = eye(length(phi))/0.5; % (matrice de covariance des paramètres)

[posterior,out] = VBA_NLStateSpaceModel(choiceLeftOptionSimu, var, f_fname, g_fname, options.dim, options);

%% save the data
toBeSaved.n_trials = n_trials;
toBeSaved.theta = theta;
toBeSaved.phi = phi;
toBeSaved.u = var;
toBeSaved.simulated_data_left_opt = choiceLeftOptionSimu;
toBeSaved.posterior = posterior;
toBeSaved.out = out;

save([results_folder,'simulation_data.mat'],'toBeSaved');