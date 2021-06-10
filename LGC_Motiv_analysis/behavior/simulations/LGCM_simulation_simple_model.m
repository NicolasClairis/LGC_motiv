

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
kF1 = ;
kF2 = ;
kR = ;
kP = ;
kEp = ;
kEm = ;

% define functions to use
f_fname = simpleModel_f_evolution;
g_fname = simpleModel_g_evolution;

% define parameters to simulate
theta = [kF1 kF2];
n_F_prm = length(theta);
phi = [kR kP kEp kEm];
n_G_prm = length(phi);
alpha = Inf;
sigma = Inf;
% set initial physical and mental fatigue at zero (but could be higher depending on the baseline)
x0 = [0 0];
n_hiddenStates = length(x0);

% enter design matrix
u = NaN(4, n_trials);
for iSession = 1:n_sessions
    [choice_opt] = choice_option_design(n_R_levels, n_E_levels, punishment_yn, n_trials, R_money);
    deltaR = choice_opt.R.left - choice_opt.R.right;
    deltaE = choice_opt.E.left - choice_opt.E_right;
    R_or_P = choice_opt.R_or_P;
    RP_trials = NaN(1,n_trials);
    for iTrial = 1:n_trials
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
    
    sessionTrial_idx = (1:n_trials) + n_trials*(iSession - 1);
    u(:,sessionTrial_idx) = [deltaR, deltaE, RP_trials, Ep_or_Em_trials];
end % loop through blocks

options = struct;
options.sources.type = 1; % binary data
options.DisplayWin  = 1; % display figure during inversion
options.verbose     = 1; % display text during inversion (1) or not (0)
% define number of parameters to estimate
options.dim = struct('n', n_hiddenStates,... % number of hidden states
    'n_t', n_trials,... % number of trials across runs
    'n_theta',n_F_prm,... % number of evolution parameters
    'n_phi',n_G_prm); % number of observation parameters

% use multisession to pool 3 runs together
options.multisession.split = repmat(n_trialsPerSession, 1, n_sessions);
% fix the parameters (only the hidden states should vary in each
% session, in theory)
options.multisession.fixed.theta = 1:n_F_prm; % fix the evolution parameters across sessions (learning rate, R/P weights)
options.multisession.fixed.phi = 1:n_G_prm; % fix the observation parameters across sessions (temperature and side bias)

[choiceLeftOptionSimu, x, x0, eta, e, u] = VBA_simulate(n_t, f_fname, g_fname,...
    theta, phi, u, alpha, sigma, options, x0);