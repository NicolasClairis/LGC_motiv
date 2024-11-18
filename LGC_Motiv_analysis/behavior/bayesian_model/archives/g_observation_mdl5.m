function [gx] = g_observation5( x, phi, var, inG )
% [gx] = g_observation5( x, phi, var, inG )
% Currently the following parameters are described as 
% phi is our parameters of interest, our sensitivities = [kR kP kEp kEm kFp kFm];
% var is our variables = [deltaR; deltaE; RP_trials; Ep_or_Em_trials;T-1;efforts_performed];
% x is empty.

%% load parameters
% reward sensitivity
kR = log(1+exp(phi(1)));
%punishment sensitivity
kP = log(1+exp(phi(2)));
% Physical effort sensitivity
kEp = log(1+exp(phi(3)));
% Mental effort sensitivity
kEm = log(1+exp(phi(4)));
%bias for physical task
bias = phi(5); % no positivity constraint on the bias
%physical fatigue
kFp = log(1+exp(phi(6)));
% mental learning
kLm = log(1+exp(phi(7)));

%% load variables
% Difference of incentive between the two choices
deltaR = var(1);
deltaP = var(2);
% Difference of effort between the two choices
deltaE = var(3);
% Are we in a physical or mental block
Ep_or_Em_trials = var(4);
% linear increasing value for fatigue
trial_T = var(5);
% time effect: physical fatigue for physical and current/previous
% efficiency for mental
time_effect = var(9);

%% Subjective value computation
SV = (kR.*deltaRP.*RP_trials + kP.*deltaRP.*(1-RP_trials))...
    - (deltaE.*(kEp + kFp*time_effect).*(Ep_or_Em_trials) +...
    deltaE.*(kEm - kLm*time_effect).*(1-Ep_or_Em_trials));

%% sigmoid transformation
gx = 1/(1+exp(-(SV+bias)));

end
