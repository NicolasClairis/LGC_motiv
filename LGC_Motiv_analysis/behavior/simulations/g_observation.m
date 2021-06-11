function [gx] = g_observation( x, phi, var, inG )
% Currently the following parameters are described as 
% phi is our parameters of interest, our sensitivities = [kR kP kEp kEm kFp kFm];
% var is our variables = [deltaR; deltaE; RP_trials; Ep_or_Em_trials;T-1];
% x is not used in our case as it would be our sensitivities from the last trial
% InG I don't know is function

kR = phi(1);
kP = phi(2);
kEp = phi(3);
kEm = phi(4);
kFp = phi(5);
kFm = phi(6);
deltaRP = var(1);
deltaE = var(2);
RP_trials = var(3);
Ep_or_Em_trials = var(4);
trial_T = var(5);

% compute equation of our model without evaluation, only observation
gx = sigmoid(kR*deltaRP*RP_trials + kP*deltaRP*(1-RP_trials) - kEp*deltaE*(1 + kFp*trial_T)*(Ep_or_Em_trials) - kEm*deltaE*(1 + kFm*trial_T)*(1-Ep_or_Em_trials));
% gx = sigmoid(kR*deltaRP*RP_trials + kP*deltaRP*(1-RP_trials) - kEp*deltaE*(Ep_or_Em_trials) - kEm*deltaE*(1-Ep_or_Em_trials));
end

