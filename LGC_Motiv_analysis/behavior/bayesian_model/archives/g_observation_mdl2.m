function [gx] = g_observation1( x, phi, var, inG )
% Currently the following parameters are described as 
% phi is our parameters of interest, our sensitivities = [kR kP kEp kEm kFp kFm];
% var is our variables = [deltaR; deltaE; RP_trials; Ep_or_Em_trials;T-1];
% x is not used in our case as it would be our sensitivities from the last trial
% InG I don't know is function

% MODEL 1, ONLY KRKP AGAINST KEP KEM DEFINED POSITIVE

kR = phi(1);
kP = phi(2);
kEp = phi(3);
kEm = phi(4);

deltaRP = var(1);
deltaE = var(2);
RP_trials = var(3);
Ep_or_Em_trials = var(4);
trial_T = var(5);

% compute equation of our model without evaluation, only observation
SV = (log(1+exp(kR))*(deltaRP)*RP_trials + log(1+exp(kP))*(deltaRP)*(1-RP_trials) - (deltaE)*log(1+exp(kEp)) *(Ep_or_Em_trials) - (deltaE)*log(1+exp(kEm))*(1-Ep_or_Em_trials));
gx = 1/(1+exp(-(SV)));

end
