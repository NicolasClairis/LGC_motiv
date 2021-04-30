function [gx] = LGCM_g_observation(x,P,u,inG)
% [gx] = LGCM_g_observation(x,P,u,inG)
% LGCM_g_observation= observation function for 
%
% INPUTS
% x: ?
%
% P: parameters to estimate
%
% u: task parameters
%
% inG: structure with additional variables of interest
%
% OUTPUTS
% gx: output of the observation function

%% load parameters to estimate
iP = 1;
kR = P(iP);
iP = iP + 1;
kE = P(iP);
iP = iP + 1;
kF = P(iP);

% note: possible to pool mental and physical effort tasks together to
% estimate one single kR parameter for both tasks while estimating kE and
% kF separately for each

%% load task parameters
trialN = u(1);
Rleft = u(2);
Rright = u(3);
Eleft = u(4);
Eright = u(5);
sumEprev = u(6);
intEprev = u(7);

%% compute variables of interest
DeltaR = Rleft - Rright;
DeltaE = Eleft - Eright;
% define fatigue component to consider
switch inG.fatigue
    case 'trialN' % use trial number as as simple account of the cumulated level of fatigue
        FatigueComponent = trialN;
    case 'sumEffort' % use the sum of all the efforts done until now
        FatigueComponent = sumEprev;
    case 'integralEffort' % integral of the effort performed instead of sum of effort levels
        % mental effort: sum of the number of correct answers provided
        % physical effort: integral of the efforts performed until the
        % current trial
        FatigueComponent = intEprev;
end

%% output
% gx = 1./( 1 + exp(- ( (1 + kR.*DeltaR) - (1 + kE.*DeltaE).*(1+kF.*trialN) ) ));
gx = 1./( 1 + exp(- ( (kR.*DeltaR) -kE.*DeltaE.*(1+kF.*FatigueComponent) ) ));

end % function