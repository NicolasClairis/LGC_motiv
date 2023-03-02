function [w,wP] = get_wMDP0(R,alpha,beta,gamma,kappa,pow,v0,T)
% Derive MDP-optimal stopping rule for MCD
% function [w,dms] = get_wMDP0(R,alpha,beta,gamma,kappa,pow,v0,T)
% This function derives the optimal threshold w(t) for the discounted
% benefit of effort investment, given parameters that control decision
% confidence and effort cost.
% IN:
%   - R: decision importance
%   - alpha: unitary effort cost
%   - beta: type#1 effort efficiency
%   - gamme: type #2 effort efficiency
%   - kappa: effort intensity (per unit of time)
%   - pow: power of effort cost function
%   - v0: initial value variance
%   - T: temporal horizon
% OUT:
%   - w: 1xT vector of discounted benefit thresholds
%   - wP: 1xT vector of confidence thresholds

% transition matrix: P(dmu(t)|dmu(t-1))
maxdmu = 4*sqrt(2*gamma*kappa*T); % bound = 4 standard deviations
ddmu = maxdmu*1e-3;
dmut = vec(-maxdmu:ddmu:maxdmu);
i0 = find(dmut==0);
Pt = NaN(length(dmut),length(dmut));
for i=1:length(dmut)
    Pt(i,:) = exp(-0.5*(dmut'-dmut(i)).^2./(4*gamma*kappa));
    Pt(i,:) = Pt(i,:)./sum(Pt(i,:));
end

% initialize oMCD variables
w = NaN(1,T); % discounted benefit threshold
EQ = NaN(length(dmut),T); % expected optimal discounted benefit

% At t=T
w(T) = R/2 - alpha*(kappa*T)^pow; % minimum discounted benefit
vt = 1/((1/v0)+beta*kappa*(T-1)); % value variance at T-1
EQ(:,T) = R*getPc_dmu(kappa,beta,vt,dmut,gamma,0) - alpha*(kappa*T)^pow;

% discounted benefits' threshold: backward induction
for t=1:T-1 % loop through decision times (backward)
    vt = 1/((1/v0)+beta*kappa*(T-t)); % value variance at T-t
    Q0t = R*getPc_dmu(0,0,vt,dmut,0,0) - alpha*(kappa*(T-t))^pow;
    EQ(:,T-t) = Pt*max([Q0t,EQ(:,T-t+1)],[],2); % Bellman recurrence equation
    iu = find(Q0t(i0:end)>=EQ(i0:end,T-t+1),1); % crossing index
    if ~isempty(iu) % if the crossing happens whithin the dmu grid
        w(T-t) = Q0t(i0+iu-1); % optimal threshold
    else % the crossing must happen outside the dmu grid
        w(T-t) = R - alpha*kappa*(T-t); % set to max[discounted benefit]
    end
end

% confidence threshold
wP = (w+alpha*(kappa*(1:T)).^pow)./R;

