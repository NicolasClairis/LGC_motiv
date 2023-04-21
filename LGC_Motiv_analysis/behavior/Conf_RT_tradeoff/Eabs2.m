function Eax = Eabs2(mu,V)
% approximates E[|x|] is x~N(mu,V)
t = mu./sqrt(V);
t(isnan(t)) = 0;
Eax = sqrt((2.*V)./pi).*exp(-0.5*t.^2) + mu.*(2*sig(pi*t./sqrt(3))-1);