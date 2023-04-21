function [Pc,Eadm,Vadm] = getPc_dmu(z,beta,v0,dm0,gamma,bias)
% approximates expected confidence level
v = 1./((1./mean(v0)) + beta.*z); % V[dv]
% dm0 = -diff(m0); % initial dm
Eadm = Eabs2(dm0,2*gamma*z); % E[|dm|]
Vadm = abs(dm0).^2 + 2*gamma.*z - Eadm.^2; % V[|dm|]
a = 1/2;%0.5066; % !!! set to 0 to remove the variance term !!!
b = 3/4;%0.7527; % !!! This was neglected in Lee&Daunizeau (2020) !!!
lambda = pi./sqrt(6*v);
Pc = sig(bias+(lambda.*Eadm./sqrt(1+a*((lambda.^2).*Vadm).^b)));