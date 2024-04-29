function [predictedForce] = Emax_morpho(pli_a,pli_p,circ,length)
% [predictedForce] = Emax_morpho(pli_a,pli_p,circ,length)
%% Estimate maximal theoretical force based on muscular volume of the forearm (ie. morphometric measurements)
%
% INPUTS:
% pli_a: pli cutané antérieur de l'avant-bras (mm)
%
% pli_p: pli cutané postérieur de l'avant-bras (mm)
%
% circ: circonférence de l'avant-bras (cm)
%
% length: longueur de l'avant-bras (cm)
%
% OUTPUT:
% predictedForce: maximal theoretical force (Newton)
 
%% physiological cross sectional area (PCSA):
% sum of inner and outer fat skin fold = fat section area
pli = pli_a + pli_p;

% arm section area minus fat section area
pcsa_a = pi*((circ/(2*pi)) - (pli/40))^2;

% Remove 12.2% of bone + muscle PCSA = approximate size of bone and marrow 
% based on (Forbes et al, 1988, Am.J.Clin.Nutr). Needs to be removed to
% focus on the muscle
pcsa = pcsa_a - pcsa_a*12.2/100;
 
% force is proportional to pcsa + correction for arm length
predictedForce = pcsa*(2.45 + 0.288*length); % Lionel Rigoux formula based on (Neu et al, 2001)
% predictedForce = pcsa*7.52; % Lionel Rigoux formula based on participants for Paris ICM grip
% predictedForce = pcsa.*2.385; % value based on LGC study 1 average

end