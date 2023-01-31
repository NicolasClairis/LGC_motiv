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
% arm section area minus fat section area minus bone section area
% bone section from litterature (Hsu et al, 1993 J.Biomechanics)
pli = pli_a + pli_p;
pcsa_a = pi*(circ/(2*pi) - (pli)/40)^2;
% Remove 12.2% of the total PCSA = approximate size of bone and marrow 
% based on (Forbes et al, 1988, Am.J.Clin.Nutr). Needs to be ignored to
% focus on the muscle
pcsa = pcsa_a - pcsa_a*12.2/100;
 
% force is proportional to pcsa + correction form arm length
% constants from Neu et al. 2001, Am. J Physiol Endocrinol Metab
 
predictedForce = pcsa * (2.45 + .288*length);
% predictedForce = pcsa*7;
end