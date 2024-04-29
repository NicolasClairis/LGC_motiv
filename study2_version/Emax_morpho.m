function [predictedForce] = Emax_morpho(pli_a,pli_p,circ,length)
%% Estimate maximal theoretical force based on muscular volume of the forearm (ie. morphometric measurements)
% % INPUTS:
% - pli_a: pli cutané antérieur de l'avant-bras (mm)
% - pli_p: pli cutané postérieur de l'avant-bras (mm)
% - circ: circonférence de l'avant-bras (cm)
% - length: longueur de l'avant-bras (cm)
% % OUTPUT:
%     - predictedForce: maximal theoretical force (Newton)
 
% physiological cross sectional area: arm section area minus fat section area minus bone section area
% bone section from litterature, eg Hsu et al, 1993 J.Biomechanics
pli = pli_a + pli_p; 
pcsa = pi*(circ/(2*pi) - (pli)/40)^2 - (.82 + .98);
 
% force is proportional to pcsa + correction form arm length
% constants from Neu et al. 2001, Am. J Physiol Endocrinol Metab
 
predictedForce = pcsa * (2.45 + .288*length);

end