function [metabolite_nm_bis] = metab_div_rnm(metabolite_nm)
%[metabolite_nm_bis] = metab_div_rnm(metabolite_nm)
% metab_div_rnm will rename the metabolites that are based on a ratio so
% that the name is properly displayed on the figures.
%
% For example, 'Glu_div_GSH' will become 'Glu/GSH'.
% INPUTS
% metabolite_nm: metabolite name
%
% OUTPUTS
% metabolite_nm_bis: metabolite name with / if based on ratio

metabolite_nm_bis = strrep(metabolite_nm,'_div_','/');
metabolite_nm_bis = strrep(metabolite_nm_bis,'_plus_','+');
end % function