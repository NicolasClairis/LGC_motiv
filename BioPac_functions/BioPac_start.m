function [stim] = BioPac_start()
%[stim] = BioPac_start()
% BioPac_start will start BioPac acquisition system for grip acquisition
%
% OUTPUTS
% stim: structure with BioPac ref
% 

stim.u_out = udp('127.0.0.1', 2012, 'LocalPort', 15010);
fopen(stim.u_out);

end