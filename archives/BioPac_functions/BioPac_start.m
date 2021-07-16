function [u_out] = BioPac_start()
%[u_out] = BioPac_start()
% BioPac_start will start BioPac acquisition system for grip acquisition
%
% OUTPUTS
% u_out: BioPac reference
% 

u_out = udp('127.0.0.1', 2012, 'LocalPort', 15010);

% open BioPac channel
fopen(u_out);

end % function