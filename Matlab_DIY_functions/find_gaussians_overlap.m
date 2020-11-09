function [ X1, X2 ] = find_gaussians_overlap(muA, sigmaA, muB, sigmaB )
%[ X1, X2 ] = find_gaussians_overlap(muA, sigmaA, muB, sigmaB )
% find_gaussians_overlap will give you the one (or two) solutions for which
% there can be an overlap between two gaussian curves A(muA, sigmaA) and B(muB, sigmaB) following a normal distribution.
%
% INPUTS
% muA: mean of distribution A
%
% sigmaA: variance of distribution A
%
% muB: mean of distribution B
%
% sigmaB: variance of distribution B
%
% OUTPUTS
% X1, X2: 2 solutions to see for which values of X, the two gaussians
% overlap
% Note: X1 should be smaller than X2
%
% Note: when sigmaA=sigmaB, only one solution, so that output X1=X2
% 
% Developped by Nicolas Clairis with the help of Quentin Feltgen - march
% 2020

if sigmaA ~= sigmaB
%% define three big parameters of the equation
a = (1/(sigmaA^2)) - (1/(sigmaB^2));
b = -2*( (muA/(sigmaA^2)) - (muB/(sigmaB^2)) );
c = (muA/sigmaA)^2 - (muB/sigmaB)^2 + log( (sigmaA^2)/(sigmaB^2) );

%% extract 2 possible solutions
X1 = (-b - sqrt(b^2 - 4*a*c))/(2*a + eps); % a can be negative in that case, X1 will be higher than X2 => that's why reordered afterwards
X2 = (-b + sqrt(b^2 - 4*a*c))/(2*a + eps);

%% order them in ascending order
[Xorder] = sort([X1, X2]);
X1 = Xorder(1);
X2 = Xorder(2);

elseif sigmaA == sigmaB
    X1 = (muA^2 - muB^2)/(2*(muA - muB));
    X2 = X1;
end

end % function