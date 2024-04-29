%% script to show how mediation works
% Designed by N.Clairis - 2023
%
% See also mediation.m

%% Demo 1 with only indirect effect, no direct effect
n_samples = 20;
noise_range = 2;

% initial vector
X = 1:n_samples;

% mediator = transformation of X + some noise
noise_M = rand(1,n_samples).*noise_range;
path_a = 1.5;
M = 0.1 + path_a.*X + noise_M;

% build Y based on transformation of M
noise_Y = rand(1,n_samples).*noise_range;
path_b = 2.5;
Y = 0.3 + path_b.*M + noise_Y;

%% look at the resulting mediation analysis
[a_fit,b_fit,c_fit,c_prime_fit, pval] = mediation(X,M,Y);

disp('As you can see, the direct path from X to Y is significant, only when not taking M into account.');
disp('Including M shows that M is closer to the output Y than X');
disp('');
disp(['In principle, ',...
    'a_fit (=',num2str(round(a_fit,3)),') and ',...
    'b_fit (=',num2str(round(b_fit,3)),') ',...
    ' should be ~= to path_a (=',num2str(round(path_a,3)),') ',...
    'and path b (=',num2str(round(path_b,3)),')']);
disp(['while c'' (=',num2str(round(c_prime_fit,3)),') should be ~= 0.']);

%% Demo 2 with both indirect and direct paths significant
n_samples = 20;
noise_range = 2;

% initial vector
X = 1:n_samples;

% mediator = transformation of X + some noise
noise_M = rand(1,n_samples).*noise_range;
path_a = 1.5;
M = 0.1 + path_a.*X + noise_M;

% build Y based on transformation of M + X
noise_Y = rand(1,n_samples).*noise_range;
path_b = 2.5;
path_c_prime = 4.5;
Y = 0.3 + path_b.*M + path_c_prime.*X + noise_Y;

%% look at the resulting mediation analysis
[a_fit,b_fit,c_fit,c_prime_fit, pval] = mediation(X,M,Y);

disp('As you can see, the direct path from X to Y is significant, only when not taking M into account.');
disp('Including M shows that M is closer to the output Y than X');
disp('');
disp(['In principle, ',...
    'a_fit (=',num2str(round(a_fit,3)),'), ',...
    'b_fit (=',num2str(round(b_fit,3)),') ',...
    'c_prime_fit (=',num2str(round(c_prime_fit,3)),') ',...
    ' should be ~= to path_a (=',num2str(round(path_a,3)),'), ',...
    'path b (=',num2str(round(path_b,3)),') and ',...
    'path c'' (=',num2str(round(path_c_prime,3)),') and ']);
