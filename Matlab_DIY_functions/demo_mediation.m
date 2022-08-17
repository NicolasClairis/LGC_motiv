%% script to show how mediation works

% initial vector
X = 1:20;

% mediator = transformation of X + some noise
noise_M = rand(1,20);
M = 1.5.*X + noise_M;

% build Y based on transformation of M
noise_Y = rand(1,20);
Y = 2.5.*M + noise_Y;

%% look at the resulting mediation analysis
[a,b,c,c_prime, pval] = mediation(X,M,Y);

disp('As you can see, the direct path from X to Y is significant, only when not taking M into account.');
disp('Including M shows that M is closer to the output Y than X');