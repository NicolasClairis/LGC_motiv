%% script to compare mediation and SEM in R

%% working dir
SEM_path = fullfile('P:','boulot','postdoc_CarmenSandi','results','SEM');

%% general parameters
n_samples = 100;
noise_range = 100;

% initial vector
X0 = 1:n_samples;

%% Demo 1: no direct effect
% mediator 1 = transformation of X0 + some noise
noise_X1 = rand(1,n_samples).*noise_range;
path_X0_to_X1 = 1.5;
X1 = 2.3 + path_X0_to_X1.*X0 + noise_X1;

% principal mediator based on transformation of X1
noise_M = rand(1,n_samples).*noise_range;
path_X1_to_M = 3.2;
M = -4.5 + path_X1_to_M.*X1 + noise_M;

% build Y based on transformation of M
noise_Y = rand(1,n_samples).*noise_range;
path_M_to_Y = 3.2;
Y = -4.5 + path_M_to_Y.*M + noise_Y;

%% look at the resulting mediation analysis
[a,b,c,c_prime, pval] = mediation(X1,M,Y);

%% create matrix
mtrx_THE = [X0', X1', M', Y'];
% save data to load in R
save([SEM_path,filesep,'demo_SEM_1.mat'],'mtrx_THE');


%% Demo 2: include both indirect and direct effects on Y
% Direct effect = X1 to Y

% build Y based on transformation of M + X1
noise_Y = rand(1,n_samples).*noise_range;
path_M_to_Y = 3.2;
path_X1_to_Y = 5.4;
Y = -4.5 + path_M_to_Y.*M + path_X1_to_Y.*X1 + noise_Y;
% create matrix
mtrx_THE = [X0', X1', M', Y'];
% save data to load in R
save([SEM_path,filesep,'demo_SEM_2.mat'],'mtrx_THE');

%% Demo 3: include both indirect and direct effects on Y
% Direct effect = X1 to Y + X0 to Y

% build Y based on transformation of M + X1 + X0
noise_Y = rand(1,n_samples).*noise_range;
path_M_to_Y = 3.2;
path_X1_to_Y = 5.4;
path_X0_to_Y = 2.7;
Y = -4.5 + path_M_to_Y.*M + path_X1_to_Y.*X1 + path_X0_to_Y.*X0 + noise_Y;
% create matrix
mtrx_THE = [X0', X1', M', Y'];
% save data to load in R
save([SEM_path,filesep,'demo_SEM_3.mat'],'mtrx_THE');
