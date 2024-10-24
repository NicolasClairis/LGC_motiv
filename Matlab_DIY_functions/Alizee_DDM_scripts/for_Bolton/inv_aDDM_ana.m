% function param=inv_aDDM_ana(model_n,dispo,sub)
% %
clear
close all
clc

addpath(genpath('/Users/user/Documents/MATLAB_toolbox/VBA-toolbox-master'))


dispo=1;
sub=1;

verb=1;
load('data_example_for_inv.mat')

% Constraint paramters values:
constr=0;

if constr==1
    model_n=1;
else
    model_n=2;
end


% Mod�les � tester :

models={'aDDM_model_analytic_to_inv','aDDM_model_analytic_to_inv_constr'};

% bound =param(1); % .15;  %upper bound (correct answer)
% drift =param(2)/100;% .001*DV(trial)/4; %drift rate (units/sec) 0.1
% start=param(3); %0.001*bias(trial)*2;
% p.s =param(4); %.1;  %standard deviation of drift (units/sec)
% ter=param(5);
% ntheta1=param(6);
% % theta2=1-theta1;%param(7);
% ntheta2=param(7);

model_name=models{model_n};

% disp(model_n);

% bound=unbound_parameter(0.8,0,1); %0.23;
% drift=2;
% start=unbound_parameter(0.01,0,1);
% % var_drift=35;
% var_drift=20;
% ter=log(0.3);




if constr==1
    drift     = 0.01; % drift value depends on the range of the input DV (here, between -100 and 100 so small value of drift)
    start     = 0;
    dev_drift = 0.1;
    theta     = 0.8;
else
    drift=log(0.01);
    start=0;
    dev_drift=0.1;
    theta=unbound_parameter(0.5,0,1);
end
    var_PRIOR=1;

switch model_name
    case 'aDDM_model_analytic_to_inv'
        model_fun = @aDDM_model_analytic_to_inv;
    case 'aDDM_model_analytic_to_inv_constr'
        model_fun = @aDDM_model_analytic_to_inv_constr;
end
prior=[drift start dev_drift theta]';
var_prior=eye(length(prior))*var_PRIOR;

modeles={model_fun};
parameter=length(prior);

n_test=1:432;

row=0;

for subj =sub
    
    row=row+1;
    
    N=repmat(n_test,1,1);
    
    
    u=[V1'; V2'];
    
    
    Ntrials=max(N);
    
    length_fix=1050;
    
    u_r=NaN(length_fix,Ntrials);
    for mod=1:length(modeles)
        
        
        g_name=modeles{mod};
        nb_param=parameter(mod);
        
        u_r(1:2,:)=u;
        for t=1:432
            u_r(3:length_fix,t)=[fixations_V1{t}'; NaN(length_fix-2-length(fixations_V1{t}),1)];
        end
        
        
        f_name=[];
        
        %% Definition of model options
        
        dim = struct('n',0,...  % number of hidden states
            'p',2,... % total output dimension
            'n_theta',0,... % number of evolution parameters
            'n_phi', nb_param,... % number of observation parameters
            'n_t',Ntrials); % number of trials
        
        options.DisplayWin = dispo; % Display setting
        options.GnFigs = 0; % Plotting option
        options.isYout = zeros(2,Ntrials); % vector of the size of y, 1 if trial out
        options.verbose=verb;
        sources(1)=struct('out',1,'type',1);% 1 if binary data, 0 if continuous data
        sources(2)=struct('out',2,'type',0);% 1 if binary data, 0 if continuous data
        
        options.sources=sources;
        
        
        %% Definition of priors
        
        % Priors on parameters (mean and Covariance matrix)
        
        % Observation parameters :
        priors.muPhi = prior;
        priors.SigmaPhi = var_prior;
        
        options.priors = priors;
        
        
        y=[choice' ; RT'];
        
        
        [options.priors.a_sigma(2), options.priors.b_sigma(2)] = getHyperpriors(y(2,:),0.1,0.9);
        [options.priors.a_sigma(1), options.priors.b_sigma(1)] = getHyperpriors(y(1,:),0.5,1);
        
        
        %% Performing the inversion
        options.figName='RT : VB-Laplace inversion of static model';
        [posteriorr,outr] = VBA_NLStateSpaceModel(y,u_r,f_name,g_name,dim,options);
        
    end
    
end
