function [gx]=aDDM_model_analytic_to_inv(x,param,u,xx)

%% VARIABLES
FIXATIONS=u(3:end,:);
uu=u(1:2,:);

%% PARAM
bound   = 1; % .15;  %upper bound (correct answer)
drift   = param(1);% .001*DV(trial)/4; %drift rate (units/sec) 0.1
start   = param(2); %0.001*bias(trial)*2;
p.s     = param(3); %standard deviation of drift (units/sec)
ntheta1 = param(4);
ntheta2 = ntheta1;

ter2=0.2;


%% ========================= PARAMETER TRANSFORMATION

% %% BOUNDS BETWEEN 0 AND 0.5
% p.a=bound_parameter(bound,0,1); % pour prior p.a � 0.15 Prior bound = -1;
% start=p.a/2;
% p.start=start;
%
% %% START BETWEEN -p.a AND p.a
% % p.start=-p.a+(p.a+p.a)/(1+exp(-start)); % % pour prior p.start � 0.0 Prior bound = 0;
%
% %% THETA MUST BE BETWEEN 0 AND 1
% theta1=bound_parameter(ntheta1,0,1);
% theta2=bound_parameter(ntheta2,0,1);
%
% %% TER MUST BE POSITIVE
%
% ter=exp(ter2);

p.a     = bound;
p.start = p.a/2+start;
theta1  = ntheta1;
theta2  = ntheta2;
ter     = ter2;


%% ================= FIXATION INTEGRATION

FIXATIONS(isnan(FIXATIONS))=0;
DV=((uu(1,:)-theta1*uu(2,:)).*sum(FIXATIONS==1) + (-(uu(2,:)-theta2*uu(1,:))).*sum(FIXATIONS==-1))./((sum(FIXATIONS==-1)+sum(FIXATIONS==1)));

%% P STRUCTURE DEFINITION
    p.u=drift*DV+0.0000001;


    Pexp=1-(EZ2PE(p.u, p.start, p.a, p.s));
    
    % v is drift
    % x is starting point
    % a is boundary
    % s is random disturbance (standard deviation of noise)
    
    
    gx(1,:)=Pexp;
    if Pexp>=0.5 && Pexp<=1
        % gx(2,:)=ter+vRT(1);
        gx(2,:)=EZ2cmrt(p.u, p.start, p.a, ter, p.s) ;
    elseif Pexp<0.5 && Pexp>=0
        % gx(2,:)=ter+vRT(2);
        gx(2,:) = EZ2emrt(p.u, p.start, p.a, ter, p.s);
    elseif Pexp>1
        gx(1,:)=1;
    elseif Pexp<0
        gx(1,:)=0;
    else
        error('Pexp is anormal')
    end




