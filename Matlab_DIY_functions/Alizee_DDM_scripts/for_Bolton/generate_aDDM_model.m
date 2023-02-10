clear
close all

load('data_example.mat')

% u_r is a 2 cells variable.
% First cell is a 2*432 matrix. Each column is a trial, first line is the
% subjective value of the left option, second line is the subjective value
% of the right option.
% Second cell comprises fixations data. It's a 1050*432 matrix. Each column
% is a trial; Each line is a fixation point (16ms per line). If 1, fixation
% is for the left option, if -1 for the right option. 0 is for somewhere on
% the screen but not on any option. NaN is for missing data.

%% Please put a debug stop on line 118 in aDDM_model_dynamic to visualize each trial dynamics.


% aDDM parameters
drift     = 0.01; % drift value depends on the range of the input DV (here, between -100 and 100 so small value of drift)
start     = 0;
dev_drift = 0.1;
theta     = 0.8;

param     = [drift start dev_drift theta];

% call aDDM dynamic model
gx_dyn(:,:)    =  aDDM_model_dynamic([],param,u_r,[],0);


% call aDDM analytic model
gx_ana(:,:)    =  aDDM_model_analytic([],param,u_r,[]);


figure

subplot(2,2,1)
scatter(u_r{1}(1,:)-u_r{1}(2,:),gx_dyn(1,:),'k','filled')
title('Simulation on one sub data with dynamic aDDM')
xlabel('DV')
ylabel('Proportion of choice for option 1')
subplot(2,2,2)
scatter(abs(u_r{1}(1,:)-u_r{1}(2,:)),gx_dyn(2,:),'k','filled')
xlabel('|DV|')
ylabel('RT')
ylim([0 4])

subplot(2,2,3)
scatter(u_r{1}(1,:)-u_r{1}(2,:),gx_ana(1,:),'k','filled')
title('Simulation on one sub data with analytic aDDM')
xlabel('DV')
ylabel('Proportion of choice for option 1')
subplot(2,2,4)
scatter(abs(u_r{1}(1,:)-u_r{1}(2,:)),gx_ana(2,:),'k','filled')
ylim([0 4])
xlabel('|DV|')
ylabel('RT')