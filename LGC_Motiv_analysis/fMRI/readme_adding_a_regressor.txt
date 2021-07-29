***If you want to add a regressor:

1) go to which_GLM.m
a) Depending on the moment of the task when you want to add it,
you should add it accordingly to the parameters set to 0 by default.

For example, if during the moment of the display of the choice options,
add GLMprm.choice.(EpEm_nm).(RP_nm).New_Regressor = 0; to the list
(this will mean that by default the regressor is not included in the model)

The same logic applies for the other phases of a trial,
ie for the moment of display of the chosen option, write
GLMprm.chosen.(EpEm_nm).(RP_nm).New_regressor = 0;

for the moment of effort performance write
GLMprm.Eperf.(EpEm_nm).(RP_nm).New_regressor = 0;

for the moment of feedback write
GLMprm.fbk.(EpEm_nm).(RP_nm).New_regressor = 0;


2) you should open GLM_details.m and define the order when your regressor will be added
you can also state whether you want different forms of your regressor according to the
value of your parameter.

For example, let's say the regressor you add is the "reward chosen" value. Then you might
want to either enter it as a level of reward (1/2/3) or as a monetary amount (1.60/1.75/1.80)
or even by pooling both reward and punishment levels (-3/-2/-1/1/2/3) or using the absolute value
(1/2/3) across reward and punishment trials.
Within GLM_details.m you need to define when to place your regressor and also what form of it you
want to use

For a given value of the parameter, you should enter the following lines for your regressor of interest:

n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
% will update the regressor index and total number of regressors included in the model

reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG time-period-studied - new-regressor-name';
% ONSET or REG to indicate whether this is an onset or a parametric modulator
% time-period-studied to indicate which time period of the task is being modulated by the parametric modulator you are entering
% new-regressor-name: name of the new regressor you are entering. Be careful to not reuse a name already in use

disp([num2str(n_regs.(task_id_nm)),') time-period-studied- new-regressor-name ']);
% this is completely useless for the analysis, but it's important to help us to know what is contained within each model
% the idea is that launching this script already gives a complete description of the GLM
% time-period-studied: same as previous paragraph
% new-regressor_name: same as previous paragraph

n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
% will update the total number of regressors to the fact that temporal and/or spatial derivative have been included in the model


3) go to LGCM_contrasts.m
