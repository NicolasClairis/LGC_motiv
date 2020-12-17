function[simu_choices, simu_stim] = LGCM_simulate_data(k_sensitivity_prm,...
    stim, nTrials)
%[simu_choices, simu_stim] = LGCM_simulate_data(k_sensitivity_prm,...
%     stim, nTrials)
% LGCM_simulate_data will simulate data for a simplistic model of
% motivation inside the LGCM task.
%
% SV = kR*(1+R_level) - kC*(1+kF*trial number)*E_levels
% SV = -kP*(1+P_level) - kC*(1+kF*trial number)*E_levels
% => choice(A) = sigm(SV(A)-SV(B))
%
% INPUTS
% k_sensitivity_prm:
% - kR, kP, kC, kF: average sensitivity to reward, punishment
% cost and fatigue across participants => add noise to each parameter
% - beta_choice: average choice stochasticity
%
% stim: structure with information about stimuli
% .fixed_or_evolving:
%   'fix': fixed stimuli
%   'evolving': simulate a staircase procedure
% .R_levels/.E_levels .low/high: details of the reward and effort levels
% define as vectors for each of the two possible options
%
% nTrials: number of total trials
%
% OUTPUTS
% simu_choices: vector with 0/1 for choices
% Rewards:
% 0= low effort low reward
% 1 = high effort high reward
% Punishments:
% 0 = low punishment high effort
% 1 = high punishment low effort
%
% simu_stim: structure with the stimulus of each trial in terms of reward
% and effort level

%% motivational parameters
kR = k_sensitivity_prm.kR + rand;
kP = k_sensitivity_prm.kP + rand;
kC = k_sensitivity_prm.kC + rand;
kF = k_sensitivity_prm.kF + rand;

%% choice stochasticity
beta_choice = k_sensitivity_prm.beta_choice + rand;

%% stimuli type = task design
task_trials_type = stim.fixed_or_evolving;
switch task_trials_type
    case 'fix'
        R_levels = stim.R_levels;
        E_levels = stim.E_levels;
end

%% prepare for simu and outputs
simu_choices = NaN(1,nTrials);



% for iTrial = 1:nTrials
%     R_hRhE = ;
%     R_lRlE = ;
%     E_hRhE = ;
%     E_lRlE = ;
%     VA = LGCM_simu_get_SV(kR, kP, kC, kF, R, E, iTrial, model_n);
%     VB = LGCM_simu_get_SV(kR, kP, kC, kF, R, E, iTrial, model_n);
%     simu_choices(iTrial) = sigmo( (VA-VB), 1/beta_choice);
% end % trial loop

end % function