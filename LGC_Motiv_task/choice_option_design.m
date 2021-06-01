function[choice_opt] = choice_option_design(n_R_levels, n_E_levels, punishment_yn, n_trials, R_money)
%[choice_opt] = choice_option_design(n_R_levels, n_E_levels, punishment_yn, n_trials)
% choice_option_design will create the design matrix of the trials to
% be presented at each trial:
% condition: punishment or reward?
% reward and effort level associated to each option
% side: which option will appear on the left and which option will appear
% on the right side of the screen
%
% INPUTS
% n_R_levels: number of reward levels
%
% n_E_levels: number of effort levels
%
% punishment_yn
% 'yes': include punishment trials
% 'no': no punishment
%
% n_trials: total number of trials for the current session
%
% R_money: structure with monetary amount corresponding to each reward and
% punishment level
%
% OUTPUTS
% choice_opt: structure with the details of the choice options

%% punishment included or not?
switch punishment_yn
    case 'yes'
        n_trials_per_condition = n_trials/2;
        % rewards
        choice_design_RE = choice_design_RE_levels(n_R_levels, n_E_levels, n_trials_per_condition, 'R');
        % punishment: be careful matrix is not defined the same way (high
        % punishment left is like low reward left)
        choice_design_PE = choice_design_RE_levels(n_R_levels, n_E_levels, n_trials_per_condition, 'P');
        % pool all together
        choice_design.allOptions = [choice_design_RE.allOptions, choice_design_PE.allOptions];
        choice_design.reward.left   = [choice_design_RE.reward.left, choice_design_PE.punishment.left];
        choice_design.reward.right  = [choice_design_RE.reward.right, choice_design_PE.punishment.right];
        choice_design.effort.left   = [choice_design_RE.effort.left, choice_design_PE.effort.left];
        choice_design.effort.right  = [choice_design_RE.effort.right, choice_design_PE.effort.right];
        
        % variable to keep track of reward vs punishment trials
        R_or_P_tmp = [repmat({'R'},1,n_trials_per_condition), repmat({'P'},1,n_trials_per_condition)]; % (0 for punishment and 1 for reward)
        
    case 'no'
        choice_design = choice_design_RE_levels(n_R_levels, n_E_levels, n_trials, 'R');
end

%% randomize order of the trials
rdm_order = randperm(n_trials);
choice_opt.R.left = choice_design.reward.left(rdm_order);
choice_opt.R.right = choice_design.reward.right(rdm_order);
choice_opt.E.left = choice_design.effort.left(rdm_order);
choice_opt.E.right = choice_design.effort.right(rdm_order);

% if punishments included, should also apply randomization to
% reward/punishment case
switch punishment_yn
    case 'yes'
        choice_opt.R_or_P = R_or_P_tmp(rdm_order);
    case 'no' % all trials are rewarding
        choice_opt.R_or_P = repmat({'R'},1,n_trials);
end

%% transform reward (and punishment) levels into monetary amounts
[choice_design.reward_amount.left, choice_design.reward_amount.right] = deal(NaN(1,n_trials));
for iTrial = 1:n_trials
    switch choice_opt.R_or_P{iTrial}
        case 'R'
            choice_design.monetary_amount.left(iTrial) = R_money.(['R_',num2str(choice_opt.R.left(iTrial))]);
            choice_design.monetary_amount.right(iTrial) = R_money.(['R_',num2str(choice_opt.R.right(iTrial))]);
        case 'P'
            choice_design.monetary_amount.left(iTrial) = R_money.(['P_',num2str(choice_opt.R.left(iTrial))]);
            choice_design.monetary_amount.right(iTrial) = R_money.(['P_',num2str(choice_opt.R.right(iTrial))]);
    end
end

end % function