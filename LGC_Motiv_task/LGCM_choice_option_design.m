function[choice_opt] = LGCM_choice_option_design(n_R_levels, n_E_levels, punishment_yn, n_trials)
%[choice_opt] = LGCM_choice_option_design(n_R_levels, n_E_levels, punishment_yn, n_trials)
% LGCM_choice_option_design will create the design matrix of the trials to
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
% OUTPUTS
% choice_opt: structure with the details of the choice options

%% initialize variables of interest
[choice_opt.LR_side,...
    choice_opt.R_levels,...
    choice_opt.E_levels] = deal(NaN(2, n_trials));

if strcmp(punishment_yn,'yes')
    choice_opt.R_or_P = cell(1, n_trials); % variable to know which trials are rewarding and which are punishments
end

%% punishment included or not?
switch punishment_yn
    case 'yes'
        n_trials_per_condition = n_trials/2;
        % rewards
        choice_design_RE = LGCM_choice_design_RE_levels(n_R_levels, n_E_levels, n_trials_per_condition);
        % punishment
        choice_design_PE_tmp = LGCM_choice_design_RE_levels(n_R_levels, n_E_levels, n_trials_per_condition);
        % for punishment, reward left and reward right should be inverted
        % otherwise the options will be strange (more high punishments with
        % high levels of effort)
        choice_design_PE.effort = choice_design_PE_tmp.effort;
        choice_design_PE.allOptions = [choice_design_PE_tmp.allOptions(:,2),...
            choice_design_PE_tmp.allOptions(:,1),...
            choice_design_PE_tmp.allOptions(:,3),...
            choice_design_PE_tmp.allOptions(:,4)]; % reverse Rleft and Rright
        choice_design_PE.punishment.left = choice_design_PE_tmp.reward.right;
        choice_design_PE.punishment.right = choice_design_PE_tmp.reward.left;
        % pool all together
        choice_design.allOptions = [choice_design_RE.allOptions, choice_design_PE.allOptions];
        choice_design.reward.left   = [choice_design_RE.reward.left, choice_design_PE.punishment.right];
        choice_design.reward.right  = [choice_design_RE.reward.right, choice_design_PE.punishment.left];
        choice_design.effort.left   = [choice_design_RE.effort.left, choice_design_PE.effort.left];
        choice_design.effort.right  = [choice_design_RE.effort.right, choice_design_PE.effort.right];
        
        % variable to keep track of reward vs punishment trials
        R_or_P_tmp = [repmat({'R'},1,n_trials_per_condition), repmat({'P'},1,n_trials_per_condition)]; % (0 for punishment and 1 for reward)
        
    case 'no'
        choice_design = LGCM_choice_design_RE_levels(n_R_levels, n_E_levels, n_trials);
end

%% randomize order of the trials
rdm_order = randperm(n_trials);
choice_opt.R.left = choice_design.reward.left(rdm_order);
choice_opt.R.right = choice_design.reward.right(rdm_order);
choice_opt.E.left = choice_design.effort.left(rdm_order);
choice_opt.E.right = choice_design.effort.right(rdm_order);

% if punishments included, should also apply randomization to
% reward/punishment case
if strcmp(punishment_yn,'yes')
    choice_opt.R_or_P = R_or_P_tmp(rdm_order);
end

end % function