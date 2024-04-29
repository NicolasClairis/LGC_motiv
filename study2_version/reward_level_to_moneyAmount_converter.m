function [monetary_amounts] = reward_level_to_moneyAmount_converter(R_levels, Rlevel_to_money, R_or_P)
%[monetary_amounts] = reward_level_to_moneyAmount_converter(R_levels, Rlevel_to_money)
% reward_level_to_moneyAmount_converter will convert R_levels vector of
% reward level values into monetary amounts for each trial based on
% Rlevel_to_money structure
%
% INPUTS
% R_levels: vector of reward (or punishment) level to use on each trial)
%
% Rlevel_to_money: structure with correspondence to know how to convert
% reward levels into monetary amounts
%
% R_or_P: vector indicating if the trial is a reward trial ('R') or a
% punishment trial ('P')
%
% OUTPUTS
% monetary_amounts: vector with monetary amount to use on each trial

n_trials = length(R_levels);
monetary_amounts = NaN(1,n_trials);
for iTrial = 1:n_trials
    switch R_or_P{iTrial}
        case 'R'
            monetary_amounts(iTrial) = Rlevel_to_money.(['R_',num2str(R_levels(iTrial))]);
        case 'P'
            monetary_amounts(iTrial) = Rlevel_to_money.(['P_',num2str(R_levels(iTrial))]);
    end
end

end % function