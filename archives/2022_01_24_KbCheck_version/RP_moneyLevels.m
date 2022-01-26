function [choiceOptions] = RP_moneyLevels(choiceOptions, R_money)
%[choiceOptions] = RP_moneyLevels(choiceOptions, R_money)
% RP_moneyLevels will add a field in choiceOptions structure with the
% corresponding monetary levels for the left/right option
%
% INPUTS
% choiceOptions: structure with reward/punishment and effort levels +
% reward or punishment trials
%
% R_money: structure with information about reward/punishment level and
% monetary equivalent
%
% OUTPUTS
% choiceOptions: structure updated with one additional field for the actual
% monetary amounts

%% identify reward/punishment trials
RP_trials = choiceOptions.R_or_P;
n_trials = length(RP_trials);

%% add monetary amount
[moneyLeft,...
    moneyRight] = deal( NaN(1,n_trials) );
for iTrial = 1:n_trials
    incentiveLeft_tmp = choiceOptions.R.left(iTrial);
    incentiveRight_tmp = choiceOptions.R.right(iTrial);
    switch RP_trials{iTrial}
        case 'P'
            moneyLeft(iTrial)   = R_money.(['P_',num2str(incentiveLeft_tmp)]);
            moneyRight(iTrial)  = R_money.(['P_',num2str(incentiveRight_tmp)]);
        case 'R'
            moneyLeft(iTrial)   = R_money.(['R_',num2str(incentiveLeft_tmp)]);
            moneyRight(iTrial)  = R_money.(['R_',num2str(incentiveRight_tmp)]);
    end
end % trial loop

%% add the field
choiceOptions.monetary_amount.left     = moneyLeft;
choiceOptions.monetary_amount.right    = moneyRight;

end % function