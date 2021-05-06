function[R_money] = R_amounts(n_R_levels)
%[R_money] = R_amounts(n_R_levels)
% R_amounts will create a structure with the planned amount for each
% reward level
%
% INPUTS
% n_R_levels: number of reward levels
%
% OUTPUTS
% R_money: structure with 1 subfield for each reward level

switch n_R_levels
    case 3
        R_money.R_1 = 0.2;
        disp(['Reward level 1 = ',num2str(R_money.R_1),' chf']);
        R_money.R_2 = 0.5;
        disp(['Reward level 2 = ',num2str(R_money.R_2),' chf']);
        R_money.R_3 = 1.0;
        disp(['Reward level 3 = ',num2str(R_money.R_3),' chf']);
    case 4
        R_money.R_1 = 0.2;
        disp(['Reward level 1 = ',num2str(R_money.R_1),' chf']);
        R_money.R_2 = 0.5;
        disp(['Reward level 2 = ',num2str(R_money.R_2),' chf']);
        R_money.R_3 = 1.0;
        disp(['Reward level 3 = ',num2str(R_money.R_3),' chf']);
        R_money.R_4 = 2.0;
        disp(['Reward level 4 = ',num2str(R_money.R_4),' chf']);
    otherwise
        error(['Please prepare Reward level - Money mapping for ',...
            num2str(n_R_levels),' reward levels.']);
end

% also define the amount of money for the loss trials
R_money.trialFail = 5;
disp(['Loss amount = ',num2str(R_money.trialFail),' chf for failure trials']);

end % function