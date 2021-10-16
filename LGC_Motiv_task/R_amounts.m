function[R_money] = R_amounts(n_R_levels, punishment_yn)
%[R_money] = R_amounts(n_R_levels, punishment_yn)
% R_amounts will create a structure with the planned amount for each
% reward level
%
% INPUTS
% n_R_levels: number of reward levels
%
% punishment_yn: 'yes'/'no': does the script include punishments as well?
%
% OUTPUTS
% R_money: structure with 1 subfield for each reward level

%% rewards
switch n_R_levels
    case 3
        R_money.R_1 = 1.60;
        R_money.R_2 = 1.75;
        R_money.R_3 = 1.90;
    case 4
        R_money.R_1 = 0.20;
        R_money.R_2 = 0.50;
        R_money.R_3 = 1.00;
        R_money.R_4 = 2.00;
    otherwise
        error(['Please prepare Reward level - Money mapping for ',...
            num2str(n_R_levels),' reward levels.']);
end

% display level of reward assigned to each amount for tracking for the
% experimenter in case of modification
for iR = 1:n_R_levels
    disp(['Reward level ',num2str(iR),' = ',num2str(R_money.(['R_',num2str(iR)])),' chf']);
end

%% punishments
if strcmp(punishment_yn,'yes')
    switch n_R_levels
        case 3
            R_money.P_1 = 1.10;
            R_money.P_2 = 1.25;
            R_money.P_3 = 1.40;
        case 4
            R_money.P_1 = 0.2;
            R_money.P_2 = 0.5;
            R_money.P_3 = 1.0;
            R_money.P_4 = 2.0;
        otherwise
            error(['Please prepare Reward level - Money mapping for ',...
                num2str(n_R_levels),' reward levels.']);
    end
    
    % display level of punishment assigned to each amount for tracking for the
    % experimenter in case of modification
    for iP = 1:n_R_levels
        disp(['Punishment level ',num2str(iP),' = ',num2str(R_money.(['P_',num2str(iP)])),' chf']);
    end
end

%% also define the amount of money for the failed trials
R_money.trialFail = 2;
disp(['Loss amount = ',num2str(R_money.trialFail),' chf for failure trials']);

end % function