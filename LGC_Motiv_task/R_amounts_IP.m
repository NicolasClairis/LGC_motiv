function[R_money] = R_amounts_IP(n_R_levels, punishment_yn, delta_IP)
%[R_money] = R_amounts_IP(n_R_levels, punishment_yn, IP_R)
% R_amounts will create a structure with the planned amount for each
% reward level
%
% INPUTS
% n_R_levels: number of reward levels
%
% punishment_yn: 'yes'/'no': does the script include punishments as well?
%
% IP_R: indifference point = amount of reward for which the medium effort
% level is equivalent to the default low option for rewards
%
% OUTPUTS
% R_money: structure with 1 subfield for each reward level

%% baselineR: baseline added by default to all gains and also starting point
% for losses.
% for gains, means you will win at least this
% for punishments, means you will lose not more than this
baselineR = 0.50;

%% rewards
R_money.R_default = baselineR;
% delta = IP_R - baselineR;
IP_R = baselineR - delta_IP;
if delta_IP < 0.01
    error('indifference point is too low');
end
switch n_R_levels
    case 2
        R_money.R_1 = IP_R - delta_IP/2;
        R_money.R_3 = IP_R + delta_IP/2;
    case 3
        R_money.R_1 = IP_R - delta_IP/2;
        R_money.R_2 = IP_R;
        R_money.R_3 = IP_R + delta_IP/2;
    case 4
        R_money.R_1 = IP_R - delta_IP;
        R_money.R_2 = IP_R - delta_IP/2;
        R_money.R_3 = IP_R + delta_IP/2;
        R_money.R_4 = IP_R + delta_IP;
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
R_money.P_default = baselineR;
IP_P = baselineR - delta_IP;
if strcmp(punishment_yn,'yes')
    switch n_R_levels
        case 2
            R_money.P_1 = IP_P - delta_IP/2;
            R_money.P_3 = IP_P + delta_IP/2;
        case 3
            R_money.P_1 = IP_P - delta_IP/2;
            R_money.P_2 = IP_P;
            R_money.P_3 = IP_P + delta_IP/2;
        case 4
            R_money.P_1 = IP_P - delta_IP/2;
            R_money.P_2 = IP_P - delta_IP/4;
            R_money.P_3 = IP_P + delta_IP/4;
            R_money.P_4 = IP_P + delta_IP/2;
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

end % function