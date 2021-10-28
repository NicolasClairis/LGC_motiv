function[R_money] = R_amounts_IP(n_R_levels, punishment_yn, IPdata, effort_type)
%[R_money] = R_amounts_IP(n_R_levels, punishment_yn, IPdata, effort_type)
% R_amounts will create a structure with the planned amount for each
% reward level
%
% INPUTS
% n_R_levels: number of reward levels
%
% punishment_yn: 'yes'/'no': does the script include punishments as well?
%
% IPdata: structure containing:
% - the indifference point = amount of reward for which the medium effort
% level is equivalent to the default low option for rewards
% - information about amount of reward 
%
% effort_type: 'mental' or 'physical' effort task
%
% OUTPUTS
% R_money: structure with 1 subfield for each reward level

%% load baseline added by default to all monetary amounts.
baselineR = IPdata.baselineR;
baselineP = IPdata.baselineP;

%% load delta between default option and medium effort level for which the
% two options are perceived as equivalent (ie indifference point IP)
delta_IP = IPdata.([effort_type,'DeltaIP']);

%% rewards
% extract value for default option
R_money.R_0 = baselineR;
% extract value for indifference point (corresponding to middle reward
% level)
IP_R = baselineR + delta_IP;
if delta_IP < 0.01
    error('indifference point is too low');
end
switch n_R_levels
    case 3
        R_money.R_1 = IP_R - delta_IP/2;
        R_money.R_2 = IP_R + delta_IP/2;
    case 4
        R_money.R_1 = IP_R - delta_IP/2;
        R_money.R_2 = IP_R;
        R_money.R_3 = IP_R + delta_IP/2;
    case 5
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
    disp(['Reward level ',num2str(iR),' = ',num2str(R_money.(['R_',num2str(iR-1)])),' chf']);
end

%% punishments
% extract value for default option
R_money.P_0 = baselineP;
% extract value for indifference point (corresponding to middle punishment
% level)
IP_P = baselineP - delta_IP;
if strcmp(punishment_yn,'yes')
    switch n_R_levels
        case 3
            R_money.P_1 = IP_P - delta_IP/2;
            R_money.P_2 = IP_P + delta_IP/2;
        case 4
            R_money.P_1 = IP_P - delta_IP/2;
            R_money.P_2 = IP_P;
            R_money.P_3 = IP_P + delta_IP/2;
        case 5
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
        disp(['Punishment level ',num2str(iP),' = ',num2str(R_money.(['P_',num2str(iP-1)])),' chf']);
    end
end

end % function