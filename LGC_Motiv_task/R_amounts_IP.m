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
% round the delta IP
delta_IP = round(delta_IP,2);

% half delta IP
half_delta_IP = round(delta_IP/2, 2);

%% rewards
% extract value for default option
R_money.R_0 = baselineR;
% extract value for indifference point (corresponding to middle reward
% level)
IP_R = round(baselineR + delta_IP,2);
if delta_IP < 0.01
    error('indifference point is too low');
end
switch n_R_levels
    case 3
        R_money.R_1 = round(IP_R - half_delta_IP, 2);
        R_money.R_2 = round(IP_R + half_delta_IP, 2);
    case 4
        R_money.R_1 = round(IP_R - half_delta_IP, 2);
        R_money.R_2 = round(IP_R, 2);
        R_money.R_3 = round(IP_R + half_delta_IP, 2);
    case 5
        R_money.R_1 = round(IP_R - delta_IP, 2);
        R_money.R_2 = round(IP_R - half_delta_IP, 2);
        R_money.R_3 = round(IP_R + half_delta_IP, 2);
        R_money.R_4 = round(IP_R + delta_IP, 2);
    otherwise
        error(['Please prepare Reward level - Money mapping for ',...
            num2str(n_R_levels),' reward levels.']);
end

%% fix weird situations

% any level is equal to the next one: increase the distance between reward
% levels
for iR = 1:(n_R_levels - 1)
    if R_money.(['R_',num2str(iR-1)]) == R_money.(['R_',num2str(iR)])
        R_money.(['R_',num2str(iR)]) = R_money.(['R_',num2str(iR)]) + 0.01;
        % increase reward values by the same amount
        if iR < n_R_levels - 1
            for iR_bis = (iR + 1):(n_R_levels - 1)
                R_money.(['R_',num2str(iR_bis)]) = R_money.(['R_',num2str(iR_bis)]) + 0.01;
            end
        end
    end
end

%% display level of reward assigned to each amount for tracking for the
% experimenter in case of modification
for iR = 1:n_R_levels
    disp(['Reward level ',num2str(iR),' = ',num2str(R_money.(['R_',num2str(iR-1)])),' chf']);
end

%% punishments
% extract value for default option
R_money.P_0 = baselineP;
% extract value for indifference point (corresponding to middle punishment
% level)
IP_P = round(baselineP - delta_IP,2);
if strcmp(punishment_yn,'yes')
    switch n_R_levels
        case 3
            R_money.P_1 = round(IP_P - half_delta_IP, 2);
            R_money.P_2 = round(IP_P + half_delta_IP, 2);
        case 4
            R_money.P_1 = round(IP_P - half_delta_IP, 2);
            R_money.P_2 = round(IP_P, 2);
            R_money.P_3 = round(IP_P + half_delta_IP, 2);
        case 5
            R_money.P_1 = round(IP_P - half_delta_IP, 2);
            R_money.P_2 = round(IP_P - delta_IP/4, 2);
            R_money.P_3 = round(IP_P + delta_IP/4, 2);
            R_money.P_4 = round(IP_P + half_delta_IP, 2);
        otherwise
            error(['Please prepare Punishment level - Money mapping for ',...
                num2str(n_R_levels),' punishment levels.']);
    end
    
    %% fix weird situations
    % in case default option is a smaller or equal punishment to the
    % biggest punishment, reduce all by 0.01
    if R_money.(['P_',num2str(n_R_levels - 1)]) >= R_money.P_0
        while R_money.(['P_',num2str(n_R_levels - 1)]) >= R_money.P_0
            for iP = (n_R_levels - 1):(-1):1
                R_money.(['P_',num2str(iP)]) = R_money.(['P_',num2str(iP)]) - 0.01;
            end
        end
    end
    
    
    % any level is equal to the next one: increase the distance between
    % punishment levels
    for iP = (n_R_levels - 1):(-1):2
        if R_money.(['P_',num2str(iP-1)]) == R_money.(['P_',num2str(iP)])
            R_money.(['P_',num2str(iP-1)]) = R_money.(['P_',num2str(iP-1)]) - 0.01;
            % decrease punishment values by the same amount
            if iP < n_R_levels - 1
                for iP_bis = (iP-1):(-1):2
                    R_money.(['P_',num2str(iP_bis)]) = R_money.(['P_',num2str(iP_bis)]) - 0.01;
                end
            end
        end
    end
    
    %% display level of punishment assigned to each amount for tracking for the
    % experimenter in case of modification
    for iP = 1:n_R_levels
        disp(['Punishment level ',num2str(iP),' = ',num2str(R_money.(['P_',num2str(iP-1)])),' chf']);
    end
end

end % function