function[R_money] = R_amounts_IP(n_R_levels, punishment_yn, lowR, delta)
%[R_money] = R_amounts_IP(n_R_levels, punishment_yn, lowR, delta)
% R_amounts will create a structure with the planned amount for each
% reward level
%
% INPUTS
% n_R_levels: number of reward levels
%
% punishment_yn: 'yes'/'no': does the script include punishments as well?
%
% lowR: baseline monetary value which has been used to find the indifference
% point with the staircase procedure
%
% delta: difference between monetary value of both options after reaching
% the indifference point with the staircase procedure
%
% OUTPUTS
% R_money: structure with 1 subfield for each reward level

%% rewards
R_money.R_default = 0.05;
switch n_R_levels
    case 3
        R_money.R_1 = lowR;
        R_money.R_2 = lowR + delta/2;
        R_money.R_3 = lowR + delta;
    case 4
        R_money.R_1 = lowR;
        R_money.R_2 = lowR + delta/2;
        R_money.R_3 = lowR + delta;
        R_money.R_4 = lowR + delta*(3/2);
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
R_money.P_default = 0.05;
if strcmp(punishment_yn,'yes')
    switch n_R_levels
        case 3
            R_money.P_1 = lowR;
            R_money.P_2 = lowR + delta/2;
            R_money.P_3 = lowR + delta;
        case 4
            R_money.P_1 = lowR;
            R_money.P_2 = lowR + delta/2;
            R_money.P_3 = lowR + delta;
            R_money.P_4 = lowR + delta*(3/2);
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
R_money.trialFail = 0.5;
disp(['Loss amount = ',num2str(R_money.trialFail),' chf for failure trials']);

end % function