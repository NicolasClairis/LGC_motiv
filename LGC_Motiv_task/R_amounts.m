function[R_money] = R_amounts(n_R_levels, punishment_yn)
%[R_money] = R_amounts(n_R_levels, punishment_yn)
% R_amounts will create a structure with the planned amount for each
% reward level
%
% INPUTS
% n_R_levels: number of reward levels (including default option)
%
% punishment_yn: 'yes'/'no': does the script include punishments as well?
%
% OUTPUTS
% R_money: structure with 1 subfield for each reward level

%% rewards
switch n_R_levels
    case 3
        R_money.R_0 = 0.50;
        R_money.R_1 = 0.65;
        R_money.R_2 = 0.95;
    case 4
        R_money.R_0 = 0.50;
        R_money.R_1 = 0.65;
        R_money.R_2 = 0.80;
        R_money.R_3 = 0.95;
    otherwise
        error(['Please prepare Reward level - Money mapping for ',...
            num2str(n_R_levels),' reward levels.']);
end

% display level of reward assigned to each amount for tracking for the
% experimenter in case of modification
for iR = 1:n_R_levels
    disp(['Reward level ',num2str(iR-1),' = ',num2str(R_money.(['R_',num2str(iR-1)])),' chf']);
end

%% punishments
if strcmp(punishment_yn,'yes')
    switch n_R_levels
        case 3
            R_money.P_0 = 0.50;
            R_money.P_1 = 0.65;
            R_money.P_2 = 0.95;
        case 4
            R_money.P_0 = 0.50;
            R_money.P_1 = 0.35;
            R_money.P_2 = 0.20;
            R_money.P_3 = 0.05;
        otherwise
            error(['Please prepare Reward level - Money mapping for ',...
                num2str(n_R_levels),' reward levels.']);
    end
    
    % display level of punishment assigned to each amount for tracking for the
    % experimenter in case of modification
    for iP = 1:n_R_levels
        disp(['Punishment level ',num2str(iP-1),' = ',num2str(R_money.(['P_',num2str(iP-1)])),' chf']);
    end
end

end % function