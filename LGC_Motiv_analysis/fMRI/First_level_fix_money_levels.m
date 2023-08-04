function[money_level_left, money_level_right] = First_level_fix_money_levels(study_nm, sub_nm, task_fullName,...
    money_level_left_v0, money_level_right_v0, RP_var)
% [money_level_left, money_level_right] = First_level_fix_money_levels(study_nm, sub_nm, task_fullName,...
%     money_level_left_v0, money_level_right_v0, RP_var)
% First_level_fix_money_levels re-arranges reward and punishment levels to
% both be within range 1 (least motivating/default option) and 4 (more
% motivating reward or punishment for the varying effort option). The
% script also fixes the data for 2 subjects where there were only 3 levels
% instead of 4.
%
% INPUTS
% study_nm: study name (to know whether data needs fixing or not)
% 
% sub_nm: subject name (to know whether data needs fixing or not)
%
% task_fullName: run task name ('mental'/'physical') (to know whether data needs fixing or not)
%
% money_level_left_v0: level of money for the left option before fixing
%
% money_level_right_v0: level of money for the right option before fixing
%
% RP_var: variable equal to (+1) for rewards and to (-1) for punishments
%
% OUTPUTS
% money_level_left: money level for the left option between 1 (default)
% and 4 (more motivating reward/punishment)
%
% money_level_right: money level for the right option between 1 (default)
% and 4 (more motivating reward/punishment)

%% initialize vars of interest
money_level_left = money_level_left_v0;
money_level_right = money_level_right_v0;

%% increase all reward levels to have a range between 1 and 4 instead of
% between 0 and 3:
% R1 = default low R;
% R4 = highest more rewarding R
money_level_left(money_level_left_v0 == 3 & RP_var == 1) = 4;
money_level_left(money_level_left_v0 == 2 & RP_var == 1) = 3;
money_level_left(money_level_left_v0 == 1 & RP_var == 1) = 2;
money_level_left(money_level_left_v0 == 0 & RP_var == 1) = 1;
money_level_right(money_level_right_v0 == 3 & RP_var == 1) = 4;
money_level_right(money_level_right_v0 == 2 & RP_var == 1) = 3;
money_level_right(money_level_right_v0 == 1 & RP_var == 1) = 2;
money_level_right(money_level_right_v0 == 0 & RP_var == 1) = 1;

%% increase all punishment levels to have a range between 1 and 4 instead 
% of between 0 and 3 + swap original P1 and original P3 so that the range 
% goes from more aversive (default P1) to least aversive (P4).
% P1 = default high P;
% P4 = lowest and least aversive P (and therefore more motivating P)

% least aversive option P1 becomes P4
money_level_left(money_level_left_v0 == 1 & RP_var == -1) = 4;
money_level_right(money_level_right_v0 == 1 & RP_var == -1) = 4;
% second least aversive option P2 becomes P3
money_level_left(money_level_left_v0 == 2 & RP_var == -1) = 3;
money_level_right(money_level_right_v0 == 2 & RP_var == -1) = 3;
% more aversive varying option P3 becomes P2
money_level_left(money_level_left_v0 == 3 & RP_var == -1) = 2;
money_level_right(money_level_right_v0 == 3 & RP_var == -1) = 2;
% default option P0 becomes P1
money_level_left(money_level_left_v0 == 0 & RP_var == -1) = 1;
money_level_right(money_level_right_v0 == 0 & RP_var == -1) = 1;


%% fix data for subjects where IP had a bug and two options were based on
% the same amount
if strcmp(study_nm,'study1')
    if strcmp(sub_nm,'064') && strcmp(task_fullName,'mental')
        
        % R4 is equal to R3 (rewards only go from R1 to R3 for this
        % subject)
        money_level_left(money_level_left == 4 & RP_var == 1) = 3;
        money_level_right(money_level_right == 4 & RP_var == 1) = 3;
        
        % P4 is also equal to P3 (punishments only go from P1 to P3 for
        % this subject)
        money_level_left(money_level_left == 4 & RP_var == -1) = 3;
        money_level_right(money_level_right == 4 & RP_var == -1) = 3;
    elseif strcmp(sub_nm,'090') && strcmp(task_fullName,'physical')
        
        % R4 is equal to R3 (rewards only go from R1 to R3 for this
        % subject)
        money_level_left(money_level_left == 4 & RP_var == 1) = 3;
        money_level_right(money_level_right == 4 & RP_var == 1) = 3;
        
        % P4 is also equal to P3 (punishments only go from P1 to P3 for
        % this subject)
        money_level_left(money_level_left == 4 & RP_var == -1) = 3;
        money_level_right(money_level_right == 4 & RP_var == -1) = 3;
    end
end

end % function