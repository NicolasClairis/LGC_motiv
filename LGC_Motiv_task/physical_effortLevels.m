function [Ep_time_levels] = physical_effortLevels(n_E_levels)
%[Ep_time_levels] = physical_effortLevels(n_E_levels)
% physical_effortLevels will determine the duration (in seconds)
% corresponding to each difficulty level according to the total number of
% difficulty conditions (n_E_levels) that has been defined.
%
% INPUTS
% n_E_levels: number of effort levels
%
% OUTPUTS
% Ep_time_levels: structure with duration of effort to perform for each
% effort difficulty
% 


switch n_E_levels
    case 3
        Ep_time_levels.level_1 = 1.5;
        Ep_time_levels.level_2 = 3;
        Ep_time_levels.level_3 = 4.5;
    case 4
        Ep_time_levels.level_1 = 0.5;
        Ep_time_levels.level_2 = 1.5;
        Ep_time_levels.level_3 = 2.5;
        Ep_time_levels.level_4 = 3.5;
end

end % function