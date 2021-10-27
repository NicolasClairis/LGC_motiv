function[n_to_reach] = mental_N_answersPerLevel(n_E_levels, NmaxPerf)
% [n_to_reach] = mental_N_answersPerLevel(n_E_levels, NmaxPerf)
% mental_N_answersPerLevel will determine the number of correct
% answers to provide for each difficulty level depending on the total
% number of difficulty levels that you want to implement
%
% INPUTS
% n_E_levels: number of effort levels to use in the task
%
% NmaxPerf: number of correct answers reached during calibration
%
% OUTPUTS
% n_to_reach: structure with the corresponding number of correct answers to
% provide for each difficulty level
%
% See also choice_task_main.m

%% define max level as a given percentage of their actual maximum
E_max = floor(NmaxPerf*(80/100)) + 1;

%% define different difficulty levels based on the calibration
switch n_E_levels
    case 3
        n_to_reach.E_level_0 = floor(E_max*(1/9)) + 1;
        n_to_reach.E_level_1 = floor(E_max*(1/2)) + 1;
        n_to_reach.E_level_2 = E_max;
    case 4
        n_to_reach.E_level_0 = floor(E_max*(1/9)) + 1;
        n_to_reach.E_level_1 = floor(E_max*(1/3)) + 1;
        n_to_reach.E_level_2 = floor(E_max*(2/3)) + 1;
        n_to_reach.E_level_3 = E_max;
    case 5
        n_to_reach.E_level_0 = floor(E_max*(1/8)) + 1;
        n_to_reach.E_level_1 = floor(E_max*(1/4)) + 1;
        n_to_reach.E_level_2 = floor(E_max*(1/2)) + 1;
        n_to_reach.E_level_3 = floor(E_max*(3/4)) + 1;
        n_to_reach.E_level_4 = E_max;
    otherwise
        error([num2str(n_E_levels),' effort levels not ready yet.']);
end

%% check that levels of effort are not identical
for iEffort = 0:(n_E_levels - 1)
    if iEffort ~=  n_E_levels - 1
        if n_to_reach.(['E_level_',num2str(iEffort)]) == n_to_reach.(['E_level_',num2str(iEffort+1)])
            error(['Mental effort levels ',num2str(iEffort),' and ',num2str(iEffort + 1),' are identical. You have to find a way to solve this']);
        end
    end
end % effort loop

end % function