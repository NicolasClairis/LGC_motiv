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

%% define different difficulty levels based on the calibration
switch n_E_levels
    case 3
        n_to_reach.E_level_0 = floor(NmaxPerf*(1/9)) + 1;
        n_to_reach.E_level_1 = floor(NmaxPerf*(1/2)) + 1;
        n_to_reach.E_level_2 = floor(NmaxPerf) + 1;
    case 4
        n_to_reach.E_level_0 = floor(NmaxPerf*(1/9)) + 1;
        n_to_reach.E_level_1 = floor(NmaxPerf*(1/3)) + 1;
        n_to_reach.E_level_2 = floor(NmaxPerf*(2/3)) + 1;
        n_to_reach.E_level_3 = floor(NmaxPerf) + 1;
    case 5
        n_to_reach.E_level_0 = floor(NmaxPerf*(1/8)) + 1;
        n_to_reach.E_level_1 = floor(NmaxPerf*(1/4)) + 1;
        n_to_reach.E_level_2 = floor(NmaxPerf*(1/2)) + 1;
        n_to_reach.E_level_3 = floor(NmaxPerf*(3/4)) + 1;
        n_to_reach.E_level_4 = floor(NmaxPerf) + 1;
    otherwise
        error([num2str(n_E_levels),' effort levels not ready yet.']);
end
    

end % function