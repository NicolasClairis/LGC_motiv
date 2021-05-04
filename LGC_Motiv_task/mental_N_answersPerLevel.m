function[n_to_reach] = LGCM_mental_N_answersPerLevel(n_E_levels)
% [n_to_reach] = LGCM_mental_N_answersPerLevel(n_E_levels)
% LGCM_mental_N_answersPerLevel will determine the number of correct
% answers to provide for each difficulty level depending on the total
% number of difficulty levels that you want to implement
% 
% INPUTS
% n_E_levels: number of effort levels to use in the task
%
% OUTPUTS
% n_to_reach: structure with the corresponding number of correct answers to
% provide for each difficulty level
%
% See also LGCM_choice_task_main.m

if n_E_levels < 3
    error('not enough difficulty levels!');
else
    switch n_E_levels
        case 3
            n_to_reach.E_level_1 = 2;
            n_to_reach.E_level_2 = 4;
            n_to_reach.E_level_3 = 6;
        case 4
            n_to_reach.E_level_1 = 2;
            n_to_reach.E_level_2 = 4;
            n_to_reach.E_level_3 = 6;
            n_to_reach.E_level_4 = 8;
        case 5
            n_to_reach.E_level_1 = 2;
            n_to_reach.E_level_2 = 4;
            n_to_reach.E_level_3 = 6;
            n_to_reach.E_level_4 = 8;
            n_to_reach.E_level_5 = 10;
        case 6
            n_to_reach.E_level_1 = 2;
            n_to_reach.E_level_2 = 4;
            n_to_reach.E_level_3 = 6;
            n_to_reach.E_level_4 = 8;
            n_to_reach.E_level_5 = 10;
            n_to_reach.E_level_6 = 12;
        otherwise
            error('not ready yet with so many levels of effort');
    end % difficulty levels
    
end % filter when not enough levels of difficulty

end % function