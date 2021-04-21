function [choice_R, choice_E] = LGCM_choice_design_RE_levels(n_R_levels, n_E_levels, n_trials)
%[choice_R, choice_E] = LGCM_choice_design_RE_levels(n_R_levels, n_E_levels, n_trials)
%
% INPUTS
% n_R_levels: number of reward levels
%
% n_E_levels: number of effort levels
%
% n_trials: total number of trial of the current session
%
% OUTPUTS
% choice_R: reward levels for all options
%
% choice_E: effort levels for all options
%
% See also LGCM_choice_option_design.m

%%
if n_R_levels == n_E_levels
   if mod(n_R_levels, 2) == 0 % pair number of conditions (easier to split)
       
       %%
       error('Check Jules method to optimize the space sampled');
       
       %% vary the distance but keep the mean constant
       % 1) same reward and effort level option 1 starts at lowRlowE and
       % increases while the other starts at highRhighE and decreases
       vec1_R_a = [1:n_R_levels; n_R_levels:(-1):1];
       vec1_E_a = [1:n_E_levels; n_E_levels:(-1):1];
       
       % 2) Reward and effort level have reverse patterns
       % option 1 starts at lowR-highE and goes toward highR-lowE while
       % option2 does the opposite
       vec1_R_b = [1:n_R_levels;      n_R_levels:(-1):1];
       veca_E_b = [n_E_levels:(-1):1; 1:n_E_levels];
       
       %% 3) vary the mean, but keep distance fixed
       % same reward and effort level for each option
       % distance of 1 between R and E of the 2 options
       vec2_R_a = [1:(n_R_levels - 1); 2:n_R_levels];
       vec2_E_a = [1:(n_E_levels - 1); 2:n_E_levels];
       % distance = median to lowest option (constant)
       vec2_R_b = [1:(n_R_levels/2); ((n_R_levels/2) + 1):n_R_levels];
       vec2_E_b = [1:(n_E_levels/2); ((n_E_levels/2) + 1):n_E_levels];
       
       % reward and effort are not "equal" anymore
       % distance of 1 (but reverse pattern between R and E)
       vec2_R_c = [1:(n_R_levels - 1);    2:n_R_levels];
       vec2_E_c = [2:n_E_levels;          1:(n_E_levels - 1)];
       % distance of 2 (but reverse pattern between R and E)
       vec2_R_d = [1:(n_R_levels - 2);    3:n_R_levels];
       vec2_E_d = [3:n_E_levels;          1:(n_E_levels - 2)];
       
       % pool together
       vec2_R = [vec2_R_a, vec2_R_b, vec2_R_c, vec2_R_d];
       vec2_E = [vec2_E_a, vec2_E_b, vec2_E_c, vec2_E_d];
       
       %% mean and distance constant but reverse R and E values
       % distance of 1 between options
       % 2 levels difference between R and E
       vec3_R_a = [ [1, n_R_levels];        [2, (n_R_levels - 1)]];
       vec3_E_a = [ [(n_E_levels - 1), 2];  [n_E_levels, 1]];
       % distance of 1 between options
       % 1 level difference between R and E
       % low mean
       
       
       % distance of 2 between options
       vec3_R_c = [ [1, n_R_levels];        [3, (n_R_levels - 2)]];
       vec3_E_c = [ [(n_E_levels - 2), 3];  [n_E_levels, 1]];
       
       
       %% pool altogether
       choice_R_tmp = [vec1_R, vec2_R];
       choice_E_tmp = [vec1_E, vec2_E];
       
       %% repeat to match total number of trials
       n_repeats = n_trials/size(choice_R_tmp, 2);
       if floor(n_repeats) ~= n_repeats
           error('Total number of trials does not match ');
       end
       choice_R = repmat(choice_R_tmp, 1, n_repeats);
       choice_E = repmat(choice_E_tmp, 1, n_repeats);
   else
       error('please adapt design for conditions where not pair number of conditions');
   end
    
else
    error('not ready yet for different number of reward and effort levels');
end

end % function