function[col_m1, col_f1, col_m2, col_f2] = col_per_sex()
% [col_m1, col_f1] = col_per_sex()
% col_per_sex attributes a RGB color code for each sex.
%
% OUTPUTS
% col_m1: [r g b] vector for males (light color for scatter)
%
% col_f1: [r g b] vector for females (light color for scatter)
%
% col_m2: [r g b] vector for males (darker color for fit curve)
%
% col_f2: [r g b] vector for females (darker color for fit curve)

% scatter color (lighter)
col_m1 = [166,189,219]./255;
col_f1 = [252,146,114]./255;

% fit color (more pronounced)
col_m2 = [43,140,190]./255;
col_f2 = [222, 45, 38]./255;

end % function