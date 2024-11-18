function[txt1_hdl, txt2_hdl, txtSize] = place_r_and_pval_2(r_corr1, pval1, col1, r_corr2, pval2, col2)
%[txt1_hdl, txt2_hdl, txtSize] = place_r_and_pval_2(r_corr1, pval1, col1, r_corr2, pval2, col2)
% place_r_and_pval_2 serves to display correlations info on graph. Similar
% to place_r_and_pval.m but will show 2 correlation coeficients instead of
% 1.
%
% INPUTS
% r_corr1: correlation coefficient for correlation 1
%
% pval1: p.value for correlation 1
%
% col1: color for correlation 1
%
% r_corr2: correlation coefficient for correlation 2
%
% pval2: p.value for correlation 2
%
% col2: color for correlation 2
%
% OUTPUTS
% txt1_hdl: text handle for first correlation
%
% txt2_hdl: text handle for second correlation
%
% txtSize: text size to use for handle
%
% See also place_r_and_pval and demo_place_r_and_pval_2

%% round the values in case they are too long
% round correlation coefficients (if below 0.001, will be equal to 0.001 or
% 0)
% round first correlation coefficient
if r_corr1 > 0.01
    r_corr1 = round(r_corr1, 2);
else
    r_corr1 = round(r_corr1,3);
end
% round second correlation coefficient
if r_corr2 > 0.01
    r_corr2 = round(r_corr2, 2);
else
    r_corr2 = round(r_corr2,3);
end

% include a number of stars depending on if p.values are significant or not
% for correlation 1
if pval1 > 0.05
    pval1_stars = '';
elseif pval1 > 0.01 && pval1 <= 0.05
    pval1_stars = '*';
elseif pval1 > 0.005 && pval1 <= 0.01
    pval1_stars = '**';
elseif pval1 > 0.001 && pval1 <= 0.005
    pval1_stars = '***';
elseif pval1 <= 0.001
    pval1_stars = '****';
elseif strcmp(pval1,'NaN')
    pval1_stars = 'NaN';
end
r_corr1_str = [num2str(r_corr1),pval1_stars];

% for correlation 2
if pval2 > 0.05
    pval2_stars = '';
elseif pval2 > 0.01 && pval2 <= 0.05
    pval2_stars = '*';
elseif pval2 > 0.005 && pval2 <= 0.01
    pval2_stars = '**';
elseif pval2 > 0.001 && pval2 <= 0.005
    pval2_stars = '***';
elseif pval2 <= 0.001
    pval2_stars = '****';
elseif strcmp(pval2,'NaN')
    pval2_stars = 'NaN';
end
r_corr2_str = [num2str(r_corr2),pval2_stars];

% for p.value need to track how small p.value is if below 0.001
if pval1 >= 0.001
    pval1_str = ['p = ',num2str(round(pval1, 3))];
else
    pval1_str = 'p < 0.001';
end
if pval2 >= 0.001
    pval2_str = ['p = ',num2str(round(pval2, 3))];
else
    pval2_str = 'p < 0.001';
end

%% define x coordinate for the text
xlim_vals = xlim();
xlim_dims = xlim_vals(2) - xlim_vals(1);
if (r_corr1 < 0 && r_corr2 < 0) || (r_corr1 >= 0 && r_corr2 >= 0) % right of the plot
    x_val_txt = xlim_vals(2);
    x_alignment = 'right';
else % center
    x_val_txt = (xlim_dims/2) + xlim_vals(1);
    x_alignment = 'center';
end
%% define y coordinate depending on if correlations are negative or positive (to avoid overlapping relevant values)
ylim_vals = ylim();
if r_corr1 < 0 && r_corr2 < 0 % negative correlations (top of the screen)
    y_val_txt1 = ylim_vals(2);
    y_alignement = 'cap';
elseif r_corr1 >= 0 && r_corr2 >= 0 % positive correlations (bottom of the screen)
    % ylim_dims = ylim_vals(2) - ylim_vals(1);
    y_val_txt1 = ylim_vals(1);
    y_alignement = 'bottom';
else % one is positive and the other is negative: top of the screen
    y_val_txt1 = ylim_vals(2);
    y_alignement = 'cap';
end

%% text size
txtSize = 20;

%% add the text
txt1_hdl = text(x_val_txt,y_val_txt1,...
    {['r = ',r_corr1_str],pval1_str},...
    'FontSize',txtSize,'Color',col1,...
    'HorizontalAlignment',x_alignment,...
    'VerticalAlignment',y_alignement);
% text.Extent = [left bottom width height] rectangle containing text infos
txt1_y_size = txt1_hdl.Extent(4);
if y_val_txt1 == ylim_vals(2) % top
    y_val_txt2 = y_val_txt1 - txt1_y_size;
elseif y_val_txt1 == ylim_vals(1) % bottom
    y_val_txt2 = y_val_txt1 + txt1_y_size;
end
txt2_hdl = text(x_val_txt,y_val_txt2,...
    {['r = ',r_corr2_str],pval2_str},...
    'FontSize',txtSize,'Color',col2,...
    'HorizontalAlignment',x_alignment,...
    'VerticalAlignment',y_alignement);

end % function