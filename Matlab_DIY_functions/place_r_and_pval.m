function[txt_hdl, txtSize] = place_r_and_pval(r_corr, pval)
% place_r_and_pval serves to display correlation info on graph for the
% mediation script (but can also be used elsewhere)
%
% INPUTS
% r_corr: correlation coefficient
%
% pval: p.value
%
% See also mediation
%
% OUTPUTS
% txt_hdl: text handle
%
% txtSize: text size to use for handle

%% round the values in case they are too long
% round correlation coefficient (if below 0.001, will be equal to 0.001 or
% 0)
if r_corr > 0.01
    r_corr = round(r_corr, 2);
else
    r_corr = round(r_corr,3);
end
% include a number of stars depending on if pval significant or not
if pval > 0.05
    pval_stars = '';
elseif pval > 0.01 && pval <= 0.05
    pval_stars = '*';
elseif pval > 0.005 && pval <= 0.01
    pval_stars = '**';
elseif pval > 0.001 && pval <= 0.005
    pval_stars = '***';
elseif pval <= 0.001
    pval_stars = '****';
end
r_corr_str = [num2str(r_corr),pval_stars];

% for p.value need to track how small p.value is if below 0.001
if pval > 0.001
    pval_str = num2str(round(pval, 3));
else
    pval_str = num2str(round(pval,1,'significant'));
end

%% define x coordinate for the text
xlim_vals = xlim();
xlim_dims = xlim_vals(2) - xlim_vals(1);
x_val_txt = xlim_vals(2) - xlim_dims/7;

%% define y coordinate depending on if correlation is negative or positive (to avoid overlapping relevant values)
if r_corr < 0 % negative correlation (top right of the screen)
    ylim_vals = ylim();
    y_val_txt = ylim_vals(2);
elseif r_corr >= 0 % positive correlation (bottom right of the screen)
    ylim_vals = ylim();
    ylim_dims = ylim_vals(2) - ylim_vals(1);
    y_val_txt = ylim_vals(1) + ylim_dims/7;
end

%% text size
txtSize = 17;

%% add the text
txt_hdl = text(x_val_txt,y_val_txt,...
    {['r = ',r_corr_str],...
    ['p = ',pval_str]},...
    'FontSize',txtSize);

end % function