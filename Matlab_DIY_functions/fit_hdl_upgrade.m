function[] = fit_hdl_upgrade(fit_hdl, color_to_use)
% [] = fit_hdl_upgrade(fit_hdl, color_to_use)
% fit_hdl_upgrade improves display for fit on graphs.
%
% INPUTS
% fit_hdl: handle of fit curve to be improved
%
% color_to_use: color to use for the fit. In case it's not defined
% (because there is only one curve to show for example), it will be
% displayed in grey by default

%% load general parameters
[~, lWidth, col] = general_fig_prm;

%% improve display
fit_hdl.LineStyle = '-';
if ~exist('color_to_use','var') || isempty(color_to_use)
    fit_hdl.Color = col.black;
else
    fit_hdl.Color = color_to_use;
end
fit_hdl.LineWidth = lWidth;

end % function