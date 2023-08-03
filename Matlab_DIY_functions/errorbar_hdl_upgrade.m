function[] = errorbar_hdl_upgrade(errorbar_hdl, color_to_use)
% [] = errorbar_hdl_upgrade(errorbar_hdl, color_to_use)
% errorbar_hdl_upgrade improves display for errorbar graphs.
%
% INPUTS
% errorbar_hdl: handle of errorbar curve to be improved
%
% color_to_use: color to use for the errorbar. In case it's not defined
% (because there is only one errorbar to show for example), it will be
% displayed in black by default

%% load general parameters
[~, lWidth, col, mSize] = general_fig_prm;

%% improve display
errorbar_hdl.LineStyle = 'none';
if ~exist('color_to_use','var') || isempty(color_to_use)
    errorbar_hdl.Color = col.black;
else
    errorbar_hdl.Color = color_to_use;
end
errorbar_hdl.Marker = 'o';
errorbar_hdl.MarkerSize = mSize;
errorbar_hdl.LineWidth = lWidth;

end % function