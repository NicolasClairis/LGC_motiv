function[] = scat_hdl_upgrade(scat_hdl, color_to_use)
% [] = scat_hdl_upgrade(scat_hdl, color_to_use)
% scat_hdl_upgrade improves display for scatter plots.
%
% INPUTS
% scat_hdl: handle of dots to be improved
%
% color_to_use: color to use for the errorbar. In case it's not defined
% (because there is only one errorbar to show for example), it will be
% displayed in black by default

%% load general parameters
[~, lWidth, col, mSize] = general_fig_prm;

%% improve display
if ~exist('color_to_use','var') || isempty(color_to_use)
    scat_hdl.MarkerEdgeColor = col.black;
else
    scat_hdl.MarkerEdgeColor = color_to_use;
end
scat_hdl.Marker = 'o';
scat_hdl.SizeData = 70;
scat_hdl.LineWidth = 1.5;

end % function