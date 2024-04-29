function[] = scat_hdl_upgrade(scat_hdl, markerEdge_color, markerFill_color)
% [] = scat_hdl_upgrade(scat_hdl, markerEdge_color, markerFill_color)
% scat_hdl_upgrade improves display for scatter plots.
%
% INPUTS
% scat_hdl: handle of dots to be improved
%
% markerEdge_color: color to use for the edges of the marker. In case it's not defined
% (because there is only one errorbar to show for example), it will remain
% empty (no edges).
%
% markerFill_color: color to use to fill inside the marker. In case it's
% not defined, it will by shown in grey by default.

%% load general parameters
[~, lWidth, col, mSize] = general_fig_prm;

%% improve display
% marker edges color
if ~exist('markerEdge_color','var') || isempty(markerEdge_color)
    scat_hdl.MarkerEdgeColor = 'none';
else
    scat_hdl.MarkerEdgeColor = markerEdge_color;
end
% marker filling
if ~exist('markerFill_color','var') || isempty(markerFill_color)
    scat_hdl.MarkerFaceColor = col.grey;
else
    scat_hdl.MarkerFaceColor = markerFill_color;
end
scat_hdl.Marker = 'o';
scat_hdl.SizeData = 70;
scat_hdl.LineWidth = 1.5;

end % function