function [ handle_signif ] = place_signif_curve( fig_hdl, x_signif, col, ypos_perc, signif_size  )
%[ handle ] = place_signif_curve( fig_hdl, x_signif, col, ypos_perc, signif_size  )
% place_signif_curve places stars on the top part of the graph based on the
% size of the graph.
%
% INPUTS
%
% fig_hdl: figure handle
%
% x_signif: vector with values of the significant clusters in x.axis
%
% col: colour for significant clusters representation
%
% ypos_perc: percentage of the y.axis on which you want to represent the
% significant curve
%
% signif_size: size of the dots you want to represent
%
% OUTPUTS
%
% handle_signif: handle of the curve representing the significant data
%
% Written by N.Clairis 11/01/2019

%% extract size of y axis
ylim_range = ylim;
ylim_size = abs( ylim_range(2) - ylim_range(1) );

%% define position on the y. axis
if ~exist('ypos_perc','var') || isempty(ypos_perc) || isnan(ypos_perc)
    ypos = 1/100*ylim_size;
else
    ypos = ypos_perc*ylim_size;
end
y_index_signifClust = ylim_range(2) - ypos;
y_axis_data = ones(length(x_signif),1).*y_index_signifClust;

%% colour
if ~exist('col','var') || isempty(col) || isnan(col)
    col = 'k';
end

%% size of dots
if ~exist('signif_size','var') || isempty(signif_size) || isnan(signif_size)
    signif_size = 5;
end

%% plot the data
figure(fig_hdl);
handle_signif = plot(x_signif, y_axis_data,...
    '*','Color',col,'LineWidth',signif_size);

end

