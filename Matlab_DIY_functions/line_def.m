function [ line_hdl,  ref_line_vals] = line_def( x_lims, y_lims,...
    col_line, line_size)
% [ line_hdl,  ref_line_vals] = line_def( x_lims, y_lims,...
%     col_line, line_size)
% Even if x_vals and y_vals are different, will draw a straight line for
% the perfect correlation y = x
%
% INPUTS
% x_lims, y_lims: values for the x and y axis where you want to draw the
% line (if left empty, the script will use the extremities of the graph
%
% col_line: color identification, if left empty, will be black ('k') by
% default
%
% line_size: line size, if left empty, will be set to 3 by default
%
% OUTPUTS
% line_hdl: line handle
%
% ref_line_vals: value of the x and y axis extremities for the line
%
% See also refline & line

%% line parameters

% by default, if no x_lims and y_lims entered, take the whole graph values
if ~exist('x_lims','var') || isempty(x_lims)
    x_lims = xlim();
end
if ~exist('y_lims','var') || isempty(y_lims)
    y_lims = ylim();
end

% black line by default
if ~exist('col_line','var') || isempty(col_line)
    col_line = 'k';
end

% size of the line by default
if ~exist('line_size','var') || isempty(line_size)
    line_size = 3;
end

%% extract limits of the graph
ref_line_vals = [min(x_lims,y_lims),...
    max(x_lims, y_lims)];

%% draw graph
line_hdl = line(ref_line_vals, ref_line_vals,...
    'Color',col_line,'LineWidth',line_size);

end % function