function[pSize, lWidth, col, mSize] = general_fig_prm
% [pSize, lWidth, col, mSize] = general_fig_prm
% general_fig_prm will set up general figure parameters
%
% OUTPUTS
% pSize: legend size
%
% lWidth: line width for curves
%
% col: structure with several colours
%
% mSize: size for markers in scatter plot

%% legend size
pSize = 30;
%% line width
lWidth = 3;

%% marker size
mSize = 20;

%% initialize colours
col.white = [1 1 1];
col.black = [0 0 0];
col.grey = [174 174 174]./255;
% 4 colours compatible with daltonian perception
% based on https://colorbrewer2.org/#type=diverging&scheme=Spectral&n=4
col.red = [215 25 28]./255;
col.orange = [253 174 97]./255;
col.green = [171 221 164]./255;
col.blue = [43 131 186]./255;

end % function