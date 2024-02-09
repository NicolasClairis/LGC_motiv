function [ fig_hdl ] = fig()
%[ fig_hdl ] = fig()
% fig creates a figure occupying all the screen.
%
% OUTPUT
% fig_hdl: figure handle
%
% See also figure
%
% Written by Nicolas Clairis - august 2019 (in Matlab 2017a)

%% set default settings for figure layout
[pSize, lWidth, col, mSize, fontName] = general_fig_prm; % load general figure parameters
set(0,'defaultfigurecolor',col.white); % change matlab font from ugly grey to white
set(0,'defaultLineLineWidth',lWidth); % set the default line width to lWidth
set(0,'defaultLineMarkerSize',mSize); % set the default line marker size to mSize
% set text features
set(0,'defaultTextFontName',fontName); % set the default police font name
set(0,'defaultTextFontSize',pSize); % set the default police font size
% set axes features
set(0,'defaultAxesFontName',fontName); % set the default police font name
set(0,'defaultAxesFontSize',pSize); % set the default police font size
% set legend features
set(0,'defaultLegendFontName',fontName); % set the default police font size for legends
set(0,'defaultLegendFontSize',pSize); % set the default police font size for legends
set(0,'defaultLegendBox','off'); % remove legend box by default
% check default_matlab_figure_fields.m if you want to check more fields

%% create figure + extract handle
fig_hdl = figure;

%% force "hold on" to add multiple plots eventually
hold on;

%% set up color palette
f = getframe(fig_hdl);
colormap(f.colormap);

%% maximize window size
matlabVersion = version('-release'); % extract matlab version
matlabYearVersion = str2double(matlabVersion(1:4)); % extract the year
if matlabYearVersion > 2018 % later versions of Matlab can work with this
    fig_hdl.WindowState = 'maximized'; % maximize window size
else % alternatively, use this code to maximize the window size
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
end

end % function