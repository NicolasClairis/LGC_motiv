function [ hdl ] = fig()
%[ hdl ] = fig()
% fig creates a figure occupying all the screen.
%
% OUTPUT
% hdl: figure handle
%
% See also figure
%
% Written by Nicolas Clairis - august 2019 (in Matlab 2017a)

%% change matlab font from ugly grey to white
set(0,'defaultfigurecolor',[1 1 1]);
%% create figure
hdl = figure;
%% force "hold on" to add multiple plots eventually
hold on;
%% maximize window size
matlabVersion = version('-release'); % extract matlab version
matlabYearVersion = str2double(matlabVersion(1:4)); % extract the year
if matlabYearVersion > 2018 % later versions of Matlab can work with this
    hdl.WindowState = 'maximized'; % maximize window size
else % alternatively, use this code to maximize the window size
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
end

end % function