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

hdl = figure; % create figure
if ismember(version('-release'),{'2018a','2018b','2019a','2019b'})
    hdl.WindowState = 'maximized'; % maximize window size
else
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
end

end % function