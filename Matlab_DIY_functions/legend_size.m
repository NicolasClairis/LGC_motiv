function [  ] = legend_size( pSize )
%legend_size( pSize )
% sets the scale of the xlabel, ylabel and legends at pSize value
%
% INPUTS
% pSize: size of the text you want to use
% Originally written by N.Clairis - 31/01/2020

% change font
set(gca,'FontName','Calibri');
set(gca,'fontsize',pSize,'FontWeight','normal');
set(findall(gcf,'type','text'),'FontSize',pSize,'fontWeight','normal');
% change matlab font from ugly grey to white
set(0,'defaultfigurecolor',[1 1 1]);

end % function