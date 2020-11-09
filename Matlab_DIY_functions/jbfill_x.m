function[fillhandle,mean_hdl,msg]=jbfill_x(xpoints,upper,lower,ypoints,color,edge,add,transparency)
%USAGE: [fillhandle,mean_hdl,msg]=jbfill_x(xpoints,upper,lower,ypoints,color,edge,add,transparency)
%This function will fill a region with a color between the two vectors provided
%using the Matlab fill command.
%
%fillhandle is the returned handle to the filled region in the plot.
%mean_hdl is the handle for the mean curve in the middle of the plot.
%xpoints= The horizontal data points (ie frequencies). Note length(Upper)
%         must equal Length(lower)and must equal length(xpoints)!
%upper = the upper curve values (data can be less than lower)
%lower = the lower curve values (data can be more than upper)
%ypoints = the vertical y points
%color = the color of the filled area and of the middle curve
%edge  = the color around the edge of the filled area
%add   = a flag to add to the current plot or make a new one.
%transparency is a value ranging from 1 for opaque to 0 for invisible for
%the filled color only.
%
%John A. Bockstege November 2006;
%Example:
%     a=rand(1,20);%Vector of random data
%     b=a+2*rand(1,20);%2nd vector of data points;
%     x=1:20;%horizontal vector
%     [ph,msg]=jbfill(x,a,b,rand(1,3),rand(1,3),0,rand(1,1))
%     grid on
%     legend('Datr')
%
% NICOLAS CLAIRIS and JULES BROCHARD EDIT (8/11/16): all the NaN values are removed from the
% "upper" points, the "lower" points and the corresponding "xpoint" index
%
%
% N.Clairis edit (29/05/20): added 'Color' in
% plot(xpoints,middle,'Color',color) (line 66) to be able to define color
% values
%
% alternative to jbfill to draw a curve with shaded area in the x axis
% while jbfill.m is made to do this only on the y-axis
%
% See also jbfill

if nargin<8;transparency=.5;end %default is to have a transparency of .5
if nargin<7;add=1;end     %default is to add to current plot
if nargin<6;edge='k';end  %dfault edge color is black
if nargin<5;color='b';end %default color is blue

upper_valid = ~isnan(upper);
lower_valid = ~isnan(lower);

if any(~upper_valid) || any(~lower_valid)
    warning('Function jbfill.m: NaN values were detected in the variable ''upper'' or ''lower'', they have been removed for display.')
end

if length(upper) == length(lower) && length(lower) == length(xpoints) %&& any(upper_valid ~= lower_valid)
    msg = '';
    filled = [lower(lower_valid), fliplr(upper(upper_valid))];
    ypointsfill = [ypoints(lower_valid),fliplr(ypoints(upper_valid))];
    if add
        hold on
    end
    fillhandle = fill(filled,ypointsfill,color);%plot the data
    set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
    if add
        hold off
    end
else
    msg='Error: Must use the same number of points in each vector';
end

if nargin >= 4
    hold on;
    mean_hdl = plot(xpoints,ypoints,'Color',color);
end
