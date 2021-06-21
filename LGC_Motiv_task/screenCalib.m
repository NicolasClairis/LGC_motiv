function[RectCoordinates] = screenCalib()
%[RectCoordinates] = screenCalib()
% script to define the window seen by the subject in the crappy CIBM 7T scanner
%
% OUTPUTS
% RectCoordinates: structure with left/right and lower/upper coordinates
%

%% open a PTB window
IRM = 1;
testing = 1;
[scr, xScreenCenter, yScreenCenter, window] = ScreenConfiguration(IRM, testing);
white = scr.colours.white;

%% define keys to check
KbName('UnifyKeyNames');
increaseHeight  = KbName('8');
reduceHeight     = KbName('2');
increaseWidth   = KbName('6');
reduceWidth      = KbName('4');
moveToRight     = KbName('RightArrow');
moveToLeft      = KbName('LeftArrow');
moveToTop       = KbName('UpArrow');
moveToBottom    = KbName('DownArrow');
escape_press    = KbName('escape');

%% draw a rectangle the size of the screen
xStart = 0;
yStart = 0;
xEnd = xScreenCenter*2;
yEnd = yScreenCenter*2;
drawRectangle(window, white, xStart, yStart, xEnd, yEnd);

%% resize the rectangle until it fits the visible screen
pxlReSize = 1;
endReSize = false;

while endReSize == false
   [keyIsDown, ~, keyCode] = KbCheck;
   
   if keyIsDown == 1
       if keyCode(increaseHeight)
           yStart  = yStart - pxlReSize; % start upper
           yEnd    = yEnd + pxlReSize; % finish lower
       elseif keyCode(reduceHeight)
           yStart  = yStart + pxlReSize; % start lower
           yEnd    = yEnd - pxlReSize; % finish upper
       elseif keyCode(moveToBottom) % move rectangle to bottom
           yStart  = yStart + pxlReSize;
           yEnd    = yEnd + pxlReSize;
       elseif keyCode(moveToTop) % move rectangle up
           yStart  = yStart - pxlReSize;
           yEnd    = yEnd - pxlReSize;
       elseif keyCode(increaseWidth)
           xStart = xStart - pxlReSize; % start more on the left
           xEnd = xEnd + pxlReSize; % finish more on the right
       elseif keyCode(reduceWidth)
           xStart = xStart + pxlReSize; % start more on the right
           xEnd = xEnd - pxlReSize; % finish more on the left
       elseif keyCode(moveToLeft) % move rectangle to the left
           xStart = xStart - pxlReSize;
           xEnd = xEnd - pxlReSize;
       elseif keyCode(moveToRight) % move rectangle to the right
           xStart = xStart + pxlReSize;
           xEnd = xEnd + pxlReSize;
       elseif keyCode(escape_press)
           endReSize = true;
       end
        
        drawRectangle(window, white, xStart, yStart, xEnd, yEnd);
    end
end % loop

%% store coordinates of the rectangle
RectCoordinates.xStart = xStart;
RectCoordinates.yStart = yStart;
RectCoordinates.xEnd = xEnd;
RectCoordinates.yEnd = yEnd;
% rectangle center
RectCoordinates.xCenter = mean([xStart, xEnd]);
RectCoordinates.yCenter = mean([yStart, yEnd]);

%% close PTB screen
sca;

end % function

function[]=drawRectangle(window, colour, xStart, yStart, xEnd, yEnd)
lWidth = 3;
Screen('FrameRect',window, colour, [xStart yStart xEnd yEnd], lWidth);
% Screen('FillOval', window, colour, [xScreenCenter, yScreenCenter]);
Screen('Flip', window);
end