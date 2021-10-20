function [timeDispWait] = waitSpace(langage, window, yScreenCenter, scr, keys)
% [timeDispWait] = waitSpace(langage, window, yScreenCenter, scr, keys)
% function to wait for space press before moving on
%
% INPUTS
% langage: 'fr' french/'engl' english
%
% window: PTB window index
%
% yScreenCenter: coordinates of y to put text
%
% scr: screen colours, etc
%
% keys: code for key space stored in this structure
%
% OUTPUTS
% timeDispWait: time when message to wait for space press is displayed on
% screen

%% display wait for space on screen
switch langage
    case 'fr'
        DrawFormattedText(window,...
            'Appuyez sur espace quand vous etes pret(e) a demarrer.',...
            'center', yScreenCenter*(5/3), scr.colours.white, scr.wrapat);
    case 'engl'
        DrawFormattedText(window,...
            'Press space key when you are ready to start.',...
            'center', yScreenCenter*(5/3), scr.colours.white, scr.wrapat);
end
[~, timeDispWait] = Screen(window,'Flip');
disp('Please press space');

%% wait for space press
[~, ~, keyCode] = KbCheck();
while(keyCode(keys.space) ~= 1)
    % wait until the key has been pressed
    [~, ~, keyCode] = KbCheck();
end

end % function