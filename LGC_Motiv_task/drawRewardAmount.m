function [] = drawRewardAmount(scr, stim, R_amount, R_or_P, xyCoord)
% [] = drawRewardAmount(scr, stim, R_amount, R_or_P, xyCoord)
% drawRewardAmount draws the money on screen depending on the input
% parameters
%
% INPUTS
% scr: structure with screen parameters (and baseline text display size)
%
% stim: structure with info about money text display size
% 
% R_or_P: reward or punishment case? Will need to adapt the sign and the
% colour of the stimulus accordingly
%
% xyCoord: vector with x coordinate (1) and y coordinate (2) where to start
% displaying the monetary amount
%
%

%% main parameters
window = scr.window;
moneyTextSize = stim.reward.textSizeForPTB;
baselineTextSize = scr.textSize.baseline;
% calibrate reward number to be in the "+X.XX CHF" format
switch R_or_P
    case 'R'
        moneySign = '+';
        moneyColour = stim.reward.text.colour;
    case 'P'
        moneySign = '-';
        moneyColour = stim.punishment.text.colour;
end
trialMoneyObtained = sprintf('%0.2f',R_amount );

%% adapt text size for rewards to appear bigger
Screen('TextSize', window, moneyTextSize);

%% display money won/lost
DrawFormattedText(window, [moneySign, trialMoneyObtained,' CHF'],...
    xyCoord(1),...
    xyCoord(2),...
    moneyColour);

%% reset baseline text size
Screen('TextSize', window, baselineTextSize);

end % function