function drawRewardAmount(scr, stim, R_amount, R_or_P, xyCoordField)
% [] = drawRewardAmount(scr, stim, R_amount, R_or_P, xyCoordField)
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
% xyCoordField: field name within stim.reward.amount to use for the x and y
% coordinates where text should be displayed
%
%

%% main parameters
window = scr.window;
moneyTextSize = stim.reward.textSizeForPTB;
baselineTextSize = scr.textSize.baseline;
% calibrate reward number to be in the "+X.XX CHF" format
switch R_or_P
    case 'R'
        RP_type = 'reward';
        moneySign = '+';
    case 'P'
        RP_type = 'punishment';
        moneySign = '-';
end
moneyColour = stim.(RP_type).text.colour;
trialMoneyObtained = sprintf('%0.2f',R_amount );

%% adapt text size for rewards to appear bigger
Screen('TextSize', window, moneyTextSize);

%% display money won/lost
DrawFormattedText(window, [moneySign, trialMoneyObtained,' CHF'],...
    stim.(RP_type).text.(xyCoordField)(1),...
    stim.(RP_type).text.(xyCoordField)(2),...
    moneyColour);

%% reset baseline text size
Screen('TextSize', window, baselineTextSize);

end % function