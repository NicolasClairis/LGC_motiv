function[] = drawMoneyProportional(scr, stim, R_chosen, R_or_P)
%[] = drawMoneyProportional(scr, R_chosen, R_or_P)
% drawMoneyProportional will draw the max and min amount of money to obtain
% (or to lose for punishments) in the current trial.
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with stimuli parameters
%
% R_chosen: reward chosen
%
% R_or_P: 'R' reward or 'P' punishment
%

%% extract screen parameters
window = scr.window;
white = scr.colours.white;

%% display the money
switch R_or_P
    case 'R'
        DrawFormattedText(window, ['+',sprintf('%0.2f',R_chosen)], stim.leftMoneyWinEperf.x, stim.leftMoneyWinEperf.y, white);
        DrawFormattedText(window, '+0.00', stim.rightMoneyWinEperf.x, stim.rightMoneyWinEperf.y, white);
    case 'P'
        DrawFormattedText(window, ['-',sprintf('%0.2f',R_chosen)],     stim_initialize.leftMoneyLoseEperf.x, stim.leftMoneyLoseEperf.y, white);
        DrawFormattedText(window, ['-',sprintf('%0.2f',R_chosen*2)],   stim_initialize.rightMoneyLoseEperf.x, stim.rightMoneyLoseEperf.y, white);
end



end % function