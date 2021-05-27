function[time_dispChoice, R_chosen, E_chosen] = choice_task_dispChosen(scr, stim, R_chosen, E_chosen, R_or_P)
% [time_dispChoice, R_chosen, E_chosen] = choice_task_dispChosen(scr, stim, choice_opt, choice,...
%     R_or_P, iTrial)
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with stimuli parameters (reward and effort display
% informations are stored in here)
%
% R_chosen: reward amount of the chosen option
%
% E_chosen: effort level of the chosen option
%
% R_or_P: character indicating the nature of the current trial
% 'R': reward trial
% 'P': punishment trial
%
% OUTPUTS
% time_dispChoice: onset of the display of the chosen option on the screen
%
% R_chosen: numeric variable reflecting the level of the reward associated
% to the chosen option
%
% E_chosen:numeric variable reflecting the level of the effort associated
% to the chosen option
%
% See also choice_task_main.m

%% extract relevant parameters
window = scr.window;
xScreenCenter = scr.xCenter;
yScreenCenter = scr.yCenter;
white = scr.colours.white;
black = scr.colours.black;

% remind the option they chose
DrawFormattedText(window,'Vous avez choisi','center',yScreenCenter/6.5,white);

%% display reward and effort level
switch R_chosen
    case 0 % no option was selected
        DrawFormattedText(window,...
            'Trop lent!',...
            'center', yScreenCenter/2, white);
    otherwise % one option was selected
        
        % if punishment trial, add also indication to know that money is to be lost
        switch R_or_P
            case 'R'
                DrawFormattedText(window,'Gagner',...
                    stim.winRewardText.top_center,...
                    white);
                moneySign = '+';
                moneyColour = stim.reward.text.colour;
            case 'P'
                DrawFormattedText(window,'Perdre',...
                    stim.loseRewardText.top_center,...
                    white);
                moneySign = '-';
                moneyColour = stim.punishment.text.colour;
        end
        
        trialMoneyObtained = sprintf('%0.2f',R_chosen);
        DrawFormattedText(window,[moneySign, trialMoneyObtained,' CHF'],...
            stim.reward.text.top_center_start,...
            moneyColour);
        
        %% display difficulty level
        chosenStartAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen)]);
        maxCircleAngle = stim.difficulty.arcEndAngle;
        Screen('FillArc', window,...
            stim.difficulty.currLevelColor,...
            stim.chosenOption.difficulty,...
            chosenStartAngle,...
            maxCircleAngle - chosenStartAngle);% chosen option difficulty
        
        % display maximal difficulty level for each option (= full circle)
        Screen('FrameOval', window, stim.difficulty.maxColor,...
            stim.chosenOption.difficulty,...
            stim.difficulty.ovalWidth);
        
        % display
        DrawFormattedText(window,'pour',...
            stim.effort_introText.bottom_center,...
            white);
end

%% display a square on top of selected reward and effort
Screen('FrameRect', window,...
    stim.chosenOption.squareColour,...
    stim.chosenOption.squareRect,...
    stim.chosenOption.squareWidth);

%% display on screen and extract timing
[~,time_dispChoice] = Screen('Flip',window);

end % function