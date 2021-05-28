function[choice_trial, onsetDispChoiceOptions, onsetChoice, stoptask] = choice_period(scr, stim,...
    R_left, R_right, E_left, E_right, R_or_P,...
    t_choice, key)
% [choice_trial, onsetDispChoiceOptions, onsetChoice] = choice_period(scr, stim, choice_opt,...
%     R_left, R_right, E_left, E_right, R_or_P, R_amounts,...
%     t_choice, key)
% choice_period will display the choice options and then wait for the 
% choice to be made (or the time limit to be reached. Provides timings and 
% choice made in output.
%
% INPUTS
% scr: structure with screen informations
%
% stim: structure with informations about the stimuli to display
%
% R_left, R_right: reward monetary amount for left and right option
%
% E_left, E_right: level of effort for left and right option
%
% R_or_P: character indicating the nature of the current trial
% 'R': reward trial
% 'P': punishment trial
%
% t_choice: maximal time to wait for choice
%
% key: code for left/right keys
%
% OUTPUTS
% choice_trial:
% (-1): left option chosen
% (0): no choice made
% (+1): right option chosen
%
% onsetDispChoiceOptions: onset of when the options appear on screen
%
% onsetChoice: if a choice was made, displays the timing
%
% stoptask:
% (0) keep on
% (1) in this case: signal for main function to stop the task

%% initialize variables of interest
window = scr.window;
stoptask = 0;
white = scr.colours.white;
% black = scr.colours.black;
yScreenCenter = scr.yCenter;
% xScreenCenter = scr.xCenter;

%% ask question on top
DrawFormattedText(window,'Que préférez-vous?','center',yScreenCenter/4,white);
DrawFormattedText(window,'OU','center','center',white);

%% display each difficulty level
leftStartAngle = stim.difficulty.startAngle.(['level_',num2str(E_left)]);
rightStartAngle = stim.difficulty.startAngle.(['level_',num2str(E_right)]);
maxCircleAngle = stim.difficulty.arcEndAngle;
Screen('FillArc', window,...
    stim.difficulty. currLevelColor,...
    stim.difficulty.below_left,...
    leftStartAngle,...
    maxCircleAngle - leftStartAngle); % left option difficulty
 Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
    stim.difficulty.below_right,...
    rightStartAngle,...
    maxCircleAngle - rightStartAngle);% right option difficulty

% display maximal difficulty level for each option (= full circle)
Screen('FrameOval', window, stim.difficulty.maxColor,...
    stim.difficulty.below_left,...
    stim.difficulty.ovalWidth);
Screen('FrameOval', window, stim.difficulty.maxColor,...
    stim.difficulty.below_right,...
    stim.difficulty.ovalWidth);

%% display each monetary incentive level
switch R_or_P
    case 'R'
        DrawFormattedText(window,'Gagner',...
            stim.winRewardText.top_left,...
            white);
        DrawFormattedText(window,'Gagner',...
            stim.winRewardText.top_right,...
            white);
        moneySign = '+';
        moneyColour = stim.reward.text.colour;
    case 'P'
        DrawFormattedText(window,'Gagner',...
            stim.loseRewardText.top_left,...
            white);
        DrawFormattedText(window,'Gagner',...
            stim.loseRewardText.top_right,...
            white);
        moneySign = '-';
        moneyColour = stim.punishment.text.colour;
end
% extract money corresponding to left and right option in the X.XX format
moneyLeft = sprintf('%0.2f',R_left);
moneyRight = sprintf('%0.2f',R_right);
DrawFormattedText(window,...
    [moneySign, moneyLeft,' CHF'],...
    stim.reward.text.top_left_start, moneyColour);
DrawFormattedText(window,...
    [moneySign, moneyRight,' CHF'],...
    stim.reward.text.top_right_start, moneyColour);

% display corresponding effort text
DrawFormattedText(window,'pour',...
    stim.effort_introText.bottom_left,...
    white);
DrawFormattedText(window,'pour',...
    stim.effort_introText.bottom_right,...
    white);

[~,onsetDispChoiceOptions] = Screen('Flip',window);

%% wait for choice to be made or time limit to be reached
choicePeriodOver = 0;
while choicePeriodOver == 0
    %% check time
    timeNow = GetSecs;
    choice_trial = 0; % by default no choice is being made
    if timeNow > (onsetDispChoiceOptions + t_choice)
        % finish the trial
        choicePeriodOver = 1;
        onsetChoice = NaN;
    end
    
    %% check key press
    [keyisdown, secs, keycode] = KbCheck;
    
    %% some key was pressed
    if keyisdown == 1
        %% left option chosen
        if keycode(key.left) == 1 &&...
                keycode(key.right) == 0
            % record time of chosen option
            timedown = secs;
            % record side of chosen option
            choice_trial = -1;
            choicePeriodOver = 1;
            %% right option chosen
        elseif keycode(key.left) == 0 &&...
                keycode(key.right) == 1
            % record time of chosen option
            timedown = secs;
            % record side of chosen option
            choice_trial = 1;
            choicePeriodOver = 1;
        %% stop the task
        elseif keycode(key.escape) == 1
            choicePeriodOver = 1;
            stoptask = 1;
        end
    end % some key was pressed
    
    %% get time when choice was made
    if choice_trial ~= 0
        onsetChoice = timedown;
    end
end % choice period

%     warning('you can either leave it that way (no visual feedback at the moment they chose an option',...
%         ' or show some visual feedback of the selected option already',...
%         ' or require the subject to keep pressing the button until the end of the trial.');

end % function