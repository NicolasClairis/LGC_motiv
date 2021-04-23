function[choice_trial, onsetDispChoiceOptions, onsetChoice, stoptask] = LGCM_choice_period(scr, stim, choice_opt,...
    R_or_P, iTrial, t_choice, key)
% [choice_trial, onsetDispChoiceOptions, onsetChoice] = LGCM_choice_period(scr, stim, choice_opt,...
%     R_or_P, iTrial, t_choice, key)
% LGCM_choice_period will display the choice options and then wait for the 
% choice to be made (or the time limit to be reached. Provides timings and 
% choice made in output.
%
% INPUTS
% scr: structure with screen informations
%
% stim: structure with informations about the stimuli to display
%
% choice_opt: structure with info about choice options (side of each
% option, reward level, effort level, etC.)
%
% R_or_P: character indicating the nature of the current trial
% 'R': reward trial
% 'P': punishment trial
%
% iTrial: trial number
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
black = scr.colours.black;
yScreenCenter = scr.yCenter;
xScreenCenter = scr.xCenter;

%% ask question on top
DrawFormattedText(window,'Que préférez-vous?','center',yScreenCenter/4,white);
DrawFormattedText(window,'OU','center','center',white);

%% extract difficulty & reward level for each side of the screen
% effort level
E_left_tmp = choice_opt.E.left(iTrial);
E_right_tmp = choice_opt.E.right(iTrial);
% extract reward level for each side of the screen
R_left_tmp  = choice_opt.R.left(iTrial);
R_right_tmp = choice_opt.R.right(iTrial);

% E_left_tmp  = 1;
% E_right_tmp = 2;
% % extract reward level for each side of the screen
% R_left_tmp  = 1;
% R_right_tmp = 3;

%% display each difficulty level
leftStartAngle = stim.difficulty.startAngle.(['level_',num2str(E_left_tmp)]);
rightStartAngle = stim.difficulty.startAngle.(['level_',num2str(E_right_tmp)]);
maxCircleAngle = stim.difficulty.arcEndAngle;
Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
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
% extract reward levels
R_left_nm = ['reward_',num2str(R_left_tmp)];
R_right_nm = ['reward_',num2str(R_right_tmp)];

% display reward images
Screen('DrawTexture', window,...
    stim.reward.texture.(R_left_nm),...
    [],...
    stim.reward.top_left.(R_left_nm));
Screen('DrawTexture', window,...
    stim.reward.texture.(R_right_nm),...
    [],...
    stim.reward.top_right.(R_right_nm));

% if punishment trial, add also indication to know that money is to be lost
if strcmp(R_or_P,'P')
   %% 
%    lineWidth = 10;
%    % cross on left option
%    Screen('DrawLine', window, black,...
%        stim.reward.top_left.(R_left_nm)(1), stim.reward.top_left.(R_left_nm)(2),...
%        stim.reward.top_left.(R_left_nm)(3), stim.reward.top_left.(R_left_nm)(4),...
%        lineWidth);
%    Screen('DrawLine', window, black,...
%        stim.reward.top_left.(R_left_nm)(3), stim.reward.top_left.(R_left_nm)(2),...
%        stim.reward.top_left.(R_left_nm)(1), stim.reward.top_left.(R_left_nm)(4),...
%        lineWidth);
%    
%    % cross on right option
%    Screen('DrawLine', window, black,...
%        stim.reward.top_right.(R_right_nm)(1), stim.reward.top_right.(R_right_nm)(2),...
%        stim.reward.top_right.(R_right_nm)(3), stim.reward.top_right.(R_right_nm)(4),...
%        lineWidth);
%    Screen('DrawLine', window, black,...
%        stim.reward.top_right.(R_right_nm)(3), stim.reward.top_right.(R_right_nm)(2),...
%        stim.reward.top_right.(R_right_nm)(1), stim.reward.top_right.(R_right_nm)(4),...
%        lineWidth);

    %% add coloured circle on top of monetary incentives
    % left option
    Screen('FillOval', window, stim.punishment.colourOverlay,...
        stim.punishment.circleOverlay.top_left.(R_left_nm));
    % right option
    Screen('FillOval', window, stim.punishment.colourOverlay,...
        stim.punishment.circleOverlay.top_right.(R_right_nm));
end

switch R_or_P
    case 'R'
        % define coordinates
        xStart_R_left_txt = xScreenCenter/2 - stim.textRectSize.xSizeWin/2;
        xStart_R_right_txt = xScreenCenter*(3/2) - stim.textRectSize.xSizeWin/2;
        yStart_R_txt = stim.reward.top_left.(R_left_nm)(2)-stim.textRectSize.ySizeWin;
        % display
        DrawFormattedText(window,'Gagner',...
            xStart_R_left_txt,...
            yStart_R_txt,...
            white);
        DrawFormattedText(window,'Gagner',...
            xStart_R_right_txt,...
            yStart_R_txt,...
            white);
    case 'P'
        % define coordinates
        xStart_P_left_txt = xScreenCenter/2 - stim.textRectSize.xSizeLose/2;
        xStart_P_right_txt = xScreenCenter*(3/2) - stim.textRectSize.xSizeLose/2;
        yStart_P_txt = stim.reward.top_left.(R_left_nm)(2)-stim.textRectSize.ySizeLose;
        % display
        DrawFormattedText(window,'Perdre',...
            xStart_P_left_txt,...
            yStart_P_txt,...
            white);
        DrawFormattedText(window,'Perdre',...
            xStart_P_right_txt,...
            yStart_P_txt,...
            white);
end
% define coordinates
xStart_forE_left_txt = xScreenCenter/2 - stim.textRectSize.xSizeForEffort/2;
xStart_forE_right_txt = xScreenCenter*(3/2) - stim.textRectSize.xSizeForEffort/2;
yStart_forEffort_txt = stim.difficulty.below_left(2)-stim.textRectSize.ySizeForEffort;
% display
DrawFormattedText(window,'pour',...
    xStart_forE_left_txt,...
    yStart_forEffort_txt,...
    white);
DrawFormattedText(window,'pour',...
    xStart_forE_right_txt,...
    yStart_forEffort_txt,...
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