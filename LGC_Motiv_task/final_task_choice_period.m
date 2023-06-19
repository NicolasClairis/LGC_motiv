function[choice_trial, onsetDispChoiceOptions, onsetChoice, stoptask] = final_task_choice_period(scr, stim,...
    R_left, R_right, E_left, E_right, R_or_P,...
    timeParameter, key, confidenceDisp)
% [choice_trial, onsetDispChoiceOptions, onsetChoice, stoptask] = final_task_choice_period(scr, stim,...
%     R_left, R_right, E_left, E_right, R_or_P,...
%     timeParameter, key, confidenceDisp)
% final_task_choice_period will display the choice options and then wait for the 
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
% timeParameter: structure with timing information
%   .timeLimit: (false) no time limit; (true) time limit for the choice
%   period
%   .t_choice: maximal time to wait for choice
%
% key: code for left/right keys
%
% confidenceDisp: true/false to know whether the mapping of the confidence
% should be displayed or not
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

% % if only 2 button for answer or no input, no confidence mapping display
% if (key.n_buttonsChoice == 2) || ~exist('confidenceDisp','var') || (isempty(confidenceDisp))
%     confidenceDisp = false;
% end
    
%% ask question on top
DrawFormattedText(window, stim.choice.choiceQuestion.text, stim.choice.choiceQuestion.x, stim.choice.choiceQuestion.y, stim.choice.choiceQuestion.colour);
DrawFormattedText(window, stim.choice.choiceOR.text, stim.choice.choiceOR.x, stim.choice.choiceOR.y, stim.choice.choiceOR.colour);

%% display each difficulty level
leftStartAngle = stim.difficulty.startAngle.(['level_',num2str(E_left)]);
rightStartAngle = stim.difficulty.startAngle.(['level_',num2str(E_right)]);
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
switch R_or_P
    case 'R'
        DrawFormattedText(window,stim.choice.win.text,...
            stim.winRewardText.top_left(1),...
            stim.winRewardText.top_left(2),...
            white);
        DrawFormattedText(window,stim.choice.win.text,...
            stim.winRewardText.top_right(1),...
            stim.winRewardText.top_right(2),...
            white);
    case 'P'
        DrawFormattedText(window,stim.choice.lose.text,...
            stim.loseRewardText.top_left(1),...
            stim.loseRewardText.top_left(2),...
            white);
        DrawFormattedText(window,stim.choice.lose.text,...
            stim.loseRewardText.top_right(1),...
            stim.loseRewardText.top_right(2),...
            white);
end
% display monetary amount associated to each option
drawRewardAmount(scr, stim, R_left, R_or_P, 'top_left_start');
drawRewardAmount(scr, stim, R_right, R_or_P, 'top_right_start');

% display corresponding effort text
DrawFormattedText(window,stim.choice.for.text,...
    stim.effort_introText.bottom_left(1),...
    stim.effort_introText.bottom_left(2),...
    white);
DrawFormattedText(window,stim.choice.for.text,...
    stim.effort_introText.bottom_right(1),...
    stim.effort_introText.bottom_right(2),...
    white);

%% add mapping for confidence
if confidenceDisp == true
    % left sure
    DrawFormattedText(window,stim.leftSure.text,...
        stim.leftSure.x, stim.leftSure.y,...
        stim.leftSure.colour);
    % left unsure
    DrawFormattedText(window,stim.leftUnsure.text,...
        stim.leftUnsure.x, stim.leftUnsure.y,...
        stim.leftUnsure.colour);
    % right unsure
    DrawFormattedText(window,stim.rightUnsure.text,...
        stim.rightUnsure.x, stim.rightUnsure.y,...
        stim.rightUnsure.colour);
    % right sure
    DrawFormattedText(window,stim.rightSure.text,...
        stim.rightSure.x, stim.rightSure.y,...
        stim.rightSure.colour);
end

%% display everything on the screen and record the timing
[~,onsetDispChoiceOptions] = Screen('Flip',window);

%% wait for choice to be made or time limit to be reached
choicePeriodOver = 0;
while choicePeriodOver == 0
    choice_trial = 0; % by default no choice is being made
    %% check time if a time limit is set for the choice
    if timeParameter.timeLimit == true
        timeNow = GetSecs;
        if timeNow > (onsetDispChoiceOptions + timeParameter.t_choice)
            % finish the trial
            choicePeriodOver = 1;
            onsetChoice = NaN;
        end
    end
    
    %% check key press
    [keyisdown, ~, ~, lastPress, ~] = KbQueueCheck;
    
    %% some key was pressed
    if keyisdown == 1
        if key.n_buttonsChoice == 2
            %% left option chosen
            if lastPress(key.left) > onsetDispChoiceOptions &&...
                    lastPress(key.right) < onsetDispChoiceOptions
                % record time of chosen option
                timedown = lastPress(key.left);
                % record side of chosen option
                choice_trial = -1;
                choicePeriodOver = 1;
                %% right option chosen
            elseif lastPress(key.left) < onsetDispChoiceOptions &&...
                    lastPress(key.right) > onsetDispChoiceOptions
                % record time of chosen option
                timedown = lastPress(key.right);
                % record side of chosen option
                choice_trial = 1;
                choicePeriodOver = 1;
                %% stop the task
            elseif lastPress(key.escape) > 0
                choicePeriodOver = 1;
                stoptask = 1;
            end
            
        elseif key.n_buttonsChoice == 4
            %% LEFT SURE option chosen
            if lastPress(key.leftSure) > onsetDispChoiceOptions &&...
                    lastPress(key.leftUnsure) < onsetDispChoiceOptions &&...
                    lastPress(key.rightUnsure) < onsetDispChoiceOptions &&...
                    lastPress(key.rightSure) < onsetDispChoiceOptions
                % record time of chosen option
                timedown = lastPress(key.leftSure);
                % record side of chosen option
                choice_trial = -2;
                choicePeriodOver = 1;
                %% LEFT UNSURE option chosen
            elseif lastPress(key.leftSure) < onsetDispChoiceOptions &&...
                    lastPress(key.leftUnsure) > onsetDispChoiceOptions &&...
                    lastPress(key.rightUnsure) < onsetDispChoiceOptions &&...
                    lastPress(key.rightSure) < onsetDispChoiceOptions
                % record time of chosen option
                timedown = lastPress(key.leftUnsure);
                % record side of chosen option
                choice_trial = -1;
                choicePeriodOver = 1;
                %% RIGHT UNSURE option chosen
            elseif lastPress(key.leftSure) < onsetDispChoiceOptions &&...
                    lastPress(key.leftUnsure) < onsetDispChoiceOptions &&...
                    lastPress(key.rightUnsure) > onsetDispChoiceOptions &&...
                    lastPress(key.rightSure) < onsetDispChoiceOptions
                % record time of chosen option
                timedown = lastPress(key.rightUnsure);
                % record side of chosen option
                choice_trial = 1;
                choicePeriodOver = 1;
                %% RIGHT SURE option chosen
            elseif lastPress(key.leftSure) < onsetDispChoiceOptions &&...
                    lastPress(key.leftUnsure) < onsetDispChoiceOptions &&...
                    lastPress(key.rightUnsure) < onsetDispChoiceOptions &&...
                    lastPress(key.rightSure) > onsetDispChoiceOptions
                % record time of chosen option
                timedown = lastPress(key.rightSure);
                % record side of chosen option
                choice_trial = 2;
                choicePeriodOver = 1;
                %% stop the task
            elseif lastPress(key.escape) > 0
                choicePeriodOver = 1;
                stoptask = 1;
            end
        end
    end % some key was pressed
    
    %% get time when choice was made
    if choice_trial ~= 0
        onsetChoice = timedown;
    end
end % choice period

end % function