function[choice_trial, onsetDispChoiceOptions, onsetChoice, stoptask] = LGCM_choice_period(scr, stim,...
    E_list, R_list, iTrial, n_E_levels, t_choice, key)
% [choice_trial, onsetDispChoiceOptions, onsetChoice] = LGCM_choice_period(scr, stim,...
%     E_list, R_list, iTrial, n_E_levels, t_choice, key)
% will display the choice options and then wait for the choice to be made
% (or the time limit to be reached. Provides timings and choice made in
% output.
%
% INPUTS
% scr: structure with screen informations
%
% stim: structure with informations about the stimuli to display
%
% E_list: structure with list of effort levels for each side (left/right)
% for each trial
%
% R_list: same as a E_list but for rewards
%
% iTrial: trial number
%
% n_E_levels: maximum effort levels
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

%% extract difficulty level for each side of the screen
E_left_tmp  = 1; % E_list.left(iTrial)
E_right_tmp = 2; % E_list.right(iTrial)
% extract reward level for each side of the screen
R_left_tmp  = 1; % R_list.left(iTrial)
R_right_tmp = 3; % R_list.right(iTrial)

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

%% display each reward level
Screen('DrawTexture', window,...
    stim.reward.texture.(['reward_',num2str(R_left_tmp)]),...
    [],...
    stim.reward.top_left.(['reward_',num2str(R_left_tmp)]));
Screen('DrawTexture', window,...
    stim.reward.texture.(['reward_',num2str(R_right_tmp)]),...
    [],...
    stim.reward.top_right.(['reward_',num2str(R_right_tmp)]));

[~,onsetDispChoiceOptions] = Screen('Flip',window);

%% wait for choice to be made or time limit to be reached
choicePeriodOver = 0;
while choicePeriodOver == 0
    %% check time
    timeNow = GetSecs;
    if timeNow > (onsetDispChoiceOptions + t_choice)
        % finish the trial
        choicePeriodOver = 1;
        choice_trial = 0;
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