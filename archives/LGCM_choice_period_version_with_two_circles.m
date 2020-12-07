function[choice_trial, onsetDispChoiceOptions, onsetChoice] = LGCM_choice_period(scr, stim,...
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

%% initialize variables of interest
choice_trial = 0;
onsetChoice = NaN;
window = scr.window;

%% display maximal difficulty level for each option
Screen('FrameOval', window, stim.difficulty.maxColor,...
    stim.difficulty.left.(['level_',num2str(n_E_levels)]),...
    stim.difficulty.ovalWidth);
Screen('FrameOval', window, stim.difficulty.maxColor,...
    stim.difficulty.right.(['level_',num2str(n_E_levels)]),...
    stim.difficulty.ovalWidth);
% extract difficulty level for each side of the screen
E_left_tmp  = 1; % E_list.left(iTrial)
E_right_tmp = 2; % E_list.right(iTrial)
% extract reward level for each side of the screen
R_left_tmp  = 1; % R_list.left(iTrial)
R_right_tmp = 3; % R_list.right(iTrial)

%% display each difficulty level
Screen('FrameOval', window, stim.difficulty.currLevelColor,...
    stim.difficulty.left.(['level_',num2str(E_left_tmp)]),...
    stim.difficulty.ovalWidth); % left option difficulty
Screen('FrameOval', window, stim.difficulty.currLevelColor,...
    stim.difficulty.right.(['level_',num2str(E_right_tmp)]),...
    stim.difficulty.ovalWidth); % right option difficulty

%% display each reward level
Screen('DrawTexture', window,...
    stim.reward.texture.(['reward_',num2str(R_left_tmp)]),...
    [],...
    stim.reward.left.(['reward_',num2str(R_left_tmp)]));
Screen('DrawTexture', window,...
    stim.reward.texture.(['reward_',num2str(R_right_tmp)]),...
    [],...
    stim.reward.right.(['reward_',num2str(R_right_tmp)]));

[~,onsetDispChoiceOptions] = Screen('Flip',window);

%% wait for choice to be made or time limit to be reached
timeNow = onsetDispChoiceOptions;
while timeNow <= onsetDispChoiceOptions + t_choice
    timeNow = GetSecs;
    [keyisdown, secs, keycode] = KbCheck;
    
    %% left option chosen
    if keyisdown == 1 &&...
            keycode(key.left) == 1 && keycode(key.right) == 0
        timedown = secs;
        choice_trial = -1;
    %% right option chosen
    elseif keyisdown == 1 &&...
            keycode(key.left) == 0 && keycode(key.right) == 1
        timedown = secs;
        choice_trial = 1;
    end
    %% get time
    onsetChoice = timedown;
    warning('here you can either leave it that way (no visual feedback at the moment they chose an option',...
        ' or show some visual feedback of the selected option already',...
        ' or require the subject to keep pressing the button until the end of the trial.');
end % choice period

end % function