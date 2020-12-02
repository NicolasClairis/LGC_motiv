function[time_dispChoice, R_chosen, E_chosen] = LGCM_choice_task_dispChosen(scr, stim, choice, n_E_levels,...
    E_left, E_right,...
    R_left, R_right)
% [time_dispChoice, R_chosen, E_chosen] = LGCM_choice_task_dispChosen(scr, stim, choice, n_E_levels,...
%     E_left, E_right,...
%     R_left, R_right)
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with stimuli parameters (reward and effort display
% informations are stored in here)
%
% choice
% (-1) left option => left reward and effort
% (0) no option chosen => show higher effort with no associated reward
% (1) right option => right reward and effort
%
% n_E_levels: number of difficulty levels
% 
% E_left, E_right: levels of difficulty for left and right options
%
% R_left, R_right: levels of reward for left and right options
%
% OUTPUTS
% time_dispChoice: onset of the display of the chosen option on the screen
%
% See also LGCM_choice_task_main.m

%% extract relevant parameters
window = scr.window;
% extract info for max possible effort
stimDiffMaxColor = stim.difficulty.maxColor;
stimMaxDiff = stim.difficulty.middle.(['level_',num2str(n_E_levels)]);
stimDiffOvalWidth = stim.difficulty.ovalWidth;
% extract info for left and right stim difficulty and reward
diffOvalColor = stim.difficulty.currLevelColor;

%% display maximal difficulty level
Screen('FrameOval', window, stimDiffMaxColor,...
    stimMaxDiff,...
    stimDiffOvalWidth);

%% extract reward and difficulty level of the chosen option
switch choice
    %% left chosen
    case -1
        E_chosen = E_left;
        R_chosen = R_left;
    %% right chosen
    case 1
        E_chosen = E_right;
        R_chosen = R_right;
    %% by default opt for the higher effort with no reward if no option was selected
    case 0
        E_chosen = max(E_left, E_right);
        R_chosen = 0;
end % choice

%% display difficulty level
Screen('FrameOval', window, diffOvalColor,...
    stim.difficulty.middle.(['level_',num2str(E_chosen)]),...
    stimDiffOvalWidth); % left option difficulty

%% display associated reward
Screen('DrawTexture', window,...
    stim.reward.texture.(['reward_',num2str(R_chosen)]),...
    [],...
    stim.reward.middle.(['reward_',num2str(R_chosen)]));

%% display on screen and extract timing
[~,time_dispChoice] = Screen('Flip',window);

end % function