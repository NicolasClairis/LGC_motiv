function[time_dispChoice, R_chosen, E_chosen] = LGCM_choice_task_dispChosen(scr, stim, choice_opt, choice,...
    iTrial)
% [time_dispChoice, R_chosen, E_chosen] = LGCM_choice_task_dispChosen(scr, stim, choice_opt, choice,...
%     iTrial)
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with stimuli parameters (reward and effort display
% informations are stored in here)
%
% choice_opt: structure with info about choice options (side of each
% option, reward level, effort level, etc.)
%
% choice
% (-1) left option => left reward and effort
% (0) no option chosen => show higher effort with no associated reward
% (1) right option => right reward and effort
%
%
% iTrial: trial number
%
% OUTPUTS
% time_dispChoice: onset of the display of the chosen option on the screen
%
% See also LGCM_choice_task_main.m

%% extract relevant parameters
window = scr.window;

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

%% extract reward and difficulty level of the chosen option
switch choice
    %% left chosen
    case -1
        E_chosen = E_left_tmp;
        R_chosen = R_left_tmp;
    %% right chosen
    case 1
        E_chosen = E_right_tmp;
        R_chosen = R_right_tmp;
    %% by default opt for the higher effort with no reward if no option was selected
    case 0
        E_chosen = max(E_left_tmp, E_right_tmp);
        R_chosen = 0;
end % choice

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

%% display associated reward
switch R_chosen
    case 0 % no option was selected
        Screen('DrawText', window, 'Too slow!',...
            xScreenCenter, yScreenCenter/2,...
            [1 1 1]);
    otherwise % one option was selected
        Screen('DrawTexture', window,...
            stim.reward.texture.(['reward_',num2str(R_chosen)]),...
            [],...
            stim.chosenOption.reward.(['reward_',num2str(R_chosen)]));
end

%% display a square on top of selected reward and effort
Screen('FrameRect', window,...
    stim.chosenOption.squareColour,...
    stim.chosenOption.squareRect,...
    stim.chosenOption.squareWidth);

%% display on screen and extract timing
[~,time_dispChoice] = Screen('Flip',window);

end % function