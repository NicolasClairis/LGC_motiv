function[time_dispChoice, R_chosen, E_chosen] = LGCM_choice_task_dispChosen(scr, stim, choice_opt, choice,...
    R_or_P, iTrial)
% [time_dispChoice, R_chosen, E_chosen] = LGCM_choice_task_dispChosen(scr, stim, choice_opt, choice,...
%     R_or_P, iTrial)
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
% R_or_P: character indicating the nature of the current trial
% 'R': reward trial
% 'P': punishment trial
%
% iTrial: trial number
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
% See also LGCM_choice_task_main.m

%% extract relevant parameters
window = scr.window;
xScreenCenter = scr.xCenter;
yScreenCenter = scr.yCenter;
white = scr.colours.white;
black = scr.colours.black;

% remind the option they chose
DrawFormattedText(window,'Vous avez choisi','center',yScreenCenter/6.5,white);

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

%% display reward and effort level
switch R_chosen
    case 0 % no option was selected
        DrawFormattedText(window,...
            'Trop lent!',...
            'center', yScreenCenter/2, white);
    otherwise % one option was selected
        % reward level
        R_chosen_nm = ['reward_',num2str(R_chosen)];
        Screen('DrawTexture', window,...
            stim.reward.texture.(R_chosen_nm),...
            [],...
            stim.chosenOption.reward.(R_chosen_nm));
        
        % if punishment trial, add also indication to know that money is to be lost
        switch R_or_P
            case 'R'
                % define coordinates
                xStart_R_txt = xScreenCenter - stim.textRectSize.xSizeWin/2;
                yStart_R_txt = stim.reward.top_center.(R_chosen_nm)(2)-stim.textRectSize.ySizeWin/3;
                % display
                DrawFormattedText(window,'Gagner',...
                    xStart_R_txt,...
                    yStart_R_txt,...
                    white);
            case 'P'
            % version with black cross on top of the monetary incentives
%             lineWidth = 10;
%             % cross on monetary incentive
%             Screen('DrawLine', window, black,...
%                 stim.reward.top_center.(R_chosen_nm)(1), stim.reward.top_center.(R_chosen_nm)(2),...
%                 stim.reward.top_center.(R_chosen_nm)(3), stim.reward.top_center.(R_chosen_nm)(4),...
%                 lineWidth);
%             Screen('DrawLine', window, black,...
%                 stim.reward.top_center.(R_chosen_nm)(3), stim.reward.top_left.(R_chosen_nm)(2),...
%                 stim.reward.top_center.(R_chosen_nm)(1), stim.reward.top_center.(R_chosen_nm)(4),...
%                 lineWidth);

            % version with negative overlay on top of monetary incentive
            Screen('FillOval', window, stim.punishment.colourOverlay,...
                    stim.punishment.circleOverlay.top_center.(R_chosen_nm));
                
                % define coordinates
                xStart_R_txt = xScreenCenter - stim.textRectSize.xSizeLose/2;
                yStart_R_txt = stim.reward.top_center.(R_chosen_nm)(2)-stim.textRectSize.ySizeLose/3;
                % display
                DrawFormattedText(window,'Perdre',...
                    xStart_R_txt,...
                    yStart_R_txt,...
                    white);
        end
        
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
        
        % define coordinates
        xStart_forE_txt = xScreenCenter - stim.textRectSize.xSizeForEffort/2;
        yStart_forEffort_txt = stim.difficulty.below_center(2)-stim.textRectSize.ySizeForEffort/2;
        % display
        DrawFormattedText(window,'pour',...
            xStart_forE_txt,...
            yStart_forEffort_txt,...
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