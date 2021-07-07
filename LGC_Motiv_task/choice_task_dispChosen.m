function[time_dispChoice] = choice_task_dispChosen(scr, stim, R_chosen, E_chosen, R_or_P)
% [time_dispChoice] = choice_task_dispChosen(scr, stim, choice_opt, choice,...
%     R_or_P, iTrial)
% choice_task_dispChosen will display the chosen option
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
% See also choice_task_main.m

%% extract relevant parameters
window = scr.window;
white = scr.colours.white;

% remind the option they chose
DrawFormattedText(window, stim.chosenOptionMsg.text,...
    stim.chosenOptionMsg.x, stim.chosenOptionMsg.y, white);

%% display reward and effort level
switch R_chosen
    case 0 % no option was selected
        DrawFormattedText(window,...
            stim.feedback.error_tooSlow.text,...
            stim.feedback.error_tooSlow.x, stim.feedback.error_tooSlow.y, white);
    otherwise % one option was selected
        
        % if punishment trial, add also indication to know that money is to be lost
        switch R_or_P
            case 'R'
                DrawFormattedText(window,stim.choice.win.text,...
                    stim.winRewardText.top_center(1),...
                    stim.winRewardText.top_center(2),...
                    white);
            case 'P'
                DrawFormattedText(window,stim.choice.lose.text,...
                    stim.loseRewardText.top_center(1),...
                    stim.loseRewardText.top_center(2),...
                    white);
        end
        drawRewardAmount(scr, stim, R_chosen, R_or_P, 'top_center_start');

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
        DrawFormattedText(window, stim.choice.for.text,...
            stim.effort_introText.bottom_center(1),...
            stim.effort_introText.bottom_center(2),...
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