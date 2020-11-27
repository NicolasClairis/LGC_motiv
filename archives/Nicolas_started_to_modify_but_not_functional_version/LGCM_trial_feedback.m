function[onsets_fbk_trial] = LGCM_trial_feedback(scr,stim,speed,totalMoney,doWin,incentive, reward_or_punishment)
%[onsets_fbk_trial] = LGCM_trial_feedback(scr,stim,speed,totalMoney,doWin,incentive)
%LGCM_trial_feedback will show feedback of the current trial in terms of
%total amount of money accumulated since the beginning of the session and
%of trial performance
%
% See also LGCM_main_experiment.m


%% load general variables
% screen parameters
window = scr.window;
screenYpixels = scr.Ycenter*2;
% text size
big_textSize = scr.textSize.big;
middle_textSize = scr.textSize.middle;

%% show feedback
switch reward_or_punishment
    case 'R'
        
        switch doWin
            case 1 % gain
                Screen('TextSize', window, big_textSize);
                DrawFormattedText(window, ['Gain: ',num2str(amountMoney),' Fr'], 'center',...
                    screenYpixels*0.3, [0 0.8 0 ]);
            case 0
        end
    case 'P'
        
end % reward/punishment

%% Draw text total amount in the bottom of the screen
Screen('TextSize', window, middle_textSize);
DrawFormattedText(window, strcat('Montant total:\n ',num2str(totalMoney),' Fr'), 'center',...
    screenYpixels*0.7, 1);

%% flip
% speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
[~,onsets_fbk_trial] = Screen('Flip',window);
WaitSecs(t_fbk);


end % function