function showTitlesInstruction(scr, stim, instructionType, isMental)
%showTitlesInstruction(scr, stim, instructionType, isMental)
% Function to display information of the next things to do
%
% scr = screen informations
% 
% stim: structure with stimuli informations
%
% instructionType = what is the next step
% "Learning" participants are understanding how to do the effort
% "training" is when you perform the full trial and the results do not matter
% "task" is when you launch a block of trials, the core of the experiment
%
% isMental = define if we are working on physical or mental effort
%
%
%

%% initialize relevant parameters
window = scr.window;
titleTextSize = scr.textSize.taskPeriodsTitles;
baselineTextSize = scr.textSize.baseline;

%% define title settings

t_wait = 3;

% change text size
Screen('TextSize', window, titleTextSize);
% define effort type
if isMental == true
    effortType = 'Em';
elseif isMental == false
    effortType = 'Ep';
end

% Announce what is next
DrawFormattedText(window, stim.(effortType).(instructionType).title.text,...
            stim.(effortType).(instructionType).title.x,...
            stim.(effortType).(instructionType).title.y,...
            stim.(effortType).(instructionType).title.colour);

% flip information on the screen
Screen(window,'Flip');
WaitSecs(t_wait);

%% put back baseline textsize value
Screen('TextSize', window, baselineTextSize);

end % function