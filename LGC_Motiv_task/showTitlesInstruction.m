function showTitlesInstruction(scr, stim, instructionType, effortTypeLetter, keys)
%showTitlesInstruction(scr, stim, instructionType, effortTypeLetter, keys)
% Function to display information of the next things to do.
%
% INPUTS
% scr = screen informations
% 
% stim: structure with stimuli informations
%
% instructionType = what is the next step
% "Learning" participants are understanding how to do the effort
% "training" is when you perform the full trial and the results do not matter
% "task" is when you launch a block of trials, the core of the experiment
%
% effortTypeLetter = define if we are working on physical ('p') or mental ('m')
% effort.
%
% langage: 'fr' for french; 'engl' for english
%
% keys: structure with information of code of keys to check

%% initialize relevant parameters
window = scr.window;
titleTextSize = scr.textSize.taskPeriodsTitles;
baselineTextSize = scr.textSize.baseline;

%% define title settings

% change text size
Screen('TextSize', window, titleTextSize);
% define effort type
effortType = ['E',effortTypeLetter];

% announce what is next
DrawFormattedText(window, stim.(effortType).(instructionType).title.text,...
            stim.(effortType).(instructionType).title.x,...
            stim.(effortType).(instructionType).title.y,...
            stim.(effortType).(instructionType).title.colour);

% flip information on the screen
waitSpace(scr, stim, window, keys);

%% put back baseline textsize value
Screen('TextSize', window, baselineTextSize);

end % function