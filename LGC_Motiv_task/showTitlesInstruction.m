function showTitlesInstruction(scr, instructionType, isMental)
% Function to display information of the next things to do
%
% Scr = screen informations
% instructionType = what is the next step
% "Learning" participants are understanding how to do the effort
% "training" is when you perform the full trial and the results do not matter
% "task" is when you launch a block of trials, the core of the experiment
% isMental = define if we are working on physical or mental effort
%
%
%

%% initialize relevant parameters
window = scr.window;
yScreenCenter = scr.yCenter;
baselineTextSize = scr.textSize.baseline;

%% define title settings
titleSize = 80;
titleCol = scr.colours.white;

t_wait = 3;
Screen('TextSize', window, titleSize);
if isMental == true
    effortType = 'Mental';
elseif isMental == false
    effortType = 'Physique';
end

% Announce what is next
switch instructionType
    case 'learning'
        DrawFormattedText(window,...
            ['Apprentissage ',effortType],...
            'center',yScreenCenter/2,titleCol);
        
    case 'training'
        DrawFormattedText(window,...
            ['Entrainement ',effortType],...
            'center',yScreenCenter/2,titleCol);
    case 'task'
        switch effortType
            case 'Physique'
                effortType_bis = effortType;
            case 'Mental'
                effortType_bis = 'Mentale';
        end
        DrawFormattedText(window,...
            ['Tache ',effortType_bis],...
            'center',yScreenCenter/2,titleCol);
end
%flip information on the screen
Screen(window,'Flip');
WaitSecs(t_wait);
% put back baseline textsize value
Screen('TextSize', window, baselineTextSize);
end

