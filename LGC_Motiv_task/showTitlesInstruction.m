function showTitlesInstruction(scr,instructionType,isMental)
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

t_wait = 3;
Screen('TextSize', scr.window, 80);
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
            'center',yScreenCenter/2,blackCol, wrapat);
        
    case 'Training'
        DrawFormattedText(window,...
            ['Entrainement ',effortType],...
            'center',yScreenCenter/2,blackCol, wrapat);
    case 'Task'
        DrawFormattedText(window,...
            ['Tâche ',effortType],...
            'center',yScreenCenter/2,blackCol, wrapat);
end
%flip information on the screen
Screen(scr.window,'Flip');
WaitSecs(t_wait);
% put back baseline textsize value
Screen('TextSize', window, 40);
end

