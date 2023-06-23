function [timenow1] = totalGainDisp(window, totalGain)
%totalGainDisp function common to all three taskLearning75, taskGripRP and
% taskMentalRP scripts to display final gain for this run on the screen
%
% INPUTS
% window: window where will be displayed
% totalGain: numeric value with amount won in this run
%
% made by Nicolas Clairis - february 2017

Screen(window,'TextSize',40);
wrapat = 48; % go to line if text too long
vspacing = 2; % higher space between lines

if totalGain > 1
    DrawFormattedText(window,['Pour cette tache, vous avez gagne ',num2str(totalGain),' euros.'],'center','center', [255 255 255],wrapat, [], [], vspacing);
elseif totalGain == 1
    DrawFormattedText(window,'Pour cette tache, vous avez gagne 1 euro.','center','center', [255 255 255],wrapat, [], [], vspacing);
else
    DrawFormattedText(window,'Pour cette tache, vous n''avez pas gagne d''argent.','center','center', [255 255 255],wrapat, [], [], vspacing);
end

[~,timenow1,~,~,~] = Screen(window, 'Flip');

end

