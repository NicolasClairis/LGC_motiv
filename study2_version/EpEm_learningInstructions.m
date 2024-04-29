function[onsets] = EpEm_learningInstructions(scr, stim, t_instructions)
% [onsets] = EpEm_learningInstructions(scr, stim, t_instructions)
% EpEm_learningInstructions will display instructions before
% starting learning to associate each level of size of the effort circle to
% a given level of effort.
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with stimuli informations
%
% t_instructions: time to display training instructions before allowing
% participant to press button and start the training
%
% OUTPUTS
% onsets: onsets when training starts

%% load main paramaters
window = scr.window;
wrapat = scr.wrapat;

%% instruction that main task will start soon
for iTimeLoop = 1:2
    DrawFormattedText(window,...
                stim.EpEm_learning.text,...
                stim.EpEm_learning.x,...
                stim.EpEm_learning.y,...
                stim.EpEm_learning.colour, wrapat);
    if iTimeLoop == 1 % force them to read at first
        [~, onsets.learningEffortLevelsWillStart] = Screen(window, 'Flip');
        WaitSecs(t_instructions);
    elseif iTimeLoop == 2 % after t_instructions seconds, they can manually start
        % display text: Press when you are ready to start
        DrawFormattedText(window, stim.pressWhenReady.text,...
            stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour);
        [~, onsets.learningEffortLevelsWillStart_bis] = Screen(window, 'Flip');
        KbQueueWait(0,3);
    end
end % loop over forced reading/manual pass loop

end % function