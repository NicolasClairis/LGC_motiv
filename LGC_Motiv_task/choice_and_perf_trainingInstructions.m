function[onsets] = choice_and_perf_trainingInstructions(scr, stim, R_or_P_or_RP_condition, t_instructions)
% [onsets] = choice_and_perf_trainingInstructions(scr, stim, R_or_P_or_RP_condition, t_instructions)
% choice_and_perf_trainingInstructions will display instructions before
% starting training for choice and performance.
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with stimuli informations
%
% R_or_P_or_RP_condition: string indicating training condition
% 'R': pure reward training
% 'P': pure punishment training
% 'RP': mixed trials with reward and punishments
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
                stim.training.(R_or_P_or_RP_condition).text,...
                stim.training.(R_or_P_or_RP_condition).x,...
                stim.training.(R_or_P_or_RP_condition).y,...
                stim.training.(R_or_P_or_RP_condition).colour, wrapat);
    if iTimeLoop == 1 % force them to read at first
        [~, onsets.trainingWillStart] = Screen(window, 'Flip');
        WaitSecs(t_instructions);
    elseif iTimeLoop == 2 % after t_instructions seconds, they can manually start
        % display text: Press when you are ready to start
        DrawFormattedText(window, stim.pressWhenReady.text,...
            stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour);
        [~, onsets.trainingWillStart_bis] = Screen(window, 'Flip');
        KbQueueWait(0,3);
    end
end % loop over forced reading/manual pass loop

end % function