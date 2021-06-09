function[onsets] = choice_and_perf_trainingInstructions(scr, R_or_P_or_RP_condition, t_instructions)
% [onsets] = choice_and_perf_trainingInstructions(scr, R_or_P_or_RP_condition, t_instructions)
% choice_and_perf_trainingInstructions will display instructions before
% starting training for choice and performance.
%
% INPUTS
% scr: structure with screen parameters
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
white = scr.colours.white;
black = scr.colours.black;
wrapat = scr.wrapat;
yScreenCenter = scr.yCenter;

%% instruction that main task will start soon
for iTimeLoop = 1:2
    switch R_or_P_or_RP_condition
        case 'R'
            DrawFormattedText(window,...
                ['Vous allez a present choisir entre deux options associees a differents niveaux de recompense et d''effort '...
                'l''option qui vous parait la plus interessante.'],...
                'center', yScreenCenter/3, white, wrapat);
        case 'P'
            DrawFormattedText(window,...
                ['Vous allez a present choisir entre deux options associees a differents niveaux de pertes et d''effort '...
                'l''option qui vous parait la moins penible.'],...
                'center', yScreenCenter/3, white, wrapat);
        case 'RP'
            DrawFormattedText(window,...
                ['Vous allez a present choisir entre deux options associees a differents niveaux de recompenses ou de pertes et d''effort '...
                'l''option qui vous parait preferable.'],...
                'center', yScreenCenter/3, white, wrapat);
    end
    if iTimeLoop == 1 % force them to read at first
        [~, onsets.trainingWillStart] = Screen(window, 'Flip');
        WaitSecs(t_instructions);
    elseif iTimeLoop == 2 % after t_instructions seconds, they can manually start
        DrawFormattedText(window, 'Appuyez quand vous etes pret(e) a commencer la tache.',...
            'center', yScreenCenter*15/8, white);
        [~, onsets.trainingWillStart_bis] = Screen(window, 'Flip');
        KbWait;
    end
end % loop over forced reading/manual pass loop

end % function