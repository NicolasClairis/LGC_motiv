function [all] = LGCM_neutral_intrinsic_motiv_block(all,...
    scr, stim, speed, sound, totalMoney,...
    pause_dur,...
    nbTrialPerCoinType)
%[all] = LGCM_neutral_intrinsic_motiv_block(all,...
%     scr, stim, speed, sound, totalMoney,...
%     pause_dur,...
%     nbTrialPerCoinType)
% LGCM_neutral_intrinsic_motiv_block allows to perform the neutral block in
% case you want to study intrinsic motivation.
%
% INPUTS
% all: structure with all behavioral relevant informations
% scr:
% stim:
% speed:
% sound:
% totalMoney
% pause_dur: duration of the break defined in the main script
% nbTrialPerCoinType: number of repetitions of each mini-block
%
% OUTPUTS
% all: structure with all the relevant information updated after performing
% the neutral block
%
% See also LGC_main_experiment.m & LGCM_stimulusPresentation.m

%% Launch a final short break before Neutral block
speed.isShortBreak = true;
stim.isPositiveStimuli = isPositiveStimuli(pause_dur);
showBreak(scr,stim,speed,totalMoney)

%% Launch Neutral block
for iTrial = 1:nbTrialPerCoinType
    stim.imageTextureIdx = 4;
    [doWin,signal,firstT2] = LGCM_stimulusPresentation(scr,stim,speed,sound);
    all.signals{4,iTrial} = signal;
    all.win{4,iTrial}  = doWin;
    all.VCmax{4,iTrial} = max(signal);
    all.trialLength{4,iTrial} = length(signal);
    all.incentive{4,iTrial} = 0;
    all.firstT2{4,iTrial} = firstT2;
    
    % As there is no result. Wait the same amount of time before launching it again
    pause(pause_dur)
end

end % function