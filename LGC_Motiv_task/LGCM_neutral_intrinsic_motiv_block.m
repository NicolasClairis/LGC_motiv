function [all] = LGCM_neutral_intrinsic_motiv_block(all, scr, stim, speed, sound, totalMoney, nbTrialPerCoinType)
%

%% Launch a final short break before Neutral block
speed.isShortBreak = true;
stim.isPositiveStimuli = isPositiveStimuli(3);
showBreak(scr,stim,speed,totalMoney)

%% Launch Neutral block
for iTrial = 1:nbTrialPerCoinType
    stim.imageTextureIdx = 4;
    [doWin,signal,firstT2] = stimulusPresentation(scr,stim,speed,sound);
    all.signals{4,iTrial} = signal;
    all.win{4,iTrial}  = doWin;
    all.VCmax{4,iTrial} = max(signal);
    all.trialLength{4,iTrial} = length(signal);
    all.incentive{4,iTrial} = 0;
    all.firstT2{4,iTrial} = firstT2;
    % As there is no result. Wait the same amount of time before launching it again
    pause(3)
end

end % function