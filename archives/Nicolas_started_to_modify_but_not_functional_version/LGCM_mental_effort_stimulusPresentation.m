function[doWin, Stroop_perf, onsets] = LGCM_mental_effort_stimulusPresentation(scr, audio_fbk_yn, sound)
%[doWin, Stroop_perf, onsets] = LGCM_mental_effort_stimulusPresentation(scr, audio_fbk_yn, sound)
%LGCM_mental_effort_stimulusPresentation will show the incentive at play,
%then ask to perform as many numerical stroop pairs as necessary
%
% INPUTS
% scr: structure with main screen parameters
%
% audio_fbk_yn:
% 'yes': include audio feedback
% 'no': no audio feedback
%
% sound: structure with sound to use for win/lose
%
% OUTPUTS
% doWin: did the subject win the trial (1) or did he miss it (0)?
%
% Stroop_perf: detail about the nature of the stroop pairs performed and
% whether the trial was successfull or not
%
% onsets: structure with timings
%

%% extract main variables
% screen
window = scr.window;
xCenter = scr.xCenter;
yCenter = scr.yCenter;
% stimulus
signalColor = stim.signalColor;
threshold_1 = stim.threshold_1;
threshold_2 = stim.threshold_2;
% timings
t_effortPhase = t_wait.effortTotalPhase; % duration of the effort period


%% 

end % function