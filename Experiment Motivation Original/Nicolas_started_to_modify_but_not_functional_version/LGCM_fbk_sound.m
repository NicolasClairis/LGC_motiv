function [sound] = LGCM_fbk_sound(main_task_folder)
% [sound] = LGCM_fbk_sound(main_task_folder)
% LGCM_fbk_sound will create audio output for win and loose feedback in the
% version where participants get an audio feedback
%
% INPUTS
% main_task_folder: path where main script lies along with the Sounds
% folder
%
% OUTPUTS
% sound: structure
%
% See also LGCM_main_experiment.m


%% directory where sounds lie
sounds_folder = [main_task_folder, 'Sounds', filesep];

%% Initialize sound driver and prepare 2 audio input we will use. Resample the second one
InitializePsychSound(1);
[sound.audio_win,sound.Fs_win] = audioread([sounds_folder 'Win.mp3']);
[sound.audio_lose,sound.Fs_lose] = audioread([sounds_folder 'loose.mp3']);
[p,q] = rat(sound.Fs_win / sound.Fs_lose);
sound.audio_lose = resample(sound.audio_lose, p, q);
sound.numberChannels = 2;

%% Open Psych-Audio port, with the following arguements
% (1) [] = default sound device
% (2) 1 = sound playback only
% (3) 1 = default level of latency
% (4) Requested frequency in samples per second
% (5) 2 = stereo putput
sound.pahandle = PsychPortAudio('Open', [], 1, 1, sound.Fs_win, sound.numberChannels);
% sound.pahandle_lose = PsychPortAudio('Open', [], 1, 1, sound.Fs_lose, sound.nrchannels);

%% Set the volume to half
PsychPortAudio('Volume', sound.pahandle, 0.5);

end % function