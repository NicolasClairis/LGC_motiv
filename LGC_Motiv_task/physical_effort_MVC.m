function[MVC, onsets] = MVC_measurement(scr, dq, n_MVC_repeat, calibTimes)
%[MVC, onsets] = MVC_measurement(scr, dq, n_MVC_repeat, calibTimes)
% MVC_measurement will display the instructions and measure the MVC
% for the physical effort task.
%
% INPUTS
% scr: structure about main screen parameters (size, window, etc.)
%
% stim: structure with stimulus properties (size, etc.)
%
% dq: device from which the force will be recorded
%
% n_MVC_repeat: number of repetitions of the MVC measurement
%
% calibTimes: structure with timing information
%   .instructions: instructions duration
%   .effort_max: duration available for the effort performance
%   .fbk: duration for the feedback
%
% OUTPUTS
% MVC: MVC value
%
% onsets: structure with the onset value of each step of this measure
%
% See also choice_task_main.m

%% screen relevant variables
window = scr.window;
xScreenCenter = scr.xCenter;
yScreenCenter = scr.yCenter;
yScreenSize = yScreenCenter*2;
text_size_1 = 100;
text_size_2 = 150;
text_size_3 = 140;

%% extract timings
t_MVC_calib_instructions1 = calibTimes.instructions_bis;
t_MVC_calib = calibTimes.effort_max;
t_MVC_rest = calibTimes.t_MVC_rest;

%% MVC
maxVoltage = 10; % maximum voltage that can be reached by the grip (force will be normalized to this value)
MVC_perCalibSession = NaN(1,n_MVC_repeat);

%% Quick text to introduce MVC calibration
Screen('TextSize', window, text_size_1);
DrawFormattedText(window, ['Avant de commencer l''expérience, ',...
    'nous allons vous demander',...
    'de serrer la poignée de force au maximum de vos capacités plusieurs',...
    'fois d''affilée.'], 'center',yScreenSize*0.7, 1);
Screen('TextSize', window, text_size_2)
DrawFormattedText(window, 'Tenez-vous prêt à serrer', 'center', yScreenSize*0.3, 1);

[~,time_disp1,~,~,~] = Screen(window,'Flip');
onsets.initial_MVC_instructions = time_disp1;
WaitSecs(t_MVC_calib_instructions1);

%% Measure MVC
% Set screen text size
Screen('TextSize', window, text_size_3);

%% initialize onsets
[onsets.effortScale_start,...
    onsets.initial_MVC_rest] = deal(NaN(1,n_MVC_repeat));

%% Measure MVC and keep maximal value
for iCalib_MVC = 1:n_MVC_repeat
    %% store all force values
    live_force_MVC.(['calib_session_',num2str(iCalib_MVC)]) = [];
    
%     %% allow subject to decide when to start calibration
%     DrawFormattedText(window, 'Appuyez sur un des boutons quand vous êtes prêt(e)',...
%         'center', yScreenSize*0.8, [0 0.8 0 ]);
%     [~,timeNow]  = Screen(window,'Flip');
%     onset.initial_MVC_rest = timeNow;
% you need to add key in the inputs if you use these lines of code
%     [keyisdown, ~, keyCode] = KbCheck();
%     while (keyisdown == 0) ||...
%             (keyCode(key.left) == 0 && keyCode(key.right) == 0)
%         [keyisdown, ~, keyCode] = KbCheck();
%     end
    
    %% start displaying effort scale and Go signal
    DrawFormattedText(window, 'GO !', 'center', yScreenSize*0.8, 1);
    F_start = 0;
    F_threshold = 100;
    F_tolerance = 0;
    disp_realtime_force(scr, F_threshold, F_tolerance, F_start, 'calib');
    [~,timeEffortScaleStart]  = Screen(window,'Flip');
    onsets.effortScale_start(i_Calib_MVC) = timeEffortScaleStart;
    
    %% During t_MVC_calib second, show signal power and record it
    forceCalib.(['calibTrial_',num2str(iCalib_MVC)]) = [];
    timeNow = GetSecs;
    while timeNow < time_disp0 + t_MVC_calib
        timeNow = GetSecs;
        F_now_Voltage_table = read(dq);
        F_now_Voltage = F_now_Voltage_table.Variables;
        % convert in percentage of maximal voltage
        F_now = F_now_Voltage/maxVoltage;
        DrawFormattedText(window, 'GO !', 'center', yScreenSize*0.8, 1);
        disp_realtime_force(scr, F_threshold, F_tolerance, F_now, 'calib');
        Screen(window,'Flip');
        forceCalib.(['calibTrial_',num2str(iCalib_MVC)]) = [forceCalib.(['calibTrial_',num2str(iCalib_MVC)]);...
            [F_now, timeNow, F_now_Voltage]];
    end % time for the current calibration trial
    
    %% Show a rest text and give some rest
    DrawFormattedText(window, 'Reposez-vous quelques secondes', 'center', yScreenSize*0.8, [0 0.8 0 ]);
    [~,timeNow]  = Screen(window,'Flip');
    onsets.initial_MVC_rest(iCalib_MVC) = timeNow;
    WaitSecs(t_MVC_rest);
    
    %% extract max force for this session
    MVC_perCalibSession(iCalib_MVC) = nanmax(live_force_MVC.(['calib_session_',num2str(iCalib_MVC)]));
    
end % calibration loop

%% store max MVC measure in output
MVC.MVC_perCalibSession = MVC_perCalibSession;
MVC.MVC = nanmax(MVC_perCalibSession);

end % function