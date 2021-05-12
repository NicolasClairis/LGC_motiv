function[MVC, onsets] = physical_effort_MVC(scr, dq, n_MVC_repeat, calibTimes)
%[MVC, onsets] = physical_effort_MVC(scr, dq, n_MVC_repeat, calibTimes)
% MVC_measurement will display the instructions and measure the MVC
% for the physical effort task.
%
% INPUTS
% scr: structure about main screen parameters (size, window, etc.)
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
% text_size_1 = 50;
% text_size_2 = 75;
% text_size_3 = 70;
orange = [255 153 0];
bottomScaleLimit    = yScreenCenter*(3/2); % bottom limit of the scale
topScaleLimit       = yScreenCenter*(1/2); % upper limit of the scale
leftScaleLimit      = xScreenCenter*(3.5/4); % left limit of the scale
rightScaleLimit     = xScreenCenter*(4.5/4); % right limit of the scale
graphYSize = bottomScaleLimit - topScaleLimit;

%% force relevant variables
F_start = 0; % initial force level at zero
F_threshold = 100; % force threshold on top to incentivize the participants
F_tolerance = 0.25; % will determine the width of the top rectangle

%% start acquiring the data in the background (if you don't use this
% function, everytime you call the read function, it will take a
% long time to process)
start(dq,"continuous");
% will need data = read(dq) function only to read the signal


%% extract timings
t_MVC_calib_instructions1 = calibTimes.instructions_bis;
t_MVC_calib = calibTimes.effort_max;
t_MVC_rest = calibTimes.MVC_rest;
t_readWait = calibTimes.physicalReadWait;

%% MVC
maxVoltage = 10; % maximum voltage that can be reached by the grip (force will be normalized to this value)
MVC_perCalibSession = NaN(1,n_MVC_repeat);

%% Quick text to introduce MVC calibration
% Screen('TextSize', window, text_size_1);
DrawFormattedText(window, ['Avant de commencer l''expérience, ',...
    'nous allons vous demander ',...
    'de serrer la poignée de force au maximum de vos capacités plusieurs ',...
    'fois d''affilée.'], 'center',yScreenSize*0.7, 1, scr.wrapat);
% Screen('TextSize', window, text_size_2)
DrawFormattedText(window, 'Tenez-vous prêt à serrer la poignée.', 'center', yScreenSize*0.3, 1);

[~,time_disp1,~,~,~] = Screen(window,'Flip');
onsets.initial_MVC_instructions = time_disp1;
WaitSecs(t_MVC_calib_instructions1);

%% Measure MVC
% Set screen text size
% Screen('TextSize', window, text_size_3);

%% initialize onsets
[onsets.effortScale_start,...
    onsets.initial_MVC_rest] = deal(NaN(1,n_MVC_repeat));

%% initialize Force variable
F_now_Voltage_tmp = read(dq,'all','OutputFormat','Matrix');
pause(t_readWait);
if ~isempty(F_now_Voltage_tmp)
    F_now_Voltage = F_now_Voltage_tmp(end);
    % convert force level from Voltage to a percentage of MVC
    F_now = (F_now_Voltage/maxVoltage)*100;
else % record when the output of read was empty to know when the force level was kept equal because of read failure
    F_now = 0;
end

%% Measure MVC and keep maximal value
for iCalib_MVC = 1:n_MVC_repeat
    
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
    disp_realtime_force(scr, F_threshold, F_tolerance, F_start, 'calib');
    [~,timeEffortScaleStart]  = Screen(window,'Flip');
    onsets.effortScale_start(iCalib_MVC) = timeEffortScaleStart;
    
    %% During t_MVC_calib second, show signal power and record it
    forceCalib.(['calibTrial_',num2str(iCalib_MVC)]) = [];
    timeNow = GetSecs;
    while timeNow < timeEffortScaleStart + t_MVC_calib
        timeNow = GetSecs;
        F_now_Voltage_tmp = read(dq,'all','OutputFormat','Matrix');
        pause(t_readWait);
        % convert in percentage of maximal voltage
        if ~isempty(F_now_Voltage_tmp)
            F_now_Voltage = F_now_Voltage_tmp(end);
            % convert force level from Voltage to a percentage of MVC
            F_now = (F_now_Voltage/maxVoltage)*100;
            sampleOk_tmp = 1;
        else % record when the output of read was empty to know when the force level was kept equal because of read failure
            sampleOk_tmp = 0;
        end
        % store force levels in the output
        forceCalib.(['calibTrial_',num2str(iCalib_MVC)]) = [forceCalib.(['calibTrial_',num2str(iCalib_MVC)]);...
            [F_now, timeNow, F_now_Voltage, sampleOk_tmp]]; % store F in % of MVC, time and F in Volts
        DrawFormattedText(window, 'GO !', 'center', yScreenSize*0.9, 1);
        disp_realtime_force(scr, F_threshold, F_tolerance, F_now, 'calib');
        
        % for calibration trials coming after the first one, you can also
        % display the max reached until now to incentivize them to make
        % better?
        if iCalib_MVC > 1
            maxMVCuntilNow = nanmax( MVC_perCalibSession(1:(iCalib_MVC-1)));
            yThreshold = bottomScaleLimit - graphYSize*(maxMVCuntilNow/maxVoltage);
            Screen('DrawLine', window, orange, leftScaleLimit, yThreshold, rightScaleLimit, yThreshold,5);
        end
        
        Screen(window,'Flip');
        %         [lastFrameTime, timeDispNow]  = Screen('Flip', window, lastFrameTime + (0.5*ifi));
    end % time for the current calibration trial
    
    %% Show a rest text and give some rest
    DrawFormattedText(window, 'Reposez-vous quelques secondes.', 'center', yScreenSize*0.8, [0 0.8 0 ]);
    [~,timeNow]  = Screen(window,'Flip');
    onsets.initial_MVC_rest(iCalib_MVC) = timeNow;
    WaitSecs(t_MVC_rest);
    
    %% extract max force for this session (in Voltage)
    MVC_perCalibSession(iCalib_MVC) = nanmax(forceCalib.(['calibTrial_',num2str(iCalib_MVC)])(:,3));
    
end % calibration loop

%% stop acquisition of biopac handgrip
stop(dq);

%% store max MVC measure in output
MVC.forceCalib = forceCalib;
MVC.MVC_perCalibSession = MVC_perCalibSession;
MVC.MVC = nanmax(MVC_perCalibSession);

end % function