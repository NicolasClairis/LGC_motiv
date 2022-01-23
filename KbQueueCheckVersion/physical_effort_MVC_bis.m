function[MVC, onsets] = physical_effort_MVC(scr, stim, dq, n_MVC_repeat, calibTimes, MVC_or_maxPerf, keys)
%[MVC, onsets] = physical_effort_MVC(scr, stim, dq, n_MVC_repeat, calibTimes, MVC_or_maxPerf, keys)
% MVC_measurement will display the instructions and measure the MVC
% for the physical effort task.
%
% INPUTS
% scr: structure about main screen parameters (size, window, etc.)
%
% stim: structure with stimuli informations
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
% MVC_or_maxPerf: indicate whether calibration ('MVC') or maximal
% performance ('maxPerf')
%
% keys: structure with code of keys of interest to know how to check for
% space press
% 
% OUTPUTS
% MVC: MVC value
%
% onsets: structure with the onset value of each step of this measure
%
% See also choice_task_main.m

%% screen relevant variables
window = scr.window;
wrapat = scr.wrapat;
orange = scr.colours.orange;

%% force relevant variables
F_start = 0; % initial force level at zero
F_threshold = 100; % force threshold on top to incentivize the participants
F_tolerance = 0.25; % will determine the width of the top rectangle

%% extract timings
t_MVC_calib_instructions1 = calibTimes.instructions_bis;
t_MVC_calib = calibTimes.effort_max;
t_MVC_rest = calibTimes.MVC_rest;
t_readWait = calibTimes.physicalReadWait;

%% MVC
maxVoltage = 10; % maximum voltage that can be reached by the grip (force will be normalized to this value)
MVC_perCalibSession = NaN(1,n_MVC_repeat);

%% Quick text to introduce MVC calibration

switch MVC_or_maxPerf
    case 'MVC'
        for iTime = 1:2
            DrawFormattedText(window, stim.Ep.MVC.instructions.text,...
                stim.Ep.MVC.instructions.x, stim.Ep.MVC.instructions.y,...
                stim.Ep.MVC.instructions.colour, wrapat);
            DrawFormattedText(window, stim.Ep.MVC.instructions_bis.text,...
                stim.Ep.MVC.instructions_bis.x, stim.Ep.MVC.instructions_bis.y,...
                stim.Ep.MVC.instructions_bis.colour);
            if iTime == 1
                [~,time_disp1,~,~,~] = Screen(window,'Flip');
                onsets.initial_MVC_instructions = time_disp1;
                WaitSecs(t_MVC_calib_instructions1);
            elseif iTime == 2
                waitSpace_bis(scr, stim, window, keys);
            end
        end
    case 'maxPerf'
        DrawFormattedText(window, stim.Ep.MVC.instructions.text,...
            stim.Ep.MVC.instructions.x, stim.Ep.MVC.instructions.y,...
            stim.Ep.MVC.instructions.colour, wrapat);
        DrawFormattedText(window, stim.Ep.MVC.instructions_bis.text,...
            stim.Ep.MVC.instructions_bis.x, stim.Ep.MVC.instructions_bis.y,...
            stim.Ep.MVC.instructions_bis.colour);
        [~,time_disp1,~,~,~] = Screen(window,'Flip');
        onsets.initial_MVC_instructions = time_disp1;
        WaitSecs(t_MVC_calib_instructions1);
end
[bottomScaleLimit, topScaleLimit, leftScaleLimit, rightScaleLimit, graphYSize] = disp_realtime_force(scr, F_threshold, F_tolerance, F_start, 'calib');

%% Measure MVC
% Set screen text size
% Screen('TextSize', window, text_size_3);

%% initialize onsets
[onsets.effortScale_start,...
    onsets.initial_MVC_rest] = deal(NaN(1,n_MVC_repeat));

%% Measure MVC and keep maximal value
for iCalib_MVC = 1:n_MVC_repeat
    %% initialize Force variable
    F_now = 0;
    
    %% start acquiring the data in the background (if you don't use this
    % function, everytime you call the read function, it will take a
    % long time to process)
    start(dq,"continuous");
    
    %% start displaying effort scale and Go signal
    DrawFormattedText(window, stim.Ep.MVC.GO.text,...
        stim.Ep.MVC.GO.x, stim.Ep.MVC.GO.y, stim.Ep.MVC.GO.colour);
    disp_realtime_force(scr, F_threshold, F_tolerance, F_start, 'calib');
    % for calibration trials coming after the first one, you can also
    % display the max reached until now to incentivize them to make
    % better?
    if iCalib_MVC > 1
        maxMVCuntilNow = max( MVC_perCalibSession(1:(iCalib_MVC-1)));
        yThreshold = bottomScaleLimit - graphYSize*(maxMVCuntilNow/maxVoltage);
        Screen('DrawLine', window, orange, leftScaleLimit, yThreshold, rightScaleLimit, yThreshold,5);
    end
    
    [~,timeEffortScaleStart]  = Screen(window,'Flip');
    onsets.effortScale_start(iCalib_MVC) = timeEffortScaleStart;
    
    %% During t_MVC_calib second, show signal power and record it
    forceCalib.(['calibTrial_',num2str(iCalib_MVC)]) = [];
    timeNow = GetSecs;
    while timeNow < timeEffortScaleStart + t_MVC_calib
        [F_now_Voltage, timeNow, sampleOk_tmp] = F_read(dq, t_readWait);
        % convert in percentage of maximal voltage
        if sampleOk_tmp == 1
            % convert force level from Voltage to a percentage of MVC
            F_now = (F_now_Voltage/maxVoltage)*100;
        end
        % store force levels in the output
        forceCalib.(['calibTrial_',num2str(iCalib_MVC)]) = [forceCalib.(['calibTrial_',num2str(iCalib_MVC)]);...
            [F_now, timeNow, F_now_Voltage, sampleOk_tmp]]; % store F in % of MVC, time and F in Volts
        DrawFormattedText(window, stim.Ep.MVC.GO.text, stim.Ep.MVC.GO.x, stim.Ep.MVC.GO.y, stim.Ep.MVC.GO.colour);
        disp_realtime_force(scr, F_threshold, F_tolerance, F_now, 'calib');
        
        % for calibration trials coming after the first one, you can also
        % display the max reached until now to incentivize them to make
        % better?
        if iCalib_MVC > 1
            maxMVCuntilNow = max( MVC_perCalibSession(1:(iCalib_MVC-1)));
            yThreshold = bottomScaleLimit - graphYSize*(maxMVCuntilNow/maxVoltage);
            Screen('DrawLine', window, orange, leftScaleLimit, yThreshold, rightScaleLimit, yThreshold,5);
        end
        
        Screen(window,'Flip');
        %         [lastFrameTime, timeDispNow]  = Screen('Flip', window, lastFrameTime + (0.5*ifi));
    end % time for the current calibration trial
    
    %% Show a rest text and give some rest
    DrawFormattedText(window, stim.MVC_rest.text,...
        stim.MVC_rest.x, stim.MVC_rest.y, stim.MVC_rest.colour);
    [~,timeNow]  = Screen(window,'Flip');
    onsets.initial_MVC_rest(iCalib_MVC) = timeNow;
    WaitSecs(t_MVC_rest);
    
    %% extract max force for this session (in Voltage)
    MVC_perCalibSession(iCalib_MVC) = max(forceCalib.(['calibTrial_',num2str(iCalib_MVC)])(:,3));
    
    %% stop acquisition of biopac handgrip
    % stop acquiring data in the grip buffer
    stop(dq);
    
    % empty the grip buffer
    flush(dq);
    
    %% display number of trials done for the experimenter
    disp(['Physical calibration trial ',num2str(iCalib_MVC),'/',num2str(n_MVC_repeat),' done']);
end % calibration loop
    
%% store max MVC measure in output
MVC.forceCalib = forceCalib;
MVC.MVC_perCalibSession = MVC_perCalibSession;
MVC.MVC = max(MVC_perCalibSession);

end % function