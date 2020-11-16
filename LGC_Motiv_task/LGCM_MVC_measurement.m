function[MVC, onsets] = LGCM_MVC_measurement(scr, session_effort_type, speed, stim)
%[MVC, onsets] = LGCM_MVC_measurement(scr, session_effort_type, speed, stim)
% LGCM_MVC_measurement will display the instructions and measure the MVC
%
% INPUTS
% scr: structure about main screen parameters (size, window, etc.)
%
% session_effort_type: 'physical'/'mental' effort task? => adapt the MVC
% measurement accordingly
%
% speed: structure with stimulus display speed
%
% stim: structure with stimulus properties (size, etc.)
%
% OUTPUTS
% MVC: MVC value
%
% onsets: structure with the onset value of each step of this measure
%
% See also LGCM_main_experiment.m

%% screen relevant variables
window = scr.window;
xCenter = scr.xCenter;
yCenter = scr.yCenter;
screenYsize = yCenter*2;
threshold_1 = scr.threshold_1;
text_size_1 = 100;
text_size_2 = 150;
text_size_3 = 140;
red_colour = [0.8 0 0];
white = [255 255 255]./255;

%% extract timings
t_MVC_calib_instructions1 = 5;
t_MVC_calib_instructions2 = 3;
t_MVC_calib = 3;
t_MVC_rest = 7;

%% number of calibration sessions
n_calib_MVC = 3;

%% MVC 
MVC_perCalibSession = NaN(1,n_calib_MVC);

%% Quick text to introduce MVC calibration
Screen('TextSize', window, text_size_1);
DrawFormattedText(window, 'First : Maximum Power', 'center',screenYsize*0.7, 1);
Screen('TextSize', window, text_size_2)
DrawFormattedText(window, 'Get Ready', 'center', screenYsize*0.3, 1);

[~,time_disp1,~,~,~] = Screen(window,'Flip');
onsets.initial_MVC_msg1 = time_disp1;
speed.vbl  = Screen('Flip', window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
WaitSecs(t_MVC_calib_instructions1);
warning('need to add onset here');

%% Measure MVC as normal
% Reset some variable locally
VCsignal = 0;
centeredCore = CenterRectOnPointd(threshold_1, xCenter, yCenter);

% Set screen text size
Screen('TextSize', scr.window, text_size_3);

%% Measure MVC and keep maximal value
for iCalib_MVC = 1:n_calib_MVC
    
    %% store all force values
    live_force_MVC.(['calib_session_',num2str(iCalib_MVC)]) = [];
    
    %% Draw text and circle
    DrawFormattedText(window, 'Ready', 'center', screenYsize*0.8, white);
    Screen('FillOval', window, red_colour, centeredCore);
    speed.vbl  = Screen('Flip', window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
    pause(t_MVC_calib_instructions2);
    DrawFormattedText(window, 'GO !', 'center', screenYsize*0.8, 1);
    speed.vbl  = Screen('Flip', window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
    
    time_disp0 = GetSecs;
    
    %% During t_MVC_calib second, show signal power and record it
    timenow = GetSecs;
    while timenow < time_disp0 + t_MVC_calib
        
        switch session_effort_type
            case 'physical' % physical effort task
                % Read the force signal from BioPac
                VCsignal = fread(stim.u_out, 1, 'double');
                
                % store all values of force
                live_force_MVC.(['calib_session_',num2str(iCalib_MVC)]) = [live_force_MVC.(['calib_session_',num2str(iCalib_MVC)]),...
                    VCsignal];
                
            case 'mental' % mental effort task
                
        end % physical/mental effort
        
        % Make a base Rectangle of 100 pixel width
        baseRect = [xCenter - 50,...
            yCenter - 85 - VCsignal*4,...
            xCenter + 50,...
            yCenter - 85];
        
        % We set bounds to make sure our square doesn't go completely off of the screen
        for iDimRect = 1:4 % x0/y0/x1/y1
            if baseRect(iDimRect) < 0
                baseRect(iDimRect) = 0;
            end
        end % dim rectangle
        
        % Draw the rectangle and the oval with the text
        Screen('FillOval', window, red_colour, centeredCore);
        Screen('FillRect', window, red_colour, baseRect);
        DrawFormattedText(window, 'GO !', 'center', screenYsize*0.8,1);
        
        % Flip to the screen
        speed.vbl  = Screen('Flip', window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
        warning('need to add onset here');
        
        % check timing to know when to get out of the loop
        timenow = GetSecs;
    end % calibration
    
    % Show a rest text and give some rest
    DrawFormattedText(scr.window, 'Rest', 'center', screenYsize*0.8, [0 0.8 0 ]);
    speed.vbl  = Screen('Flip', scr.window, speed.vbl + (speed.waitframes - 0.5) * speed.ifi);
    warning('need to add onset here');
    WaitSecs(t_MVC_rest);
    
    % extract max force for this session
    MVC_perCalibSession(iCalib_MVC) = nanmax(live_force_MVC.(['calib_session_',num2str(iCalib_MVC)]));
    
end % number of calibrations loop

%% store max MVC measure in output
MVC.MVC_perCalibSession = MVC_perCalibSession;
MVC.MVC = nanmax(MVC_perCalibSession);

end % function