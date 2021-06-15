function[t_min_reached_duringCalib, calib_summary, calib_success] = mental_calibTime(scr, stim, key,...
    numberVector_calib, mentalE_prm, n_calibTrials, n_calibMax, calibTimes, errorLimits, langage)
%[t_min_reached_duringCalib, calib_summary, calib_success] = mental_calibTime(scr, stim, key,...
%     numberVector_calib, mentalE_prm, n_calibTrials, n_calibMax, calibTimes, errorLimits, langage)
%
% mental_calibTime will extract the maximal amount of time requested
% for the hardest level of difficulty (ie the highest number of correct
% answers to provide). Calibration will be repeated until this calibrated
% time is shorter than calibTimes.effort_max.
%
% INPUTS
% scr: structure with screen parameters
%
% stim: stucture with stimuli parameters
%
% key: structure with key code for Psychtoolbox to identify which key
% corresponds to left and right cues
%
% numberVector_calib:
%
% mentalE_prm: structure with main parameters for mental effort task:
%   .mental_n_col: structure with the colour to use for the font
%       .oddEven: colour to use for odd or even question
%       .lowHigh: colour to use for lower/higher than 5 question
%
%   .sideQuestion: structure which tells where each answer is expected to be
%       .oE.pair: -1 means left button corresponds to the pair answer and
%       .oE.impair = 1 means right button corresponds to the impair answer
%       same logic applies to .hL.low and .hL.high fields
%
%   .switchPerc: percentage of switches required (based on total number of
%   subsequent correct answers you want)
%
% n_calibTrials: number of calibration trials
%
% n_calibMax: reference number of subsequent correct answers to use for
% maximally difficult level
%
% calibTimes: structure with timings for this phase of the task
%   .instructions: instructions duration
%   .effort_max: time limit for calibration
%   .fbk: time to display feedback during calibration
%
% errorLimits: structure containing information about way to handle
% errors
%   .useOfErrorThreshold: if true, means the trial is considered a failure,
%   if the number of errors set as a threshold is reached
%   .errorThreshold: consider the trial a failure if more than this number
%   of errors are made
%   .useOfErrorMapping: if true, display the mapping where to answer and
%   type of the trial after a given number of errors has been made
%   .errorMappingLimit: display the mapping after this number of errors has
%   been reached
%
% langage: string indicating the langage with which to display the
% instructions and feedback
% 'fr': french
% 'engl': english
%
% OUTPUTS
% t_min_reached_duringCalib: minimal time necessary to perform all the
% requested questions during the calibration
%
% calib_summary: structure with relevant data of the calibration procedure
%
% calib_success:
% (true) if calibrated time was lower than the
% calibTimes.effortMax in at least one calibration trial;
% (false) if the performance never reached the top in the requested time

%% extract relevant variables
% screen parameters
window = scr.window;
wrapat = scr.wrapat;

% define main parameters
calib_time_limit = true; % time will be limited (as opposed to learning where time was infinite)

% introduce variables of interest
[t_min_calibPerf,...
    onset_fbk,...
    onset_fbk_press] = deal( NaN(1, n_calibTrials) ); % maximal perf at each calibration trial

%% display instructions for the calibration
for iInstructionsLoop = 1:2
    DrawFormattedText(window, stim.mentalCalibInstructions.text,...
        stim.mentalCalibInstructions.x, stim.mentalCalibInstructions.y, stim.mentalCalibInstructions.colour, wrapat);
    % display that participant can press to move on
    if iInstructionsLoop == 2
        DrawFormattedText(window, stim.pressWhenReady.text,...
            stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour, wrapat);
    end
    % display text on screen
    [~, timeInstru] = Screen(window, 'Flip');
    % record information and wait according to loop
    switch iInstructionsLoop
        case 1
            calib_summary.onset_instructions = timeInstru;
            WaitSecs(calibTimes.instructions);
        case 2
            calib_summary.onset_instructions_press = timeInstru;
            % allow the participant to start whenever he/she feels ready by
            % pressing a button
            KbWait; % wait for a button press to go to next phase
            KbReleaseWait; % wait button press to be off to avoid it being recorder as an answer
    end
end % instructions loop

%% perform calibration
for iCalibTrial = 1:n_calibTrials
    
    %% maybe we could add a short reminder of the colour code here to
    % maximise their performance in the calibration?
    
    %% in case you want to adapt the timing based on previous performance
%     switch iCalibTrial
%         case 1
            t_effort_max = calibTimes.effort_max;
%             t_min_reached_duringCalib = t_effort_max;
%         otherwise % next trials use the shortest time until now
%             t_min_reached_duringCalib = nanmin(t_min_calibPerf);
%             % if minimal time reached during calibration
%             if isnan(t_min_reached_duringCalib)
%                 t_min_reached_duringCalib = t_effort_max;
%             end
%     end
    
    %% calibration trial start: finish when max time reached OR when correct number of answers has been provided
    [mentalE_perf, calibTrial_success] = mental_effort_perf_Nback(scr, stim, key,...
        numberVector_calib(iCalibTrial,:),...
        mentalE_prm, n_calibMax,...
        'all','noInstructions', calib_time_limit, t_effort_max, errorLimits); % no instruction (calibration as in the real task)
    
    calib_summary.mentalE_perf(iCalibTrial) = mentalE_perf;
    % store current maximum performance
    t_min_calibPerf(iCalibTrial) = mentalE_perf.totalTime_success;
    
    %% provide feedback according to if reached the top or not + prepare for the next phase of the trial
    % force to watch feedback for a short amount of time
    
    % extract best timing
    t_min_reached_duringCalib = nanmin(t_min_calibPerf);
    t_min_reached_duringCalib_str = sprintf('%0.3f',t_min_reached_duringCalib);
    % display feedback accordingly
    switch calibTrial_success
        case true % reached the top
            switch langage
                case 'fr'
                    DrawFormattedText(window,...
                        ['Bravo vous avez tout resolu dans le temps imparti! ',...
                        'Votre meilleur temps est de ',t_min_reached_duringCalib_str,' s.'],...
                        stim.mentalCalibSuccessFbk.x, stim.mentalCalibSuccessFbk.y,...
                        stim.mentalCalibSuccessFbk.colour, wrapat);
                case 'engl'
                    DrawFormattedText(window,...
                        ['Well done, you solved everything in the allotted time! ',...
                        'Your best timing is ',t_min_reached_duringCalib_str,' s.'],...
                        stim.mentalCalibSuccessFbk.x, stim.mentalCalibSuccessFbk.y,...
                        stim.mentalCalibSuccessFbk.colour, wrapat);
            end
        case false % didn't reach the top in the dedicated time
            if iCalibTrial < n_calibTrials
                DrawFormattedText(window, stim.mentalCalibFailureFbk.text,...
                    stim.mentalCalibFailureFbk.x, stim.mentalCalibFailureFbk.y, stim.mentalCalibFailureFbk.colour, wrapat);
            elseif iCalibTrial == n_calibTrials % last trial
                if ~isnan(t_min_reached_duringCalib)
                    switch langage
                        case 'fr'
                            DrawFormattedText(window,...
                                ['Bravo! Votre meilleur temps est de ',t_min_reached_duringCalib_str,' s.'],...
                                stim.mentalCalibEnd.x, stim.mentalCalibEnd.y, stim.mentalCalibEnd.colour, wrapat);
                        case 'engl'
                            DrawFormattedText(window,...
                                ['Well done! Your best timing is ',t_min_reached_duringCalib_str,' s.'],...
                                stim.mentalCalibEnd.x, stim.mentalCalibEnd.y, stim.mentalCalibEnd.colour, wrapat);
                    end
                else
                    DrawFormattedText(window, stim.mentalCalibFailureFbk.text,...
                        stim.mentalCalibFailureFbk.x, stim.mentalCalibFailureFbk.y, stim.mentalCalibFailureFbk.colour, wrapat);
                end
            end
    end % trial is a success or not?
    if iCalibTrial == n_calibTrials % last trial
        if ~isnan(t_min_reached_duringCalib)
            calib_success = true;
        else
            calib_success = false;
        end
    end % last trial check
    [~, time_fbk] = Screen(window, 'Flip');
    onset_fbk(iCalibTrial) = time_fbk;
    WaitSecs(calibTimes.fbk);
    
    % allow the participant to restart whenever he/she feels ready by
    % pressing a button (no sense for the last trial though)
    if iCalibTrial < n_calibTrials
        DrawFormattedText(window, stim.pressWhenReady.text,...
            stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour, wrapat);
        [~, time_fbkPress] = Screen(window, 'Flip');
        onset_fbk_press(iCalibTrial) = time_fbkPress;
        KbWait; % wait for a button press to go to next phase
    end
    KbReleaseWait; % wait button press to be off to avoid it being recorder as an answer
    
end % number of tests to try to get max

%% store all relevant variables in the output
calib_summary.t_min_reached_duringCalib = t_min_reached_duringCalib;
calib_summary.t_min_calibPerf           = t_min_calibPerf;
calib_summary.onset_fbk                 = onset_fbk;
calib_summary.onset_fbk_press           = onset_fbk_press;
calib_summary.mentalE_perf              = mentalE_perf;

end % function