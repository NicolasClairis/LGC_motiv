function[n_mental_max, calib_summary] = mental_calibNumbers(scr, stim, key,...
    numberVector_calib, mentalE_prm, n_calibTrials, calibTimes, langage)
%[n_mental_max_perTrial, calib_summary] = mental_calibNumbers(scr, stim, key,...
%     numberVector_calib, mentalE_prm, n_calibTrials, calibTimes, langage)
%
% mental_calibNumbers will extract maximum number of subsequent correct
% answers participants can provide in the limited amount of time that is
% available for them to answer.
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
% calibTimes: structure with timings for this phase of the task
%   .instructions: instructions duration
%   .effort_max: time limit for calibration
%   .fbk: time to display feedback during calibration
%
% langage: 'fr' for french; 'engl' for english
%
% n_errorsThreshold: number of errors tolerated before repeating the trial
% (typically for actual calibration should be at zero, but for learning you
% can tolerate a higher number of errors if you don't want the learning to
% be too long)
%
% OUTPUTS
% n_mental_max_perTrial: maximum number of subsequent correct answers per
% trial
%
% calib_summary: structure with relevant data of the calibration procedure
%

%% extract relevant variables
% screen parameters
window = scr.window;
wrapat = scr.wrapat;

% define main parameters
calib_time_limit = true; % time will be limited (as opposed to learning where time was infinite)
instructions_disp = 0; % no instructions anymore, goal is to calibrate as if it was the actual task

% thresholds for calibration
n_errorsThreshold = 5; % number of errors allowed before repeating the trial
n_minCorrectAnswersToReach = 6; % number of correct answers required to consider the trial ok (if too low, then trial is repeated)

% introduce variables of interest
[n_max_calibPerf_perTrial,...
    onset_fbk,...
    onset_fbk_press] = deal( NaN(1, n_calibTrials) ); % maximal perf at each calibration trial

%% display instructions for the calibration
for iInstructionsLoop = 1:2
    DrawFormattedText(window, stim.mentalCalibInstructions.text,...
        stim.mentalCalibInstructions.x, stim.mentalCalibInstructions.y, stim.mentalCalibInstructions.colour, wrapat);
    if iInstructionsLoop == 1
        [~, timeInstru] = Screen(window, 'Flip');
        calib_summary.onset_instructions = timeInstru;
        WaitSecs(calibTimes.instructions);
    elseif iInstructionsLoop == 2
        % allow the participant to start whenever he/she feels ready by
        % pressing a button
        DrawFormattedText(window, stim.pressWhenReady.text,...
            stim.pressWhenReady.x, stim.pressWhenReady.y, stim.pressWhenReady.colour, wrapat);
        [~, timeInstru_bis] = Screen(window, 'Flip');
        calib_summary.onset_instructions_press = timeInstru_bis;
        KbWait; % wait for a button press to go to next phase
        KbReleaseWait; % wait button press to be off to avoid it being recorder as an answer
    end
end

% max (impossible) to reach
n_calibMax = mentalE_prm.n_maxToReachCalib;
% errors handling
errorLimits.useOfErrorThreshold = false;
errorLimits.useOfErrorMapping = false;
% no best performance at first
n_maxReachedUntilNow = [];

%% perform calibration
iCalibTrial = 1;
n_calibTrialsMissed = 0;
idx_calibTrialsMissed = [];
missedTrialData = {};
while iCalibTrial <= n_calibTrials
    
    % indicator whether calibration trial worked or not
    calibTrialSuccess = 0;
    
    %% calibration trial start: finish when max time reached OR when correct number of answers has been provided
    mentalE_perf = mental_effort_perf_Nback(scr, stim, key,...
        numberVector_calib(iCalibTrial,:),...
        mentalE_prm, n_calibMax, instructions_disp, calib_time_limit, calibTimes.effort_max, errorLimits, n_maxReachedUntilNow, n_minCorrectAnswersToReach);
    
    calib_summary.mentalE_perf(iCalibTrial) = mentalE_perf;
    % store current maximum performance
    n_max_calibPerf_perTrial(iCalibTrial) = mentalE_perf.n_correctAnswersForDisplay;
    n_errorsMade_tmp = mentalE_perf.n_errorsMade;
    
    %% provide feedback according to if reached the top or not + prepare for the next phase of the trial
    % force to watch feedback for a short amount of time
    if (n_max_calibPerf_perTrial(iCalibTrial) < n_minCorrectAnswersToReach) ||...
            (n_errorsMade_tmp >= n_errorsThreshold)
        % not enough correct answers OR too much errors
        % => consider the trial as a failure and redo it until a reasonable performance is attained
        n_calibTrialsMissed = n_calibTrialsMissed + 1;
        % store index of the calibration trial missed
        idx_calibTrialsMissed = [idx_calibTrialsMissed, iCalibTrial];
        % store the data of the missed trial because otherwise will be
        % erased
        missedTrialData{n_calibTrialsMissed} = mentalE_perf;
        if (n_max_calibPerf_perTrial(iCalibTrial) < n_minCorrectAnswersToReach)
            DrawFormattedText(window, stim.mentalCalibFailureFbk.text,...
                stim.mentalCalibFailureFbk.x, stim.mentalCalibFailureFbk.y,...
                stim.mentalCalibFailureFbk.colour, wrapat);
        elseif n_errorsMade_tmp >= n_errorsThreshold
            DrawFormattedText(window, stim.mentalCalibFailureTooManyErrorsFbk.text,...
                stim.mentalCalibFailureTooManyErrorsFbk.x, stim.mentalCalibFailureTooManyErrorsFbk.y,...
                stim.mentalCalibFailureTooManyErrorsFbk.colour, wrapat);
        end
    else % calibration worked => display feedback of the best performance until now
        calibTrialSuccess = 1;
        iCalibTrial = iCalibTrial + 1;
        % update best performance score
        n_maxReachedUntilNow = max(n_max_calibPerf_perTrial);
        switch langage
            case 'fr'
                DrawFormattedText(window, ['Bravo! Votre meilleur score jusque-la est de ',...
                    num2str(n_maxReachedUntilNow),' bonnes reponses.'],...
                    stim.mentalCalibFbk.x, stim.mentalCalibFbk.y,...
                    stim.mentalCalibFbk.colour, wrapat);
            case 'engl'
                DrawFormattedText(window, ['Well done! Your best score until now is ',...
                    num2str(n_maxReachedUntilNow),' correct answers.'],...
                    stim.mentalCalibFbk.x, stim.mentalCalibFbk.y,...
                    stim.mentalCalibFbk.colour, wrapat);
        end
    end

    [~, time_fbk] = Screen(window, 'Flip');
    onset_fbk(iCalibTrial) = time_fbk;
    WaitSecs(calibTimes.fbk);

    %% disp trial is done
    if calibTrialSuccess == 1
        disp(['Mental calibration trial ',num2str(iCalibTrial-1),'/',num2str(n_calibTrials),' - done']);
    else
        disp(['Mental calibration trial ',num2str(iCalibTrial+1),' is gonna be repeated because calibration failed.']);
    end
    
    %% allow the participant to restart whenever he/she feels ready by
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

%% get maximum for the participant
n_mental_max = max(n_max_calibPerf_perTrial);

%% store all relevant variables in the output
calib_summary.n_minCorrectAnswersToReach = n_minCorrectAnswersToReach;
calib_summary.calib_time_limit      = calib_time_limit;
calib_summary.t_effort_max          = calibTimes.effort_max;
calib_summary.errorLimits           = errorLimits;
calib_summary.n_calibTrialsToDo     = n_calibTrials;
calib_summary.n_mental_max          = n_mental_max;
calib_summary.n_max_calibPerf_perTrial = n_max_calibPerf_perTrial;
calib_summary.onset_fbk             = onset_fbk;
calib_summary.onset_fbk_press       = onset_fbk_press;
% record data for missed trials
calib_summary.n_calibTrialsMissed   = n_calibTrialsMissed;
calib_summary.idx_calibTrialsMissed = idx_calibTrialsMissed;
calib_summary.missedTrialData       = missedTrialData;

end % function