function[n_mental_max_perTrial, calib_summary] = LGCM_mental_calibNumbers(scr, stim, key,...
    numberVector_calib, mentalE_prm, n_calibTrials, n_calibMax, calibTimes)
%[n_mental_max_perTrial, calib_summary] = LGCM_mental_calibNumbers(scr, stim, key,...
%     numberVector_calib, mentalE_prm, n_calibTrials, n_calibMax, calibTimes)
%
% LGCM_mental_calibNumbers will extract maximum number of subsequent correct
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
% n_calibMax: reference number of subsequent
% correct answers to use
%
% calibTimes: structure with timings for this phase of the task
%   .instructions: instructions duration
%   .effort_max: time limit for calibration
%   .fbk: time to display feedback during calibration
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
yScreenCenter = scr.yCenter;
blackCol = scr.colours.black;
wrapat = scr.wrapat;

% define main parameters
calib_time_limit = true; % time will be limited (as opposed to learning where time was infinite)
instructions_disp = 0; % no instructions anymore, goal is to calibrate as if it was the actual task

% introduce variables of interest
[n_max_calibPerf,...
    onset_fbk,...
    onset_fbk_press] = deal( NaN(1, n_calibTrials) ); % maximal perf at each calibration trial

%% display instructions for the calibration
DrawFormattedText(window,...
    ['Désormais vous devrez répondre dans un temps limité. Essayez de compléter',...
    ' le cercle en répondant aussi vite que possible et correctement aux questions posées.',...
    ' Comme précédemment, une réponse incorrecte réinitialisera le compteur.'],...
    'center', yScreenCenter/3, blackCol, wrapat);
[~, timeInstru] = Screen(window, 'Flip');
calib_summary.onset_instructions = timeInstru;
WaitSecs(calibTimes.instructions);

% allow the participant to start whenever he/she feels ready by
% pressing a button
DrawFormattedText(window,...
    ['Désormais vous devrez répondre dans un temps limité. Essayez de compléter',...
    ' le cercle en répondant aussi vite que possible et correctement aux questions posées.',...
    ' Comme précédemment, une réponse incorrecte réinitialisera le compteur.'],...
    'center', yScreenCenter/3, blackCol, wrapat);
DrawFormattedText(window,...
    'Vous pouvez appuyer quand vous vous sentez prêt(e) à commencer.',...
    'center', yScreenCenter*(5/3), blackCol, wrapat);
[~, timeInstru_bis] = Screen(window, 'Flip');
calib_summary.onset_instructions_press = timeInstru_bis;
KbWait; % wait for a button press to go to next phase
KbReleaseWait; % wait button press to be off to avoid it being recorder as an answer

%% perform calibration
for iCalibTrial = 1:n_calibTrials
    
    %% maybe we could add a short reminder of the colour code here to
    % maximise their performance in the calibration?
    
    
    %% define number of correct sequences to reach to get to top
    switch iCalibTrial
        case 1
            n_max = n_calibMax;
        otherwise % next trials use initial participant max + 3
            max_perf_reached_duringCalib = nanmax(n_max_calibPerf);
            if max_perf_reached_duringCalib < n_calibMax
                n_max = n_calibMax;
            else % increase difficulty to push towards
                n_max = max_perf_reached_duringCalib + 2;
                n_switch = mentalE_prm.switchPerc*n_max;
                % number of switches has to be integer => increase
                % difficulty accordingly
                while n_switch ~= round(n_switch)
                    n_max = n_max + 1;
                    n_switch = mentalE_prm.switchPerc*n_max;
                end
            end
    end
    
    %% calibration trial start: finish when max time reached OR when correct number of answers has been provided
    [mentalE_perf, calibTrial_success] = LGCM_mental_effort_perf(scr, stim, key,...
        numberVector_calib(iCalibTrial,:),...
        mentalE_prm, n_max, instructions_disp, calib_time_limit, calibTimes.effort_max);
    
    calib_summary.mentalE_perf(iCalibTrial) = mentalE_perf;
    % store current maximum performance
    n_max_calibPerf(iCalibTrial) = mentalE_perf.n_questions_correct;
    
    %% provide feedback according to if reached the top or not + prepare for the next phase of the trial
    % force to watch feedback for a short amount of time
    switch calibTrial_success
        case true % reached the top
            DrawFormattedText(window,...
                'Bravo vous avez tout résolu dans le temps imparti!',...
                'center', yScreenCenter/3, blackCol, wrapat);
        case false
            DrawFormattedText(window,...
                'Essayez encore!',...
                'center', yScreenCenter/3, blackCol, wrapat);
    end
    [~, time_fbk] = Screen(window, 'Flip');
    onset_fbk(iCalibTrial) = time_fbk;
    WaitSecs(calibTimes.fbk);
    
    % allow the participant to restart whenever he/she feels ready by
    % pressing a button (no sense for the last trial though)
    if iCalibTrial < n_calibTrials
        DrawFormattedText(window,...
            'Vous pouvez appuyer quand vous vous sentez prêt(e) à recommencer.',...
            'center', yScreenCenter*(5/3), blackCol, wrapat);
        [~, time_fbkPress] = Screen(window, 'Flip');
        onset_fbk_press(iCalibTrial) = time_fbkPress;
        KbWait; % wait for a button press to go to next phase
    end
    KbReleaseWait; % wait button press to be off to avoid it being recorder as an answer
    
end % number of tests to try to get max

%% get maximum for the participant
n_mental_max_perTrial = nanmax(n_max_calibPerf);

%% store all relevant variables in the output
calib_summary.n_mental_max_perTrial = n_mental_max_perTrial;
calib_summary.n_max_calibPerf       = n_max_calibPerf;
calib_summary.onset_fbk             = onset_fbk;
calib_summary.onset_fbk_press       = onset_fbk_press;

end % function