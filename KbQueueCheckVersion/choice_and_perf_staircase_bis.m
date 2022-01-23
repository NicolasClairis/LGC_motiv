function[summary] = choice_and_perf_staircase_bis(scr, stim, key,...
    effort_type, Ep_or_Em_vars,...
    R_or_P, E_right, E_left, nTrials, timings,...
    results_folder, file_nm)
% [summary] = choice_and_perf_staircase_bis(scr, stim, key,...
%     effort_type, Ep_or_Em_vars,...
%     R_or_P, E_right, E_left, nTrials, timings,...
%     results_folder, file_nm)
% choice_and_perf_staircase_bis: script to perform choice and effort performance.
% choice_and_perf_staircase_bis is similar to choice_and_perf_staircase but
% replaces the call to all functions using KbCheck by equivalent functions
% using KbQueueCheck
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with information about the stimuli
%
% key: structure with relevant keys to use for the experiment
%
% effort_type: information about which effort to perform after the choice
% 'mental': mental effort task
% 'physical': physical effort task
%
% Ep_or_Em_vars: structure with variables specific to physical or to mental
% effort task
%   MVC: maximal voluntary contraction force in the physical effort task
% effort calibration
%   i_sub: subject number (for mental effort)
%   n_to_reach: structure telling the number of correct answers to reach to
%   for each effort level
%
% R_or_P: reward ('R') or punishment ('P') trial
%
% E_right, E_left: effort level for right and left option
%
% nTrials: number of trials
%
% timings: structure with information about the relevant timings
%
% results_folder: path where data needs to be saved
%
% file_nm: file name for results
%
% OUTPUTS
% summary: structure with most relevant variables extracted during
% the performance

%% load main paramaters
window = scr.window;
white = scr.colours.white;
black = scr.colours.black;
barTimeWaitRect = stim.barTimeWaitRect;

% confidence feedback visual display
confidenceDispChosen.display = true;
conf.lowOrHigh = NaN(1,nTrials);
% no display of the confidence mapping
confidenceChoiceDisp = false;

%% timings
timings.cross.mainTask = 0.5;
% precise if the choice and the performance periods will have a time
% constraint
choiceTimeParameters.timeLimit = false;
timeLimitPerf = true; % time limit to reach level of force required
% if there is a time limit for the choices (and/or the performance), extract the time limit you
% should use
if timeLimitPerf == true
    t_max_effort    = timings.max_effort;
else
    t_max_effort = [];
end
t_dispChoice    = timings.dispChoice;
t_fbk           = timings.feedback;

% specific variables
switch effort_type
    case 'mental'
        n_to_reach = Ep_or_Em_vars.n_to_reach;
        errorLimits = Ep_or_Em_vars.errorLimits;
    case 'physical'
        MVC = Ep_or_Em_vars.MVC;
        dq = Ep_or_Em_vars.dq;
        Ep_time_levels = Ep_or_Em_vars.Ep_time_levels;
        F_threshold = Ep_or_Em_vars.F_threshold;
        F_tolerance = Ep_or_Em_vars.F_tolerance;
end
timeRemainingEndTrial_ONOFF = Ep_or_Em_vars.timeRemainingEndTrial_ONOFF;

%% initialize onsets
[onsets.preChoiceCross,...
    onsets.dispChoiceOptions,...
    onsets.choice,...
    onsets.preChoiceCross_keyReleaseMessage,...
    onsets.preChoiceCross_after_buttonRelease,...
    onsets.preEffortCross,...
    onsets.preEffortCross_keyReleaseMessage,...
    onsets.preEffortCross_after_buttonRelease,...
    onsets.fbk, onsets.fbk_win, onsets.fbk_loss,...
    onsets.fbk_fail,...
    onsets.timeBarWait,...
    dur.preChoiceCross,...
    dur.dispChoiceOptions,...
    dur.preChoiceCross_keyReleaseMessage,...
    dur.preChoiceCross_after_buttonRelease,...
    dur.preEffortCross,...
    dur.preEffortCross_keyReleaseMessage,...
    dur.preEffortCross_after_buttonRelease,...
    dur.effortPeriod,...
    dur.fbk, dur.fbk_win, dur.fbk_loss,...
    dur.fbk_fail,...
    dur.timeBarWait] = deal(NaN(1,nTrials));
% variables during effort period should record the data for each trial
[onsets.effortPeriod,...
    perfSummary] = deal(cell(1, nTrials));

%% launch main task
choice = zeros(1,nTrials);

[R_chosen, E_chosen,...
    effortTime,...
    trial_was_successfull,...
    gain, ratioPerf,...
    totalGain] = deal(NaN(1, nTrials));
[was_a_key_pressed_bf_trial,...
    was_a_key_pressed_bf_effort] = deal(NaN(1,nTrials));
switch effort_type
    case 'mental'
        % initialize main parameters of the task
        mentalE_prm = mental_effort_parameters();
        % randomize the order of the numbers appearing on screen
        mental_nbers_per_trial = mental_numbers(nTrials);
        % number of good answers to reach at each trial
        n_max_to_reach_perTrial = NaN(1,nTrials);
end

%% define monetary incentive, effort level and reward/punishment condition

% initialize left-right values, depending on the condition
R_right = NaN(1,nTrials+1);
if strcmp('R',R_or_P)
    R_left = 0.50;
    R_right(1) = 1;
    R_right_tmp = R_right(1);
elseif strcmp('P',R_or_P)
    R_left = 0.5;
    R_right(1) = 0.05;
end
% remember initial baseline value for the right value as it will change due to staircase
R_right_baseline = R_right(1);
failed_trials = {};

for iTrial = 1:nTrials
    
    %% fixation cross period
    Screen('FillRect',window, white, stim.cross.verticalLine); % vertical line
    Screen('FillRect',window, white, stim.cross.horizontalLine); % horizontal line
    [~,onsets.preChoiceCross(iTrial)] = Screen('Flip',window); % display the cross on screen
    WaitSecs(2);
    dur.preChoiceCross(iTrial) = GetSecs - onsets.preChoiceCross(iTrial);
    
    %% check that no key is being pressed before the choice trial starts
    [was_a_key_pressed_bf_trial(iTrial),...
        onsets.keyReleaseMessage(iTrial),...
        dur.preChoiceCross_keyReleaseMessage(iTrial)] = check_keys_are_up_bis(scr, stim, key);
    
    % if a key was pressed before starting the trial => show the fixation
    % cross again with a similar amount of time
    if was_a_key_pressed_bf_trial(iTrial) == 1
        Screen('FillRect',window,white, stim.cross.verticalLine); % vertical line
        Screen('FillRect',window,white, stim.cross.horizontalLine); % horizontal line
        [~,onsets.cross_after_buttonRelease(iTrial)] = Screen('Flip',window); % display the cross on screen
        WaitSecs(1);
        dur.preChoiceCross_after_buttonRelease(iTrial) = GetSecs - onsets.preChoiceCross_after_buttonRelease(iTrial);
    end
    
    %% choice period
    % keep choice period until a choice is done
    while choice(iTrial) == 0
        [choice(iTrial),...
            onsets.dispChoiceOptions(iTrial),...
            onsets.choice(iTrial),...
            stoptask] = choice_period_bis(scr, stim,...
            R_left, R_right_tmp, E_left, E_right, R_or_P,...
            choiceTimeParameters, key, confidenceChoiceDisp);
    end % keep performing the trial until a choice is made
    
    % extract choice made
    switch choice(iTrial)
        case {-2,-1} % choice = left option
            R_chosen(iTrial) = R_left;
            E_chosen(iTrial) = E_left;
        case {1,2} % choice = right option
            R_chosen(iTrial) = R_right_tmp;
            E_chosen(iTrial) = E_right;
        case 0 % no option was selected
            R_chosen(iTrial) = 0;
            E_chosen(iTrial) = 0;
    end
    
    switch R_or_P
        
        % If it is a reward trial
        case 'R'
            % if right values becomes higher than baseline (cause they were lazy) put back baseline
            switch choice(iTrial)
                case {-2,-1}
                    R_right_tmp = R_right_tmp + (R_right_tmp - R_left)/2;
                case {1,2}
                    R_right_tmp = R_right_tmp - (R_right_tmp - R_left)/2;
                otherwise
                    error('Il y a un bug dans la partie choix');
            end
            % In case computed value is higher than baseline, put it back to baseline
            if R_right_baseline <= R_right_tmp
                R_right_tmp = R_right_baseline;
            end
            
            % if it is a punishment trial
        case 'P'
            % case it is a punishment trial
            switch choice(iTrial)
                case {-2,-1}
                    R_right_tmp = R_right_tmp - (R_left - R_right_tmp)/2;
                case {1,2}
                    R_right_tmp = R_right_tmp + (R_left - R_right_tmp)/2;
                otherwise
                    error('Il y a un bug dans la partie choix');
            end
            
            % In case computed value is lower than baseline, put it back to baseline
            if R_right_baseline >= R_right_tmp
                R_right_tmp = R_right_baseline;
            end
    end
    
    % record values of right option across trials
    R_right(iTrial + 1) = R_right_tmp;
    
    % save the Indifference point, considered as the last choice to the right after computation
    if iTrial == nTrials
        IP = R_right_tmp;
        delta_IP = abs(IP - R_left);
    end
    
    % in the case where confidence is measured, also extract confidence
    % level
    if confidenceDispChosen.display == true
        switch choice(iTrial)
            case {-2,2} % high confidence
                conf.lowOrHigh(iTrial) = 2;
            case {-1,1} % low confidence
                conf.lowOrHigh(iTrial) = 1;
            otherwise % no choice made = as if low confidence
                conf.lowOrHigh(iTrial) = 0;
        end
        confidenceDispChosen.lowOrHigh = conf.lowOrHigh(iTrial);
    end
    
    % record timing
    dur.dispChoiceOptions(iTrial) = GetSecs - onsets.dispChoiceOptions(iTrial);
    
    %% check if escape was pressed => stop everything if so
    if stoptask == 1
        % save all the data in case you still want to analyze it
        save([results_folder, file_nm,'_earlyEnd_tmp.mat']);
        break;
    end
    
    %% chosen option display period
    [time_dispChoice] = choice_task_dispChosen(scr, stim, choice(iTrial), R_chosen(iTrial), E_chosen(iTrial), R_or_P, confidenceDispChosen);
    onsets.dispChoice(iTrial) = time_dispChoice;
    WaitSecs(t_dispChoice);
    dur.dispChoice(iTrial) = GetSecs - onsets.dispChoice(iTrial);
    
    %% fixation cross pre-effort period
    Screen('FillRect',window, black, stim.cross.verticalLine); % vertical line
    Screen('FillRect',window, black, stim.cross.horizontalLine); % horizontal line
    [~,onsets.preEffortCross(iTrial)] = Screen('Flip',window); % display the cross on screen
    WaitSecs(1);
    dur.preEffortCross(iTrial) = GetSecs - onsets.preEffortCross(iTrial);
    
    %% check that no key is being pressed before the effort starts
    [was_a_key_pressed_bf_effort(iTrial),...
        onsets.preEffortCross_keyReleaseMessage(iTrial),...
        dur.preEffortCross_keyReleaseMessage(iTrial)] = check_keys_are_up_bis(scr, stim, key);
    
    % if a key was pressed before starting the trial => show the fixation
    % cross again with a similar amount of time
    if was_a_key_pressed_bf_effort(iTrial) == 1
        Screen('FillRect',window, black, stim.cross.verticalLine); % vertical line
        Screen('FillRect',window, black, stim.cross.horizontalLine); % horizontal line
        [~,onsets.preEffortCross_after_buttonRelease(iTrial)] = Screen('Flip',window); % display the cross on screen
        WaitSecs(0.5);
        dur.preEffortCross_after_buttonRelease(iTrial) = GetSecs - onsets.preEffortCross_after_buttonRelease(iTrial);
    end
    
    %% perform the effort
    tic;
    switch effort_type
        case 'physical'
            [perfSummary{iTrial},...
                trial_was_successfull(iTrial),...
                onsets.effortPeriod{iTrial}] = physical_effort_perf(scr, stim, dq,...
                MVC,...
                E_chosen(iTrial),...
                Ep_time_levels,...
                F_threshold, F_tolerance,...
                timeLimitPerf, timings, 'R', R_chosen(iTrial));
            % record duration for effort performance
            dur.effortPeriod(iTrial) = GetSecs - onsets.effortPeriod{iTrial}.effort_phase;
            
        case 'mental'
            mentalE_prm.startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen(iTrial))]); % adapt start angle to current level of difficulty
            n_max_to_reach_perTrial(iTrial) = n_to_reach.(['E_level_',num2str(E_chosen(iTrial))]);
            [perfSummary{iTrial},...
                trial_was_successfull(iTrial),...
                onsets.effortPeriod{iTrial}] = mental_effort_perf_Nback(scr, stim, key,...
                mental_nbers_per_trial(iTrial,:),...
                mentalE_prm, n_max_to_reach_perTrial(iTrial),...
                'noInstructions', timeLimitPerf, t_max_effort,errorLimits, [], [],...
                'R', R_chosen(iTrial));
            % record duration for effort performance
            dur.effortPeriod(iTrial) = GetSecs - onsets.effortPeriod{iTrial}.nb_1;
    end % effort type
    effortTime(iTrial) = toc;
    
    %% Feedback period
    ratioPerf(iTrial) = perfSummary{iTrial}.performance/100;
    gain(iTrial) = round(R_chosen(iTrial)*ratioPerf(iTrial),2);
    % reward fbk
    DrawFormattedText(window, stim.feedback.reward.text,...
        stim.feedback.reward.x, stim.feedback.reward.y,...
        stim.feedback.colour);
    % display amount of money obtained/loss
    drawRewardAmount(scr, stim, abs(gain(iTrial)), 'R', 'middle_center_start');
    [~,onsets.fbk(iTrial)] = Screen(window,'Flip');
    onsets.fbk_win(iTrial) = onsets.fbk(iTrial);
    
    % update monetary amount earned until now
    totalGain(iTrial) = sum(gain(1:iTrial));
    WaitSecs(t_fbk);
    
    %% Time waiting period
    % this period allows to de-confound effort and delay, ie even if lower
    % effort has been selected and performed quicker, a short waiting time
    % will force to wait
    if timeRemainingEndTrial_ONOFF == 1
        if effortTime(iTrial) < t_max_effort
            tic;
            onsets.timeBarWait(iTrial) = GetSecs; % record start of time bar
            
            % show a dynamic waiting bar until the timing ends
            while toc <= (t_max_effort - effortTime(iTrial))
                timeSinceStart_tmp = toc;
                % update bar with time remaining
                percTimeAchieved = (timeSinceStart_tmp + effortTime(iTrial))./t_max_effort;
                barTimeWaitRect_bis = barTimeWaitRect;
                % start on the right corner of the bar + percentage already
                % achieved and move to the left
                if percTimeAchieved > 0 && percTimeAchieved < 1
                    barTimeWaitRect_bis(3) = barTimeWaitRect(3) - percTimeAchieved*(barTimeWaitRect(3) - barTimeWaitRect(1));
                elseif percTimeAchieved > 1
                    warning('you should get out of the loop when the time spent is too long but it seems there was a bug, display of timebar was locked to zero to compensate');
                    barTimeWaitRect_bis(3) = barTimeWaitRect(1) + 1;
                end
                %
                DrawFormattedText(window, stim.remainingTime.text, stim.remainingTime.x, stim.remainingTime.y, stim.remainingTime.colour);
                % draw one global fixed rectangle showing the total duration
                Screen('FrameRect',window, stim.barTimeWait.colour, barTimeWaitRect);
                
                % draw one second rectangle updating dynamically showing the
                % time remaining
                Screen('FillRect',window, stim.barTimeWait.colour, barTimeWaitRect_bis);
                
                Screen(window,'Flip');
            end % display until time catches up with maximum effort time
        end % if all time taken, no need for time penalty
    end % if a time limit is added
    
    %% display number of trials done for the experimenter
    disp(['IP trial ',num2str(iTrial),'/',num2str(nTrials),' done']);
end % trial loop

%% extract relevant training data
summary.onsets = onsets;
if confidenceDispChosen.display == true
    summary.confidence = conf;
end
summary.confidence_bis = confidenceDispChosen;
summary.IP = IP;
summary.delta_IP = delta_IP;
summary.onsets = onsets;
summary.E_left = E_left;
summary.E_right = E_right;
summary.R_left = R_left;
summary.R_right = R_right;
summary.R_chosen = R_chosen;
summary.E_chosen = E_chosen;
summary.optionChosen = choice;
summary.perfSummary = perfSummary;
summary.gain = gain;
summary.percentagePerf = ratioPerf;
summary.totalGain = totalGain;
summary.trial_was_successfull = trial_was_successfull;
summary.effortTime = effortTime;
summary.failed_trials = failed_trials;
switch effort_type
    case 'mental'
        summary.mentalE_prm = mentalE_prm;
        summary.n_max_to_reach_perTrial = n_max_to_reach_perTrial;
        summary.n_to_reach = n_to_reach;
    case 'physical'
        summary.MVC = MVC;
        summary.dq = dq;
        summary.Ep_time_levels = Ep_time_levels;
        summary.F_threshold = F_threshold;
        summary.F_tolerance = F_tolerance;
end

%% save all the data in case of crash later on
save([results_folder, file_nm,'_behavioral_tmp.mat'],'summary');
end % function