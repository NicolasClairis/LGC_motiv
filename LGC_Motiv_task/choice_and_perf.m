function[summary] = choice_and_perf(scr, stim, key,...
    effort_type, Ep_or_Em_vars,...
    training_R_P_RP_or_mainTask, nTrials, choiceOptions, confidenceDispChoice,...
    timings,...
    results_folder, file_nm)
% [summary] = choice_and_perf(scr, stim, key,...
%     effort_type, Ep_or_Em_vars,...
%     R_or_P_or_RP_condition, nTrials, choiceOptions, confidenceDispChoice,...
%     timings,...
%     results_folder, file_nm)
% choice_and_perf: script to perform choice and effort performance.
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
% training_R_P_RP_or_mainTask:
% 'R': reward only training
% 'P': punishment only training
% 'RP': reward and punishment training
%
% nTrials: number of trials
%
% choiceOptions: reward and effort level to display for each option (left/right) for
% each trial
%
% confidenceDispChoice: (false/true) indicates if confidence mapping
% should be displayed or not during the choice period (mainly for choices)
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
% black = scr.colours.black;
% yScreenCenter = scr.yCenter;
% xScreenCenter = scr.xCenter;
barTimeWaitRect = stim.barTimeWaitRect;

switch key.n_buttonsChoice
    case 2
        confidenceDispChosen.display = false;
    case 4
        confidenceDispChosen.display = true;
        conf.lowOrHigh = NaN(1,nTrials);
end

%% timings
t_preChoiceCross = timings.preChoiceCross.(training_R_P_RP_or_mainTask);
t_preEffortCross = timings.preEffortCross.(training_R_P_RP_or_mainTask);
% precise if the choice and the performance periods will have a time
% constraint
choiceTimeParameters.timeLimit = true;
timeLimitPerf = true; % time limit to reach level of force required
% if there is a time limit for the choices (and/or the performance), extract the time limit you
% should use
if choiceTimeParameters.timeLimit == true
    choiceTimeParameters.t_choice = timings.choice;
end
if timeLimitPerf == true
    t_max_effort    = timings.max_effort;
else
    t_max_effort = [];
end
t_dispChoice    = timings.dispChoice;
t_fbk           = timings.feedback;

%% specific variables
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

failedTrials = {};
for iTrial = 1:nTrials
    
    %% fixation cross pre-choice period
    Screen('FillRect',window, white, stim.cross.verticalLine); % vertical line
    Screen('FillRect',window, white, stim.cross.horizontalLine); % horizontal line
    [~,onsets.preChoiceCross(iTrial)] = Screen('Flip',window); % display the cross on screen
    WaitSecs(t_preChoiceCross(iTrial));
    dur.preChoiceCross(iTrial) = GetSecs - onsets.preChoiceCross(iTrial);
    
    %% check that no key is being pressed before the choice trial starts
    [was_a_key_pressed_bf_trial(iTrial),...
        onsets.preChoiceCross_keyReleaseMessage(iTrial),...
        dur.preChoiceCross_keyReleaseMessage(iTrial)] = check_keys_are_up(scr, stim, key);
    
    % if a key was pressed before starting the trial => show the fixation
    % cross again with a similar amount of time
    if was_a_key_pressed_bf_trial(iTrial) == 1
        Screen('FillRect',window, white, stim.cross.verticalLine); % vertical line
        Screen('FillRect',window, white, stim.cross.horizontalLine); % horizontal line
        [~,onsets.preChoiceCross_after_buttonRelease(iTrial)] = Screen('Flip',window); % display the cross on screen
        WaitSecs(t_preChoiceCross(iTrial));
        dur.preChoiceCross_after_buttonRelease(iTrial) = GetSecs - onsets.preChoiceCross_after_buttonRelease(iTrial);
    end
    
    %% extract monetary incentive, effort level and reward/punishment condition
    R_left_tmp = choiceOptions.monetary_amount.left(iTrial);
    R_right_tmp = choiceOptions.monetary_amount.right(iTrial);
    defaultSide_tmp = choiceOptions.default_LR(iTrial);
    E_left_tmp = choiceOptions.E.left(iTrial);
    E_right_tmp = choiceOptions.E.right(iTrial);
    R_or_P_tmp = choiceOptions.R_or_P{iTrial};
    
    %% choice period
    if ~strcmp(training_R_P_RP_or_mainTask,'mainTask')
        % for training: keep choice period until a choice is done
        while choice(iTrial) == 0
            [choice(iTrial),...
                onsets.dispChoiceOptions(iTrial),...
                onsets.choice(iTrial),...
                stoptask] = choice_period(scr, stim,...
                R_left_tmp, R_right_tmp, E_left_tmp, E_right_tmp, R_or_P_tmp,...
                choiceTimeParameters, key, confidenceDispChoice);
        end % keep performing the trial until a choice is made
    else % for actual task, if they don't answer in time, consider the trial as a failure
        [choice(iTrial),...
            onsets.dispChoiceOptions(iTrial),...
            onsets.choice(iTrial),...
            stoptask] = choice_period(scr, stim,...
            R_left_tmp, R_right_tmp, E_left_tmp, E_right_tmp, R_or_P_tmp,...
            choiceTimeParameters, key, confidenceDispChoice);
    end
    
    % extract choice made
    switch choice(iTrial)
        case {-2,-1} % choice = left option (sure or unsure)
            R_chosen(iTrial) = R_left_tmp;
            E_chosen(iTrial) = E_left_tmp;
        case {1,2} % choice = right option (sure or unsure)
            R_chosen(iTrial) = R_right_tmp;
            E_chosen(iTrial) = E_right_tmp;
        case 0 % no option was selected: take the default option
            switch defaultSide_tmp
                case -1 % default on the left
                    R_chosen(iTrial) = R_left_tmp;
                    E_chosen(iTrial) = E_left_tmp;
                case 1 % default on the right
                    R_chosen(iTrial) = R_right_tmp;
                    E_chosen(iTrial) = E_right_tmp;
            end
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
        sca; % close PTB
        break; % stop for loop
    end
    
    %% chosen option display period
    [onsets.dispChoice(iTrial)] = choice_task_dispChosen(scr, stim, choice(iTrial), R_chosen(iTrial), E_chosen(iTrial), R_or_P_tmp, confidenceDispChosen);
    WaitSecs(t_dispChoice);
    dur.dispChoice(iTrial) = GetSecs - onsets.dispChoice(iTrial);
    
    %% fixation cross pre-effort period
    Screen('FillRect',window, black, stim.cross.verticalLine); % vertical line
    Screen('FillRect',window, black, stim.cross.horizontalLine); % horizontal line
    [~,onsets.preEffortCross(iTrial)] = Screen('Flip',window); % display the cross on screen
    WaitSecs(t_preEffortCross(iTrial));
    dur.preEffortCross(iTrial) = GetSecs - onsets.preEffortCross(iTrial);
    
    %% check that no key is being pressed before the choice trial starts
    [was_a_key_pressed_bf_effort(iTrial),...
        onsets.preEffortCross_keyReleaseMessage(iTrial),...
        dur.preEffortCross_keyReleaseMessage(iTrial)] = check_keys_are_up(scr, stim, key);
    
    % if a key was pressed before starting the trial => show the fixation
    % cross again with a similar amount of time
    if was_a_key_pressed_bf_effort(iTrial) == 1
        Screen('FillRect',window, black, stim.cross.verticalLine); % vertical line
        Screen('FillRect',window, black, stim.cross.horizontalLine); % horizontal line
        [~,onsets.preEffortCross_after_buttonRelease(iTrial)] = Screen('Flip',window); % display the cross on screen
        WaitSecs(t_preEffortCross(iTrial));
        dur.preEffortCross_after_buttonRelease(iTrial) = GetSecs - onsets.preEffortCross_after_buttonRelease(iTrial);
    end
    
    %% Effort period
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
                timeLimitPerf, timings, R_or_P_tmp, R_chosen(iTrial));
            % record duration for effort performance
            dur.effortPeriod(iTrial) = GetSecs - onsets.effortPeriod{iTrial}.effort_phase;
            
        case 'mental'
            mentalE_prm.startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen(iTrial))]); % adapt start angle to current level of difficulty
            n_max_to_reach_perTrial(iTrial) = n_to_reach.(['E_level_',num2str(E_chosen(iTrial))]);
            [perfSummary{iTrial},...
                trial_was_successfull(iTrial),...
                onsets.effortPeriod{iTrial}] = mental_effort_perf_Nback(scr, stim, key,...
                mental_nbers_per_trial(iTrial,:),...
                mentalE_prm,...
                n_max_to_reach_perTrial(iTrial),...
                'noInstructions', timeLimitPerf, t_max_effort, errorLimits, [], [],...
                R_or_P_tmp, R_chosen(iTrial));
            % record duration for effort performance
            dur.effortPeriod(iTrial) = GetSecs - onsets.effortPeriod{iTrial}.nb_1;
    end % effort type loop
    effortTime(iTrial) = toc;
    
    %% Feedback period
    % compute gains/losses for the current trial
%     if choice(iTrial) == defaultSide_tmp || choice(iTrial) == 0 % if default chosen, result is binary (success or lose)
%         switch trial_was_successfull(iTrial)
%             case 0 % loss
%                 switch R_or_P_tmp
%                     case 'R'
%                         gain(iTrial) = 0;
%                     case 'P'
%                         gain(iTrial) = -R_chosen(iTrial)*2;
%                 end
%             case 1 % gain
%                 switch R_or_P_tmp
%                     case 'R'
%                         gain(iTrial) = R_chosen(iTrial);
%                     case 'P'
%                         gain(iTrial) = -R_chosen(iTrial);
%                 end
%         end % trial successfull or not?
        
%     else % if non-default chosen, result is proportional to performance
        ratioPerf(iTrial) = perfSummary{iTrial}.performance/100;
        switch R_or_P_tmp
            case 'R' % gain between 0 and R_chosen depending on performance
                gain(iTrial) = round(R_chosen(iTrial)*ratioPerf(iTrial),2);
            case 'P' % loss between -R_chosen and 2*(-R_chosen) depending on performance
                gain(iTrial) = round(-R_chosen(iTrial) - (1 - ratioPerf(iTrial))*R_chosen(iTrial),2);
        end
%     end % filter if performance was for the default or the non-default option
    
    % display feedback
    switch R_or_P_tmp
        case 'R'
            DrawFormattedText(window, stim.feedback.reward.text,...
                stim.feedback.reward.x, stim.feedback.reward.y,...
                stim.feedback.colour);
        case 'P'
            DrawFormattedText(window, stim.feedback.punishment.text,...
                stim.feedback.punishment.x, stim.feedback.punishment.y,...
                stim.feedback.colour);
    end
    % display amount of money obtained/loss
    drawRewardAmount(scr, stim, abs(gain(iTrial)), R_or_P_tmp, 'middle_center_start');
    
    [~,onsets.fbk(iTrial)] = Screen(window,'Flip');
    switch R_or_P_tmp
        case 'R'
            onsets.fbk_win(iTrial) = onsets.fbk(iTrial);
        case 'P'
            onsets.fbk_loss(iTrial) = onsets.fbk(iTrial);
    end
    
    % update monetary amount earned until now
    totalGain(iTrial) = sum(gain(1:iTrial));
    WaitSecs(t_fbk);
    
    % extract feedback duration
    dur.fbk(iTrial) = GetSecs - onsets.fbk(iTrial);
    switch R_or_P_tmp
        case 'R'
            dur.fbk_win(iTrial) = GetSecs - onsets.fbk_win(iTrial);
        case 'P'
            dur.fbk_loss(iTrial) = GetSecs - onsets.fbk_loss(iTrial);
    end
    
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
            % record timing
            dur.timeBarWait(iTrial) = GetSecs - onsets.timeBarWait(iTrial);
        end % if all time taken, no need for time penalty
    end % if a time limit is added
    
    %% record information of trial when choice was failed
    if choice(iTrial) == 0
        % trial information
        failedTrials{iTrial}.R_or_P                      = R_or_P_tmp;
        failedTrials{iTrial}.conf                      = conf.lowOrHigh(iTrial);
        failedTrials{iTrial}.R_defaultToDo             = R_chosen(iTrial);
        failedTrials{iTrial}.E_defaultToDo             = E_chosen(iTrial);
        failedTrials{iTrial}.defaultSide_LR            = defaultSide_tmp;
        % timings
        failedTrials{iTrial}.onsets.preChoiceCross       = onsets.preChoiceCross(iTrial);
        failedTrials{iTrial}.dur.preChoiceCross          = dur.preChoiceCross(iTrial);
        if was_a_key_pressed_bf_trial(iTrial) == 1
            failedTrials{iTrial}.onsets.preChoiceCross_after_buttonRelease = onsets.preChoiceCross_after_buttonRelease(iTrial);
            failedTrials{iTrial}.dur.preChoiceCross_after_buttonRelease = dur.preChoiceCross_after_buttonRelease(iTrial);
        end
        failedTrials{iTrial}.onsets.dispChoiceOptions    = onsets.dispChoiceOptions(iTrial);
        failedTrials{iTrial}.dur.dispChoiceOptions       = dur.dispChoiceOptions(iTrial);
        failedTrials{iTrial}.onsets.choice               = onsets.choice(iTrial);
        failedTrials{iTrial}.onsets.dispChoice           = onsets.dispChoice(iTrial);
        failedTrials{iTrial}.dur.dispChoice              = dur.dispChoice(iTrial);
        failedTrials{iTrial}.onsets.preEffortCross       = onsets.preEffortCross(iTrial);
        failedTrials{iTrial}.dur.preEffortCross          = dur.preEffortCross(iTrial);
        if was_a_key_pressed_bf_effort(iTrial) == 1
            failedTrials{iTrial}.onsets.preEffortCross_after_buttonRelease = onsets.preEffortCross_after_buttonRelease(iTrial);
            failedTrials{iTrial}.dur.preEffortCross_after_buttonRelease = dur.preEffortCross_after_buttonRelease;
        end
        failedTrials{iTrial}.onsets.effortPeriod         = onsets.effortPeriod(iTrial);
        failedTrials{iTrial}.dur.effortPeriod            = dur.effortPeriod(iTrial);
        failedTrials{iTrial}.perfSummary                 = perfSummary{iTrial};
        failedTrials{iTrial}.onsets.fbk                  = onsets.fbk(iTrial);
        failedTrials{iTrial}.dur.fbk                     = dur.fbk(iTrial);
        if timeRemainingEndTrial_ONOFF == 1
            failedTrials{iTrial}.onsets.timeBarWait = onsets.timeBarWait(iTrial);
            failedTrials{iTrial}.dur.timeBarWait = dur.timeBarWait(iTrial);
        end
    end
    
    %% display number of trials done for the experimenter
    disp(['Trial ',num2str(iTrial),'/',num2str(nTrials),' done']);
end % trial loop

%% extract relevant training data
summary.onsets = onsets;
summary.durations = dur;
summary.choice = choice;
if confidenceDispChosen.display == true
    summary.confidence = conf;
end
summary.confidence_bis = confidenceDispChosen;
summary.choiceOptions = choiceOptions;
summary.R_chosen = R_chosen;
summary.E_chosen = E_chosen;
summary.perfSummary = perfSummary;
summary.gain = gain;
summary.percentagePerf = ratioPerf;
summary.totalGain = totalGain;
summary.trial_was_successfull = trial_was_successfull;
summary.effortTime = effortTime;
summary.failedTrials = failedTrials;
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