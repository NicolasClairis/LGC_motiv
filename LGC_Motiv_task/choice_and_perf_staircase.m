function[summary] = choice_and_perf_staircase(scr, stim, key,...
    effort_type, Ep_or_Em_vars, R_money,...
    training_R_P_RP_or_mainTask, R_or_P, E_right, E_left, nTrials, timings,...
    results_folder, file_nm)
% [summary] = choice_and_perf_staircase(scr, stim, key,...
%     effort_type, Ep_or_Em_vars, R_money,...
%     training_R_P_RP_or_mainTask, R_or_P, E_right, E_left, nTrials, timings,...
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
% R_money: structure with equivalence between reward levels and reward
% money to compute gain within the session and includes money loss amount
% for errors
%
% training_R_P_RP_or_mainTask:
% 'R': reward only training
% 'P': punishment only training
% 'RP': reward and punishment training
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
barTimeWaitRect = stim.barTimeWaitRect;

% no confidence display for the staircase pilots
confidence.display = false;

%% timings
t_cross         = timings.cross.(training_R_P_RP_or_mainTask);
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
t_fail_and_repeat_fbk = timings.fail_and_repeat_fbk;

% specific variables
switch effort_type
    case 'mental'
        i_sub = Ep_or_Em_vars.i_sub;
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
[onsets.cross,...
    onsets.dispChoiceOptions,...
    onsets.choice,...
    onsets.keyReleaseMessage,...
    onsets.cross_after_buttonRelease,...
    onsets.fbk, onsets.fbk_win, onsets.fbk_loss,...
    onsets.timeBarWait] = deal(NaN(1,nTrials));
% variables during effort period should record the data for each trial
[onsets.effortPeriod,...
    perfSummary] = deal(cell(1, nTrials));

%% launch main task
choice = zeros(1,nTrials);

[R_chosen, E_chosen,...
    effortTime,...
    trial_was_successfull,...
    gain,...
    totalGain] = deal(NaN(1, nTrials));
was_a_key_pressed_bf_trial = NaN(1,nTrials);
switch effort_type
    case 'mental'
        % initialize main parameters of the task
        mentalE_prm = mental_effort_parameters(i_sub);
        % randomize the order of the numbers appearing on screen
        mental_nbers_per_trial = mental_numbers(nTrials);
        
        % randomize the type of the first trial (odd/even or higher/lower
        % than 5)
        %         mental_taskType_trialStart = mental_task_start(nTrials);
        % number of good answers to reach at each trial
        n_max_to_reach_perTrial = NaN(1,nTrials);
end

%% define monetary incentive, effort level and reward/punishment condition

% initialize left-right values, depending on the condition
if strcmp('R',R_or_P)
    R_left = 1.5;
    R_right_tmp = 2.5;
elseif strcmp('P',R_or_P)
    R_left = 1.5;
    R_right_tmp = 0.5;
end
% remember initial baseline value for the right value as it will change due to staircase
R_right_baseline = R_right_tmp;
failed_trials = {};

for iTrial = 1:nTrials
    
    %% fixation cross period
    Screen('FillRect',window, white, stim.cross.verticalLine); % vertical line
    Screen('FillRect',window, white, stim.cross.horizontalLine); % horizontal line
    [~,onsets.cross(iTrial)] = Screen('Flip',window); % display the cross on screen
    WaitSecs(t_cross(iTrial));
    
    %% check that no key is being pressed before the choice trial starts
    [was_a_key_pressed_bf_trial(iTrial),...
        onsets.keyReleaseMessage(iTrial)] = check_keys_are_up(scr, stim, key);
    
    % if a key was pressed before starting the trial => show the fixation
    % cross again with a similar amount of time
    if was_a_key_pressed_bf_trial(iTrial) == 1
        Screen('FillRect',window,white, stim.cross.verticalLine); % vertical line
        Screen('FillRect',window,white, stim.cross.horizontalLine); % horizontal line
        [~,onsets.cross_after_buttonRelease(iTrial)] = Screen('Flip',window); % display the cross on screen
        WaitSecs(t_cross(iTrial));
    end
    
    % initialize variables in case of failure on this trial
    i_trial_failed = 0;
    trial_success = 0;
    redo_limit = 5;
    %% choice period
    if ~strcmp(training_R_P_RP_or_mainTask,'mainTask')
        % for training: keep choice period until a choice is done
        while choice(iTrial) == 0
            [choice(iTrial),...
                onsets.dispChoiceOptions(iTrial),...
                onsets.choice(iTrial),...
                stoptask] = choice_period(scr, stim,...
                R_left, R_right_tmp, E_left, E_right, R_or_P,...
                choiceTimeParameters, key);
        end % keep performing the trial until a choice is made
    else % for actual task, if they don't answer in time, consider the trial as a failure
        [choice(iTrial),...
            onsets.dispChoiceOptions(iTrial),...
            onsets.choice(iTrial),...
            stoptask] = choice_period(scr, stim,...
            R_left, R_right_tmp, E_left, E_right, R_or_P,...
            choiceTimeParameters, key);
    end
    
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
                case 0
                    % if no choice was made, keep same value
                    R_right_tmp = R_right_tmp;
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
                case 0
                    % if no choice was made, keep same value for next trial
                    R_right_tmp = R_right_tmp;
                otherwise
                    error('Il y a un bug dans la partie choix');
            end
            
            % In case computed value is lower than baseline, put it back to baseline
            if R_right_baseline >= R_right_tmp
                R_right_tmp = R_right_baseline;
            end
    end
    
    % save the Indifference point, considered as the last choice to the right after computation
    if iTrial == 5
        IP = R_right_tmp;
    end
    
    %% check if escape was pressed => stop everything if so
    if stoptask == 1
        % save all the data in case you still want to analyze it
        save([results_folder, file_nm,'_earlyEnd_tmp.mat']);
        break;
    end
    
    %% chosen option display period
    [time_dispChoice] = choice_task_dispChosen(scr, stim, R_chosen(iTrial), E_chosen(iTrial), R_or_P, confidence);
    onsets.dispChoice(iTrial) = time_dispChoice;
    WaitSecs(t_dispChoice);
    
    %% Effort period
    % perform the effort if a choice was made, otherwise punish the
    % subject for not answering to the choice
    if choice(iTrial) == 0 % no choice was made => failure
        trial_was_successfull(iTrial) = 0;
        
    elseif ismember(choice(iTrial), [-2, -1, 1, 2]) % choice done => perform the corresponding effort
        
        %% mental effort: check no key is being pressed before the start of the effort period
        % for physical effort: useless since only the grip is required
        if strcmp(effort_type,'mental')
            [was_a_key_pressed_bf_trial(iTrial),...
                onsets.keyReleaseMessage(iTrial)] = check_keys_are_up(scr, stim, key);
        end
        
        %% perform the effort
        tic;
        switch effort_type
            case 'physical'
                while trial_success == 0 && i_trial_failed <= redo_limit
                    [perfSummary{iTrial},...
                        trial_was_successfull(iTrial),...
                        onsets.effortPeriod{iTrial}] = physical_effort_perf(scr, stim, dq,...
                        MVC,...
                        E_chosen(iTrial),...
                        Ep_time_levels,...
                        F_threshold, F_tolerance,...
                        timeLimitPerf, timings);
                    trial_success = trial_was_successfull(iTrial);
                    if trial_success == 0
                        i_trial_failed = i_trial_failed +1;
                        failed_trials{i_trial_failed}.perfSummary = perfSummary{iTrial};
                        failed_trials{i_trial_failed}.trial_was_successfull = trial_was_successfull(iTrial);
                        failed_trials{i_trial_failed}.onset.effortPeriod =  onsets.effortPeriod{iTrial};
                        failed_trials{i_trial_failed}.i_trial_idx = iTrial;
                        % for the physical effort, when participant was too
                        % slow
                        
                        DrawFormattedText(window, stim.feedback.error_tooSlow.text,...
                            stim.feedback.error_tooSlow.x, stim.feedback.error_tooSlow.y, ...
                            stim.feedback.colour);
                        
                        DrawFormattedText(window,stim.feedback.error_tryAgain.text,...
                            stim.feedback.error_tryAgain.x, stim.feedback.error_tryAgain.y, ...
                            stim.feedback.colour);
                        [~,onsets.(['fbk_fail_trial_',num2str(iTrial)]).(['fail_',num2str(i_trial_failed)])] = Screen(window,'Flip');
                        WaitSecs(t_fail_and_repeat_fbk);
                    end
                end
                
            case 'mental'
                mentalE_prm.startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen(iTrial))]); % adapt start angle to current level of difficulty
                n_max_to_reach_perTrial(iTrial) = n_to_reach.(['E_level_',num2str(E_chosen(iTrial))]);
                % If participant do not finish effort in time, redo the trial.
                while trial_success == 0 && i_trial_failed <= redo_limit
                    [perfSummary{iTrial},...
                        trial_was_successfull(iTrial),...
                        onsets.effortPeriod{iTrial}] = mental_effort_perf_Nback(scr, stim, key,...
                        mental_nbers_per_trial(iTrial,:),...
                        mentalE_prm, n_max_to_reach_perTrial(iTrial),...
                        'col1', 'noInstructions', timeLimitPerf, t_max_effort,errorLimits);
                    trial_success = trial_was_successfull(iTrial);
                    % Save the data if it was an uncessfull trial, to prevent rewriting of the data
                    if trial_success == 0
                        i_trial_failed = i_trial_failed +1;
                        failed_trials{i_trial_failed}.perfSummary = perfSummary{iTrial};
                        failed_trials{i_trial_failed}.trial_was_successfull = trial_was_successfull(iTrial);
                        failed_trials{i_trial_failed}.onset.effortPeriod =  onsets.effortPeriod{iTrial};
                        failed_trials{i_trial_failed}.i_trial_idx = iTrial;
                        % for the mental effort case where too subject was
                        % too slow
                        DrawFormattedText(window, stim.feedback.error_tooSlow.text,...
                            stim.feedback.error_tooSlow.x, stim.feedback.error_tooSlow.y, ...
                            stim.feedback.colour);
                        
                        DrawFormattedText(window,stim.feedback.error_tryAgain.text,...
                            stim.feedback.error_tryAgain.x, stim.feedback.error_tryAgain.y, ...
                            stim.feedback.colour);
                        [~,onsets.(['fbk_fail_trial_',num2str(iTrial)]).(['fail_',num2str(i_trial_failed)])] = Screen(window,'Flip');
                        WaitSecs(t_fail_and_repeat_fbk);
                    end
                end
        end % effort type loop
        effortTime(iTrial) = toc;
    end % choice made or not?
    
    %% Feedback period
    switch trial_was_successfull(iTrial)
        case 0 % trial failed
            
            % trial failure = too slow (for physical effort always)
            % for mental effort: either too slow or because too many errors
            % => adapt the feedback accordingly
            % display error message
            if choice(iTrial) == 0 ||...
                    strcmp(effort_type,'physical') ||...
                    (strcmp(effort_type,'mental') &&...
                    Ep_or_Em_vars.errorLimits.useOfErrorThreshold == true &&...
                    perfSummary{iTrial}.n_errorsMade < Ep_or_Em_vars.errorLimits.errorThreshold)
                DrawFormattedText(window, stim.feedback.error_tooSlow.text,...
                    stim.feedback.error_tooSlow.x, stim.feedback.error_tooSlow.y, ...
                    stim.feedback.colour);
            elseif strcmp(effort_type,'mental') &&...
                    (Ep_or_Em_vars.errorLimits.useOfErrorThreshold == true &&...
                    perfSummary{iTrial}.n_errorsMade >= Ep_or_Em_vars.errorLimits.errorThreshold) ||...
                    (i_trial_failed > redo_limit)
                % for the mental effort case where too many errors were made,
                % adapt the error feedback accordingly
                DrawFormattedText(window, stim.feedback.error_tooManyErrors.text,...
                    stim.feedback.error_tooManyErrors.x, stim.feedback.error_tooManyErrors.y, ...
                    stim.feedback.colour);
            end
            % display amount of money lost because participant was too
            % slow/did too many mistakes
            DrawFormattedText(window, stim.feedback.error_moneyLoss.text,...
                stim.feedback.error_moneyLoss.x, stim.feedback.error_moneyLoss.y,...
                stim.feedback.error_moneyLoss.colour);
            [~,onsets.fbk_fail(iTrial)] = Screen(window,'Flip');
            
            % record loss for the current trial
            gain(iTrial) = -R_money.trialFail;
            
        case 1 % trial is a success
            % display feedback
            switch R_or_P
                case 'R'
                    DrawFormattedText(window, stim.feedback.reward.text,...
                        stim.feedback.reward.x, stim.feedback.reward.y,...
                        stim.feedback.colour);
                case 'P'
                    DrawFormattedText(window, stim.feedback.punishment.text,...
                        stim.feedback.punishment.x, stim.feedback.punishment.y,...
                        stim.feedback.colour);
            end
            drawRewardAmount(scr, stim, R_chosen(iTrial), R_or_P, 'middle_center_start');
            
            [~,onsets.fbk(iTrial)] = Screen(window,'Flip');
            switch R_or_P
                case 'R'
                    onsets.fbk_win(iTrial) = onsets.fbk(iTrial);
                case 'P'
                    onsets.fbk_loss(iTrial) = onsets.fbk(iTrial);
            end
            
            % record loss for the current trial
            switch R_or_P
                case 'R'
                    gain(iTrial) = R_chosen(iTrial);
                case 'P'
                    gain(iTrial) = (-1)*R_chosen(iTrial);
            end
        otherwise
            error(['weird behavior with trial_was_successfull variable. ',...
                'It appears as if it was not initialized for the current trial.']);
    end
    
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
    disp(['Trial ',num2str(iTrial),'/',num2str(nTrials),' done']);
end % trial loop

%% extract relevant training data
summary.IP = IP;
summary.onsets = onsets;
summary.R_chosen = R_chosen;
summary.E_chosen = E_chosen;
summary.optionChosen = choice;
summary.perfSummary = perfSummary;
summary.gain = gain;
summary.totalGain = totalGain;
summary.trial_was_successfull = trial_was_successfull;
summary.effortTime = effortTime;
summary.failed_trials = failed_trials;
switch effort_type
    case 'mental'
        summary.mentalE_prm = mentalE_prm;
        summary.n_max_to_reach_perTrial = n_max_to_reach_perTrial;
        summary.i_sub = i_sub;
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