function[summary] = choice_and_perf(scr, stim, key,...
    effort_type, Ep_or_Em_vars, R_money,...
    training_R_P_RP_or_mainTask, nTrials, choiceOptions, timings,...
    results_folder, file_nm)
% [summary] = choice_and_perf(scr, stim, key,...
%     effort_type, Ep_or_Em_vars, R_money,...
%     R_or_P_or_RP_condition, nTrials, choiceOptions, timings,...
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
% money to compute gain within the session
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
yScreenCenter = scr.yCenter;
xScreenCenter = scr.xCenter;
barTimeWaitRect = stim.barTimeWaitRect;

% timings
t_max_effort    = timings.max_effort;
t_cross         = timings.cross.(training_R_P_RP_or_mainTask);
t_choice        = timings.choice;
t_dispChoice    = timings.dispChoice;
t_fbk           = timings.feedback;

% specific variables
switch effort_type
    case 'mental'
        i_sub = Ep_or_Em_vars.i_sub;
        n_to_reach = Ep_or_Em_vars.n_to_reach;
    case 'physical'
        MVC = Ep_or_Em_vars.MVC;
        dq = Ep_or_Em_vars.dq;
        Ep_time_levels = Ep_or_Em_vars.Ep_time_levels;
        F_threshold = Ep_or_Em_vars.F_threshold;
        F_tolerance = Ep_or_Em_vars.F_tolerance;
end


%% initialize onsets
[onsets.cross,...
    onsets.dispChoiceOptions,...
    onsets.choice,...
    onsets.keyReleaseMessage,...
    onsets.cross_after_buttonRelease,...
    onsets.fbk, onsets.fbk_win, onsets.fbk_loss,...
    onsets.fbk_fail,...
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

time_limit = true; % time limit to reach level of force required
for iTrial = 1:nTrials
    
    %% fixation cross period
    Screen('FillRect',window, white, stim.cross.verticalLine); % vertical line
    Screen('FillRect',window, white, stim.cross.horizontalLine); % horizontal line
    [~,onsets.cross(iTrial)] = Screen('Flip',window); % display the cross on screen
    WaitSecs(t_cross(iTrial));
    
    %% check that no key is being pressed before the choice trial starts
    [was_a_key_pressed_bf_trial(iTrial),...
        onsets.keyReleaseMessage(iTrial)] = check_keys_are_up(scr, key);
    
    % if a key was pressed before starting the trial => show the fixation
    % cross again with a similar amount of time
    if was_a_key_pressed_bf_trial(iTrial) == 1
        Screen('FillRect',window,white, stim.cross.verticalLine); % vertical line
        Screen('FillRect',window,white, stim.cross.horizontalLine); % horizontal line
        [~,onsets.cross_after_buttonRelease(iTrial)] = Screen('Flip',window); % display the cross on screen
        WaitSecs(t_cross(iTrial));
    end
    
    %% extract monetary incentive, effort level and reward/punishment condition
    R_left_tmp = choiceOptions.R.left(iTrial);
    R_right_tmp = choiceOptions.R.right(iTrial);
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
                t_choice, key);
        end % keep performing the trial until a choice is made
    else % for actual task, if they don't answer in time, consider the trial as a failure
        [choice(iTrial),...
            onsets.dispChoiceOptions(iTrial),...
            onsets.choice(iTrial),...
            stoptask] = choice_period(scr, stim,...
            R_left_tmp, R_right_tmp, E_left_tmp, E_right_tmp, R_or_P_tmp,...
            t_choice, key);
    end
    
    % extract choice made
    switch choice(iTrial)
        case -1 % choice = left option
            R_chosen(iTrial) = R_left_tmp;
            E_chosen(iTrial) = E_left_tmp;
        case 1 % choice = right option
            R_chosen(iTrial) = R_right_tmp;
            E_chosen(iTrial) = E_right_tmp;
        case 0 % no option was selected
            R_chosen(iTrial) = 0;
            %             E_chosen(iTrial) = max(E_left_tmp, E_right_tmp); % by default opt for the higher effort with no reward if no option was selected
            E_chosen(iTrial) = 0;
    end
    
    %% check if escape was pressed => stop everything if so
    if stoptask == 1
        % save all the data in case you still want to analyze it
        save([results_folder, file_nm,'_earlyEnd_tmp.mat']);
        break;
    end
    
    %% chosen option display period
    [time_dispChoice] = choice_task_dispChosen(scr, stim, R_chosen(iTrial), E_chosen(iTrial), R_or_P_tmp);
    onsets.dispChoice(iTrial) = time_dispChoice;
    WaitSecs(t_dispChoice);
    
    %% Effort period
    % perform the effort if a choice was made, otherwise punish the
    % subject for not answering to the choice
    if choice(iTrial) == 0 % no choice was made => failure
        trial_was_successfull(iTrial) = 0;
        
    elseif ismember(choice(iTrial), [-1,1]) % choice done => perform the corresponding effort
        
        %% mental effort: check no key is being pressed before the start of the effort period
        % for physical effort: useless since only the grip is required
        if strcmp(effort_type,'mental')
            [was_a_key_pressed_bf_trial(iTrial),...
                onsets.keyReleaseMessage(iTrial)] = check_keys_are_up(scr, key);
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
                    time_limit, timings);
            case 'mental'
                mentalE_prm.startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen(iTrial))]); % adapt start angle to current level of difficulty
                n_max_to_reach_perTrial(iTrial) = n_to_reach.(['E_level_',num2str(E_chosen(iTrial))]);
                [perfSummary{iTrial},...
                    trial_was_successfull(iTrial),...
                    onsets.effortPeriod{iTrial}] = mental_effort_perf(scr, stim, key,...
                    mental_nbers_per_trial(iTrial,:),...
                    mentalE_prm, n_max_to_reach_perTrial(iTrial),...
                    'all', 'noInstructions', time_limit, t_max_effort);
        end % effort type loop
        effortTime(iTrial) = toc;
    end % choice made or not?
    
    %% Feedback period
    switch trial_was_successfull(iTrial)
        case 0 % trial failed
            % display message
            %             DrawFormattedText(window,'Too slow!', xScreenCenter, (1/5)*yScreenCenter);
            DrawFormattedText(window,'Trop lent!',...
                'center', (1/5)*yScreenCenter,...
                white);
            % display money loss for failing
            
            [~,onsets.fbk_fail(iTrial)] = Screen(window,'Flip');
            
            % record loss for the current trial
            gain(iTrial) = -R_money.trialFail;
            
        case 1 % trial is a success
            % display feedback
            switch R_or_P_tmp
                case 'R'
                    DrawFormattedText(window,'Vous avez obtenu',...
                        'center', (1/5)*yScreenCenter,...
                        white);
                    fbkSign = '+';
                    fbkColour = stim.reward.text.colour;
                case 'P'
                    DrawFormattedText(window,'Vous avez perdu',...
                        'center', (1/5)*yScreenCenter,...
                        white);
                    fbkSign = '-';
                    fbkColour = stim.punishment.text.colour;
            end
            % display money won/lost
            trialMoneyObtained = sprintf('%0.2f',R_chosen(iTrial));
            DrawFormattedText(window, [fbkSign, trialMoneyObtained,' CHF'],...
                stim.reward.text.middle_center_start(1),...
                stim.reward.text.middle_center_start(2),...
                fbkColour);
            
            [~,onsets.fbk(iTrial)] = Screen(window,'Flip');
            switch R_or_P_tmp
                case 'R'
                    onsets.fbk_win(iTrial) = onsets.fbk(iTrial);
                case 'P'
                    onsets.fbk_loss(iTrial) = onsets.fbk(iTrial);
            end
            
            % record loss for the current trial
            switch R_or_P_tmp
                case 'R'
                    gain(iTrial) = R_money.(['R_',num2str(R_chosen(iTrial))]);
                case 'P'
                    gain(iTrial) = (-1)*R_money.(['R_',num2str(R_chosen(iTrial))]);
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
            DrawFormattedText(window,'Temps restant','center',yScreenCenter*(1/2),white);
            % draw one global fixed rectangle showing the total duration
            Screen('FrameRect',window, black, barTimeWaitRect);
            
            % draw one second rectangle updating dynamically showing the
            % time remaining
            Screen('FillRect',window, black, barTimeWaitRect_bis);
            
            Screen(window,'Flip');
        end % display until time catches up with maximum effort time
    end % if all time taken, no need for time penalty
    
    %% display number of trials done for the experimenter
    disp(['Trial ',num2str(iTrial),'/',num2str(nTrials),' done']);
end % trial loop

%% extract relevant training data
summary.onsets = onsets;
summary.choiceOptions = choiceOptions;
summary.R_chosen = R_chosen;
summary.E_chosen = E_chosen;
summary.perfSummary = perfSummary;
summary.gain = gain;
summary.totalGain = totalGain;
summary.trial_was_successfull = trial_was_successfull;
summary.effortTime = effortTime;
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