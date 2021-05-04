function[trainingSummary] = LGCM_choice_and_perf_training(scr, stim, key,...
    effort_type, Ep_or_Em_vars, R_money,...
    R_or_P_or_RP_condition, R_or_P_trials, nTrials, choiceOptions, timings)
% [trainingSummary] = LGCM_choice_and_perf_training(scr, stim, key,...
%     effort_type, Ep_or_Em_vars, R_money,...
%     R_or_P_or_RP_condition, R_or_P_trials, nTrials, choiceOptions, timings)
% LGCM_choice_and_perf_training: script to perform choice and effort
% performance (learning).
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
% R_or_P_or_RP_condition:
% 'R': reward only training
% 'P': punishment only training
% 'RP': reward and punishment training
%
% R_or_P_trials: vector with 'R' and/or 'P' strings to know which condition
% is the current trial
%
% nTrials: number of trials
% 
% choiceOptions: reward and effort level to display for each option (left/right) for
% each trial
%
% timings: structure with information about the relevant timings
%
% OUTPUTS
% trainingSummary: structure with most relevant variables extracted during
% the performance

%% load main paramaters
window = scr.window;
white = scr.colours.white;
black = scr.colours.black;
yScreenCenter = scr.yCenter;
xScreenCenter = scr.xCenter;
barTimeWaitRect = stim.barTimeWaitRect;

% timings
t_max_effort = timings.t_max_effort;
t_cross = timings.t_cross;
t_instructions = timings.t_instructions;
t_choice = timings.t_choice;
t_dispChoice = timings.t_dispChoice;
t_fbk = timings.t_fbk;

% specific variables
switch effort_type
    case 'mental'
        i_sub = Ep_or_Em_vars.i_sub;
        n_to_reach = Ep_or_Em_vars.n_to_reach;
    case 'physical'
        MVC = Ep_or_Em_vars.MVC;
end

%% instruction that main task will start soon
for iTimeLoop = 1:2
    switch R_or_P_or_RP_condition
        case 'R'
            DrawFormattedText(window,...
                ['Vous allez à présent choisir entre deux options associées à différents niveaux de récompense et d''effort '...
                'l''option qui vous paraît la plus intéressante.'],...
                'center', yScreenCenter/3, black, scr.wrapat);
        case 'P'
            DrawFormattedText(window,...
                ['Vous allez à présent choisir entre deux options associées à différents niveaux de pertes et d''effort '...
                'l''option qui vous paraît la moins pénible.'],...
                'center', yScreenCenter/3, black, scr.wrapat);
        case 'RP'
            DrawFormattedText(window,...
                ['Vous allez à présent choisir entre deux options associées à différents niveaux de récompenses ou de pertes et d''effort '...
                'l''option qui vous paraît préférable.'],...
                'center', yScreenCenter/3, black, scr.wrapat);
    end
    if iTimeLoop == 1 % force them to read at first
        [~, onsets.trainingWillStart] = Screen(window, 'Flip');
        WaitSecs(t_instructions);
    elseif iTimeLoop == 2 % after t_instructions seconds, they can manually start
        DrawFormattedText(window, 'Appuyez quand vous êtes prêts à commencer la tâche.',...
            'center', yScreenCenter*15/8, black);
        [~, onsets.trainingWillStart_bis] = Screen(window, 'Flip');
        KbWait;
    end
end % loop over forced reading/manual pass loop

    
    %% initialize onsets for main task
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
            mentalE_prm = LGCM_mental_effort_parameters(i_sub);
            % randomize the order of the numbers appearing on screen
            mental_nbers_per_trial = LGCM_mental_numbers(nTrials);
            
            % randomize the type of the first trial (odd/even or higher/lower
            % than 5)
            mental_taskType_trialStart = LGCM_mental_task_start(nTrials);
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
            onsets.keyReleaseMessage(iTrial)] = LGCM_check_keys_are_up(scr, key);
        
        % if a key was pressed before starting the trial => show the fixation
        % cross again with a similar amount of time
        if was_a_key_pressed_bf_trial(iTrial) == 1
            Screen('FillRect',window,white, stim.cross.verticalLine); % vertical line
            Screen('FillRect',window,white, stim.cross.horizontalLine); % horizontal line
            [~,timeCrossNow] = Screen('Flip',window); % display the cross on screen
            onsets.cross_after_buttonRelease(iTrial) = timeCrossNow;
            WaitSecs(t_cross(iTrial));
        end
        
        %% extract reward of punishment trial condition
        R_or_P_tmp = R_or_P_trials{iTrial};
        
        %% choice period
        while choice(iTrial) == 0
            [choice(iTrial),...
                onsets.dispChoiceOptions(iTrial),...
                onsets.choice(iTrial),...
                stoptask] = LGCM_choice_period(scr, stim,...
                choiceOptions, R_or_P_tmp, iTrial, t_choice, key);
        end % keep performing the trial until a choice is made
        
        %% check if escape was pressed => stop everything if so
        if stoptask == 1
            break;
        end
        
        %% chosen option display period
        [time_dispChoice,...
            R_chosen(iTrial),...
            E_chosen(iTrial)] = LGCM_choice_task_dispChosen(scr, stim, choiceOptions,...
            choice(iTrial), R_or_P_tmp,...
            iTrial);
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
                    onsets.keyReleaseMessage(iTrial)] = LGCM_check_keys_are_up(scr, key);
            end
            
            %% perform the effort
            tic;
            switch effort_type
                case 'physical'
                    [perfSummary{iTrial},...
                        trial_was_successfull(iTrial),...
                        onsets.effortPeriod{iTrial}] = LGCM_physical_effort_perf(scr, stim, dq,...
                        MVC,...
                        E_chosen,...
                        Ep_time_levels,...
                        F_threshold, F_tolerance,...
                        time_limit, t_max_effort);
                case 'mental'
                    mentalE_prm.startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen(iTrial))]); % adapt start angle to current level of difficulty
                    n_max_to_reach_tmp = n_to_reach.(['E_level_',num2str(E_chosen(iTrial))]);
                    [perfSummary{iTrial},...
                        trial_was_successfull(iTrial),...
                        onsets.effortPeriod{iTrial}] = LGCM_mental_effort_perf(scr, stim, key,...
                        mental_nbers_per_trial(iTrial,:),...
                        mentalE_prm, n_max_to_reach_tmp,...
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
                    case 'P'
                        DrawFormattedText(window,'Vous avez perdu',...
                            'center', (1/5)*yScreenCenter,...
                            white);
                end
                
                % display monetary incentive
                R_chosen_tmp = ['reward_',num2str(R_chosen(iTrial))];
                Screen('DrawTexture', window,...
                    stim.reward.texture.(R_chosen_tmp),...
                    [],...
                    stim.chosenOption.reward.(R_chosen_tmp));
                
                % punishments: add negative  overlay on top of the monetary
                % incentive
                if strcmp(R_or_P_tmp,'P')
                    Screen('FillOval', window, stim.punishment.colourOverlay,...
                        stim.punishment.circleOverlay.top_center.(R_chosen_tmp));
                end
                
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
            while toc < (t_max_effort - effortTime(iTrial))
                % update bar with time remaining
                percTimeAchieved = (toc + effortTime(iTrial))./t_max_effort;
                barTimeWaitRect_bis = barTimeWaitRect;
                % start on the right corner of the bar + percentage already
                % achieved and move to the left
                barTimeWaitRect_bis(3) = barTimeWaitRect(3) - percTimeAchieved*(barTimeWaitRect(3) - barTimeWaitRect(1));
                %
                DrawFormattedText(window,'Temps restant','center',yScreenCenter*(1/2));
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
    trainingSummary.onsets = onsets;
    trainingSummary.choiceOptions = choiceOptions;
    trainingSummary.R_or_P = R_or_P_trials;
    trainingSummary.R_chosen = R_chosen;
    trainingSummary.E_chosen = E_chosen;
    trainingSummary.perfSummary = perfSummary;
    trainingSummary.gain = gain;
    trainingSummary.totalGain = totalGain;
    trainingSummary.trial_was_successfull = trial_was_successfull;
    trainingSummary.effortTime = effortTime;
    
end % function