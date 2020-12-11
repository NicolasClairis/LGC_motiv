function[effortDone, perf, onsetEffortPhase, onsetEffort] = LGCM_physical_effort(scr, stim,...
    E_chosen, R_chosen,...
    E_time_levels,...
    F_threshold, F_tolerance,...
    t_effort,...
    R_or_P)
% [effortDone, perf, onsetEffortPhase, onsetEffort] = LGCM_physical_effort(scr, stim,...
%     E_chosen, R_chosen,...
%     E_time_levels,...
%     F_threshold, F_tolerance,...
%     t_effort,...
%     R_or_P)
%
% INPUTS
% scr: structure with screen data
%
% stim: structure with data around the stimuli to display
%
% E_chosen: effort level of the chosen option
%
% R_chosen: reward level of the chosen option (or null if no reward chosen)
%
% E_time_levels: structure which stores the timing to reach for each effort
% level (in seconds)
%
% F_threshold: threshold of force (% of MVC) to maintain
%
% F_tolerance: percentage of MVC tolerance we accept around the threshold
%
% t_effort: maximum time until when the required effort can be performed
%
% R_or_P:
% 'R': reward trial
% 'P': punishment trial
%
% OUTPUTS
% effortDone:
% (0) effort was not achieved => the trial will be repeated with high
% effort without reward
% (1) effort was achieved => reward obtained
%
% perf: detail about the performance (detailled timing + corresponding
% force level)
%
% onsetEffortPhase: timing when the effort starts
%
% onsetEffort: structure with all the timings when the effort started to be
% above the required threshold for the first time
%
% See also LGCM_choice_task_main.m

%% initialize the variables of interest
effortDone = 0;
has_F_threshold_ever_been_reached = 0;
was_F_above_threshold = 0;
% angle values
startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen)]);
currentAngle = startAngle;
endAngle = stim.difficulty.arcEndAngle;
totalAngleDistance = endAngle - startAngle;

% effort time to keep
onsetEffort.starts = [];
t_effort_to_keep = E_time_levels.(['level_',num2str(E_chosen)]);
    
%% display all relevant variables on the screen for the effort period

% display reward/punishment level at the center of the screen
Screen('DrawTexture', window,...
    stim.reward.texture.(['reward_',num2str(R_chosen)]),...
    [],...
    stim.reward.middle_center.(['reward_',num2str(R_chosen)]));

% display effort scale on the left for the live display of the effort level
% add force threshold to maintain

% display difficulty level on top of the reward as an arc
Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
    stim.difficulty.middle_center,...
    startAngle,...
    endAngle - startAngle);

[~,onsetEffortPhase] = Screen('Flip',window);
timeNow = onsetEffortPhase;
    
%% start timer: if they don't achieve the effort in time they will have to repeat without the reward
while timeNow <= onsetEffortPhase + t_effort

    timeNow = GetSecs;
    
    %% read current level of force
    F_now;
    % display effort scale on the left for the live display of the effort level
    % add force threshold to maintain and current level of force being made
    
    %% update the center display according to if force above or below the threshold
    if F_now >= F_threshold - F_tolerance % while force > threshold, update the timer
        
        
        %% if force was not above the threshold
        % => store when the force started above
        if was_F_above_threshold == 0
            %% if force threshold had not been reached yet
            % => update this
            % => store the timing
            if has_F_threshold_ever_been_reached == 0
                has_F_threshold_ever_been_reached = 1;
                onsetEffort.start1 = timeNow;
            else
                n_starts = fieldnames(onset.Effort);
                onsetEffort.(['start',num2str(n_starts + 1)]) = timeNow;
            end
            onsetEffort.starts = [timeNow; onsetEffort.starts];
        end % identify if transition from below to above threshold or not
        
        %% update current angle based on time passed and time remaining until achieving performance
        timeEffortStart_tmp = max(onsetEffort.starts);
        percentageTimeForceAlreadyMaintained = (timeNow - timeEffortStart_tmp)/t_effort_to_keep;
        currentAngle = totalAngleDistance*percentageTimeForceAlreadyMaintained;
        
        %% display principal informations
        % display reward/punishment level at the center of the screen
        Screen('DrawTexture', window,...
            stim.reward.texture.(['reward_',num2str(R_chosen)]),...
            [],...
            stim.reward.middle_center.(['reward_',num2str(R_chosen)]));
        
        % display performance achieved level on top of the reward as an arc
        Screen('FillArc', window,...
            stim.difficulty.currLevelColor,...
            stim.difficulty.middle_center,...
            currentAngle,...
            endAngle - currentAngle);

        [~,timeNow] = Screen('Flip',window);
        
        %% update indicator if F is or not above threshold
        was_F_above_threshold = 1;
        
    else
        %% update angle restart at zero
        currentAngle = startAngle;
        
        %% display
        
        %% update indicator if F is or not above threshold
        was_F_above_threshold = 0;
    end
end % time loop

%% check whether performance was or not achieved
if currentAngle >= endAngle
    effortDone = 1;
end

end % function