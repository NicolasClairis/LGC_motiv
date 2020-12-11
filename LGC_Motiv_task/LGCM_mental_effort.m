function[perf, onsetEffortPeriod] = LGCM_mental_effort(scr, stim,...
    t_effort_max, n_correctAnswers,...
    R_chosen, R_or_P, E_chosen, n_E_levels,...
    numberVector, switchPerc, task_start)
% 
% LGCM_mental_effort sub-function of LGCM_choice_task_main.m to perform the
% mental effort task
%
% INPUTS
% scr: structure with screen parameters
%
% stim: stucture with stimuli parameters
% 
% t_effort_max: maximal time allowed to perform the effort
%
% n_correctAnswers: number of correct subsequent answers required for this
% trial
%
% R_chosen: reward level chosen for this trial
%
% R_or_P: reward or punishment trial?
%
% E_chosen: effort level chosen for this trial
%
% n_E_levels: maximum number of effort levels
%
% numberVector: big vector with list of numbers
%
% switchPerc: percentage of switch for this task
%
% task_start: string which indicates which is the first task to do
% (1) start with odd/even task
% (2) higher/lower than 5 task
% 
% OUTPUTS
% perf: vector with the information about subsequent answers to each
% question (0 when wrong and 1 when correct)
%
% onsetEffortPeriod: onset when participant can start performing the effort
%
% onsetEfforts
%
% questions: structure with detail of each question asked in the trial
% .number: vector with the number corresponding to each question asked
% .task_nature: vector with a code corresponding to the nature of the task
% asked at each question ((1) odd/even; (2): higher/lower than 5)

%% extract main variables of interest
window = scr.window;
% angle values
startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen)]);
currentAngle = startAngle;
endAngle = stim.difficulty.arcEndAngle;
totalAngleDistance = endAngle - startAngle;

%% display reward (or punishment) display on the background
switch R_or_P
    case 'R'
        Screen('DrawTexture', window,...
            stim.reward.texture.(['reward_',num2str(R_chosen)]),...
            [],...
            stim.reward.middle_center.(['reward_',num2str(R_chosen)]));
    case 'P'
        error('punishment background to be defined yet');
end

%% display amount of effort to be done on top of the incentive as an arc
Screen('FillArc', window,...
    stim.difficulty.currLevelColor,...
    stim.difficulty.middle_center,...
    startAngle,...
    endAngle - startAngle);

%% display initial information on the screen and record timing
[~,onsetEffortPeriod] = Screen(window,'Flip');

%% start the effort
effortPhaseOver = 0;
questions_solved = 0;

while effortPhaseOver == 0
    %% check timing => task stops when maximal time achieved, 
    % even if the correct amount of answers has not been attained
    timeNow = GetSecs;
    if timeNow > onsetEffortPeriod + t_effort_max
        trial_failed = 1;
        effortPhaseOver = 1;
    end
    
    %% check key press
    [keyisdown, secs, keycode] = KbCheck;
    
    %% check answers provided based on the current task
    if questions_solved < n_correctAnswers
        
        %% when error made => record + re-initialize the counter
        isCorrect_currAnswer = 0;
        questions_solved = 0;
        perf = [perf, isCorrect_currAnswer];
        onsetEfforts.(['effort_',num2str(isCorrect_currAnswer)]
        
        %% when correct answer provided
        isCorrect_currAnswer = 1;
        questions_solved = questions_solved + isCorrect_currAnswer;
        perf = [perf, isCorrect_currAnswer];
        onsetEfforts.(['effort_',num2str(isCorrect_currAnswer)]
    else
        trial_failed = 0;
        effortPhaseOver = 1;
    end
end %

end % function