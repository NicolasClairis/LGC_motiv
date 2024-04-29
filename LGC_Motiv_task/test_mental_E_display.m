% script to show different effort levels in the physical effort task

[scr, xScreenCenter, yScreenCenter,...
    window, baselineTextSize] = ScreenConfiguration(0, 1);

[trainingTimes, calibTimes, learningTimes, taskTimes, mainTimes] = timings_definition({'RP'}, 54, 4, 'mental');
langage = 'fr';

n_E_levels = 4;
[stim] = stim_initialize(scr, n_E_levels, langage);
key = relevant_key_definition('mental', 0, 4);

% initialize main parameters of the task
mentalE_prm = mental_effort_parameters();
mental_n_col = mentalE_prm.mental_n_col;
sideQuestion = mentalE_prm.sideQuestion;
N_back       = mentalE_prm.Nback;

%% main variables
NmaxCalib = 6;
E_chosen = 3;
n_to_reach = mental_N_answersPerLevel(n_E_levels, NmaxCalib);
n_max_to_reach = n_to_reach.(['E_level_',num2str(E_chosen)]);
R_or_P = 'R';
R_chosen = 0.65;
n_trials = 54;
n_questions = 8;

[mental_nbers_per_trial] = mental_numbers(n_trials);
savePath = 'P:\boulot\postdoc_CarmenSandi\experiment_figures\task\';
%% angle values
mentalE_prm.startAngle = stim.difficulty.startAngle.(['level_',num2str(E_chosen)]); % adapt start angle to current level of difficulty
startAngle_currentTrial = mentalE_prm.startAngle;
endAngle = 360;
totalAngleDistance = endAngle - startAngle_currentTrial;
% coordinates for the angle corresponding to the max until now
maxPerfUntilNowAngle = [];
% coordinates for the angle corresponding to the minimal performance to
% accomplish
minPerfAngle = [];
% define angle
currentAngle = NaN(1,n_questions);
% stable for 2 first replies before N-back starts
currentAngle(1:(1+N_back)) = startAngle_currentTrial;

%% display effort scale on the left for the live display of the effort level
nVectorDisplay = mental_nbers_per_trial(1,:);

for i_question = 1:n_questions
    mental_display_stim(scr, stim,...
        currentAngle(i_question), endAngle,...
        sideQuestion, 1, nVectorDisplay(i_question), mental_n_col,...
        'noInstructions', maxPerfUntilNowAngle, minPerfAngle, R_chosen, R_or_P);
    
    imageArray = Screen('GetImage', window);
    imwrite(imageArray,[savePath,'test_physical_E_display_E',num2str(E_chosen),'_quest',num2str(i_question),'.png']);
    
    if i_question > N_back
        currentAngle(i_question + 1) = currentAngle(i_question) + totalAngleDistance/n_max_to_reach;
    end
    WaitSecs(1);
end
sca;