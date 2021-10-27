function[perfSummary, onsets] = mental_learning(scr, stim, key, n_E_levels, n_to_reach, n_learningRepeats, timings)
%
% mental_learning will make the participant perform nXX repetitions of each
% effort level to see what each effort level means.
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with stimuli parameters
%
% key: structure with key information
%
% n_E_levels: number of effort levels to learn
%
% n_to_reach: structure with number of answers to provide for each effort
% level
%
% n_learningRepeats: number of repetitions of each effort level
%
% timings: structure with relevant timings for the experiment
%
% OUTPUTS
% perfSummary: structure with summary of performance
%
% onsets: sturcture with onsets data

%% screen parameters
window = scr.window;

%% time parameters
time_limit = timings.time_limit;
t_max = timings.t_max;
t_learning_rest = timings.learning_rest;

%% error handling
errorLimits.useOfErrorMapping = false;
errorLimits.useOfErrorThreshold = false;

%% define main parameters
mentalE_prm = mental_effort_parameters();

% define number of trials to perform learning
n_learningTrials = n_learningRepeats*n_E_levels;
% initialize the numbers to be used
numberVector = mental_numbers(n_learningTrials);
learning_instructions = 'noInstructions';

% initialize var of interest
[perfSummary.trialSummary, trialSuccess, onsets.effortPeriod] = deal(cell(1,n_learningTrials));

%% perform learning
jTrial = 0;
for iEffortRepeat = 1:n_learningRepeats
    for iEffortLevel = 1:n_E_levels
        jTrial = jTrial + 1;
        %% extract effort level for the current trial
        n_max_to_reach_trial = n_to_reach.(['E_level_',num2str(iEffortLevel - 1)]);
        %% extract start angle according to effort level of the current trial
        mentalE_prm.startAngle = (1 - (n_max_to_reach_trial/n_to_reach.(['E_level_',num2str(n_E_levels-1)])))*360;

        %% perform the effort
        [perfSummary.trialSummary{jTrial}, trialSuccess, onsets.effortPeriod{jTrial}] = mental_effort_perf_Nback(scr, stim, key,...
            numberVector, mentalE_prm, n_max_to_reach_trial,...
            learning_instructions, time_limit, t_max, errorLimits);

        %% Show a rest text and give some rest
        DrawFormattedText(window, stim.MVC_rest.text,...
            stim.MVC_rest.x, stim.MVC_rest.y, stim.MVC_rest.colour);
        [~,timeNow]  = Screen(window,'Flip');
        onsets.rest(jTrial) = timeNow;
        WaitSecs(t_learning_rest);

        %% display number of trials done for the experimenter
        disp(['Mental learning (2) trial ',num2str(jTrial),'/',num2str(n_learningTrials),' done']);
    end % effort level
end % effort repetition

%% extract relevant data
perfSummary.numberVector = numberVector;
perfSummary.onsets = onsets;
perfSummary.trialSuccess = trialSuccess;
perfSummary.n_learningTrials = n_learningTrials;
perfSummary.n_E_levels = n_E_levels;
perfSummary.n_learningRepeats = n_learningRepeats;
perfSummary.timings = timings;
perfSummary.t_learning_rest = t_learning_rest;

end % function