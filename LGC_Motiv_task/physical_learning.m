function[perfSummary, onsets] = physical_learning(scr, stim, dq, n_E_levels, Ep_time_levels,...
    F_threshold, F_tolerance, MVC,...
    n_learningForceRepeats, timings)
% [perfSummary, onsets] = physical_learning(scr, stim, dq, n_E_levels, Ep_time_levels,...
%     F_threshold, F_tolerance, MVC,...
%     n_learningForceRepeats, timings)
%physical_learning will perform a short learning for the physical effort
%task. For each level of effort, a few trials will be performed in order to
%ensure that the participant understands what corresponds to each level of
%difficulty.
%
% INPUTS
% scr: structure with screen parameters
%
% stim: structure with stimuli parameters
%
% dq: identification of grip
%
% n_E_levels: number of effort levels
%
% Ep_time_levels: structure with duration corresponding to each level of
% force
%
% F_threshold: force level to reach
%
% F_tolerance: tolerance around the threshold
%
% MVC: maximum voluntary contraction force obtained during the calibration
% process
%
% n_learningForceRepeats: number of repetitions of the learning process
% (how many time they will need to perform each level of force)
%
% timings: structure with relevant timings for learning phase
%
% OUTPUTS
% perfSummary: structure with performance summary variables
%
% onsets: structure with information about timings of the experiment

%% screen parameters
window = scr.window;
%% time parameters
time_limit = false;
t_learning_rest = timings.learning_rest;

n_learningTrials = n_learningForceRepeats*n_E_levels;
%% initialize vars of interst
[perfSummary, onsets.effortPeriod] = deal(cell(1,n_learningTrials));
onsets.learningRest = NaN(1,n_learningTrials);

%% perform learning
jTrial = 0;
for iForceRepeat = 1:n_learningForceRepeats
    for iEffortLevel = 0:n_E_levels-1
        jTrial = jTrial + 1;
        [perfSummary{jTrial},...
            ~,...
            onsets.effortPeriod{jTrial}] = physical_effort_perf(scr, stim, dq,...
            MVC,...
            iEffortLevel,...
            Ep_time_levels,...
            F_threshold, F_tolerance,...
            time_limit, timings);
        
        %% Show a rest text and give some rest
        DrawFormattedText(window, stim.MVC_rest.text,...
            stim.MVC_rest.x, stim.MVC_rest.y, stim.MVC_rest.colour);
        [~,timeNow]  = Screen(window,'Flip');
        onsets.rest(jTrial) = timeNow;
        WaitSecs(t_learning_rest);
        
        %% display number of trials done for the experimenter
        disp(['Physical learning trial ',num2str(jTrial),'/',num2str(n_learningTrials),' done']);
    end % effort level loop
end % loop of learning repetitions
        
end % function