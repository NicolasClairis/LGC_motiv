function [] = LGC_onsets_for_fMRI_RL(subject_id)
%LGC_onsets_for_fMRI_RL(subject_id)
% LGC_onsets_for_fMRI_RL extracts the onsets and relevant variables for
% performing the analysis of the reinforcement-learning task in the Motiscan2 2017 march
% study.
%
% INPUTS
% root: path where the study is stored
%
% behaviorDir: path to subject behavioral data
%
% fMRI_Dir: path to subject fMRI data
%
% subid: subject identification number (string)
%
% subject_id: subject full identification name (string)
%
% OUTPUTS
% (learn: structure containing all the relevant infos saved in each subject
% folder)

%% working directories
root = ['E:' filesep 'study2',filesep,'pilots' filesep,...
    'fMRI_pilots',filesep];
subDir = [root subject_id filesep];
behaviorDir = [subDir 'behavior' filesep];
fMRI_Dir = [subDir 'fMRI_analysis' filesep];

%% main infos
if strcmp(subject_id,'fMRI_pilot1_AC')
    subid = '1';
else
    error('other subject than fMRI_pilot1_AC not ready yet');
end

nbLearnRuns = 3; % number of learning runs

% extract specific names and run numbers for each task
% learn_runnames = ls([behaviorDir, filesep, 'MBB_battery_*','*learning_onsets_sub',subid,'_sess*','*.mat']);
learn_runnames = ls([behaviorDir, 'MBB_battery_*learning_onsets_sub',subid,'_sess*.mat']);
learn_runs = nan(1,nbLearnRuns);
for iLearnRun = 1:nbLearnRuns
    learn_runs(iLearnRun) = str2double(learn_runnames(iLearnRun,end-4));
end

% trial_types = {'RL_nonExcludedTrials',...
%     'GL_PairsTrials',...
%     'gainPairsTrials','neutralPairsTrials','lossPairsTrials',...
%     'gainPair_gainFeedbackTrials','gainPair_neutralFeedbackTrials',...
%     'lossPair_neutralFeedbackTrials','lossPair_lossFeedbackTrials'};
trial_periods = {'dispOptions','choice','dispChoice','feedback'};
trial_periods_for_loadedOnset = {'displayOptions','choice','displayChoice','feedback'};
model_names = {'Nico'};
n_mdl_types = length(model_names);

%% main parameters
% trial number different cases
n_learn_trials_per_run = 60;
mid_trialN = n_learn_trials_per_run/2;
n_learn_runs = 3;
trialN_raw = 1:n_learn_trials_per_run;
trialN_zPerRun = zscore( 1:n_learn_trials_per_run );
trialN_zARuns = zscore( 1:(n_learn_trials_per_run*n_learn_runs));
n_G_or_L_trials_per_run = 24;
n_GL_trials_per_run = n_G_or_L_trials_per_run*2;
n_nt_trials_per_run = 12;
mid_G_or_L_trials = n_G_or_L_trials_per_run/2;
mid_GL_trials = n_GL_trials_per_run/2;
mid_nt_trials = n_nt_trials_per_run/2;


%% loop through reinforcement learning runs
jLearnRun = 0; % jLearnRun: index for run learn
for iLearnRun = learn_runs % iLearnRun: real number of learn run (1,4,7)
    sessnber = num2str(iLearnRun);
    
    %% load onsets and duration
    jLearnRun = jLearnRun + 1;
    loadOnsetDur = load([behaviorDir learn_runnames(jLearnRun,:)],'onset','duration');
    onset       = loadOnsetDur.onset;
    duration    = loadOnsetDur.duration;
    % pool all feedbacks together
    [onset.feedback, duration.feedback] = deal( NaN(n_learn_trials_per_run,1) );
    % pool onsets
    % gain PAIR trials
    gainPair_trials_idx = ~isnan(onset.gainFeedback);
    onset.feedback( gainPair_trials_idx )       = onset.gainFeedback(gainPair_trials_idx);
    % neutral PAIR trials
    ntalPair_trials_idx = ~isnan(onset.neutralFeedback); % onset.neutralFeedback will only consider neutral pair trials
    onset.feedback( ntalPair_trials_idx )       = onset.neutralFeedback(ntalPair_trials_idx);
    % loss PAIR trials
    lossPair_trials_idx = ~isnan(onset.lossFeedback);
    onset.feedback( lossPair_trials_idx )       = onset.lossFeedback(lossPair_trials_idx);
    
    % pool durations
    % gain FEEDBACK trials
    gainFbk_trials_idx = ~isnan(duration.gainFeedback);
    duration.feedback( gainFbk_trials_idx )    = duration.gainFeedback(gainFbk_trials_idx);
    % neutral FEEDBACK trials
    ntalFbk_trials_idx = ~isnan(duration.neutralFeedback);
    duration.feedback( ntalFbk_trials_idx )    = duration.neutralFeedback(ntalFbk_trials_idx);
    % loss FEEDBACK trials
    lossFbk_trials_idx = ~isnan(duration.lossFeedback);
    duration.feedback( lossFbk_trials_idx )    = duration.lossFeedback(lossFbk_trials_idx);
    
    % duration choice = from RT to choice display (not done during the
    % task, needs to be added)
    duration.choice = onset.displayChoice - onset.choice;
    % load T0 for this run
    T0 = getfield( load([behaviorDir 'TTL_sub' subid '_sess' sessnber '.mat'],'T0'), 'T0');
    % load behavioral modulators: which pair, which choice, which side
    loadStruct_behavior = load([behaviorDir 'global_sub_' subid '_session_' sessnber '_learning.mat'],...
        'npair','feedback','gain','duration','side','choice'); % ,'response'
    npair       = loadStruct_behavior.npair;
    n_trials    = length(npair);
%     response = loadStruct_behavior.response;
    feedback    = loadStruct_behavior.feedback;
    gain        = loadStruct_behavior.gain; % -1/0/1 according to feedback received
    totalGain   = cumsum(gain.*10); % convert -1/0/1 gains to cumulated euros
    sideBest    = loadStruct_behavior.side;
    choice_LR   = loadStruct_behavior.choice;
    choice_BestOrWorse = (choice_LR == sideBest);
    %% exclude trials with no answer (and eventually also add those where
    % the answer was given too soon by mistake)
    omissionTrials                          = (feedback == 0);
    RL_nonExcludedTrials                    = (feedback ~= 0);
%     n_nonExcludedTrials                     = sum(RL_nonExcludedTrials);
    omissionTrials_idx                      = find(omissionTrials);
    RL_nonExcludedTrials_idx                = find(RL_nonExcludedTrials);
    % identify different type of trials
    gain_idx = npair == 1;
    ntal_idx = npair == 2;
    loss_idx = npair == 3;
    gainPairsTrials                         = (gain_idx.*(RL_nonExcludedTrials) == 1); % "==1" condition necessary to get a logical variable, otherwise rest bugs
    neutralPairsTrials                      = (ntal_idx.*(RL_nonExcludedTrials) == 1);
    lossPairsTrials                         = (loss_idx.*(RL_nonExcludedTrials) == 1);
    %% condition where gain and loss (but not neutral) pairs are pooled has to be extracted
    GL_PairsTrials                          = ((gain_idx | loss_idx).*(RL_nonExcludedTrials) == 1); % gain and loss trials but no neutral pair
    %%
    
    % identify different type of trials
    missedTrials_gainPairsTrials            = (gain_idx.*(omissionTrials) == 1); % "==1" condition necessary to get a logical variable, otherwise rest bugs
    missedTrials_neutralPairsTrials         = (ntal_idx.*(omissionTrials) == 1);
    missedTrials_lossPairsTrials            = (loss_idx.*(omissionTrials) == 1);
    % different type of feedbacks
    gainPair_gainFeedbackTrials         = find((feedback ==  1).*gainPairsTrials); % gain pair - gain feedback
    gainPair_neutralFeedbackTrials      = find((feedback == -1).*gainPairsTrials); % gain pair - neutral feedback
    lossPair_neutralFeedbackTrials      = find((feedback ==  1).*lossPairsTrials); % loss pair - neutral feedback
    lossPair_lossFeedbackTrials         = find((feedback == -1).*lossPairsTrials); % loss pair - loss feedback
    % neutral pair always gives neutral feedback => no need to split
    
%     %% condition where trials are split between first and second half and pooled across conditions is not ready yet
%     RL_firstTrials                  = (RL_nonExcludedTrials == 1).*( (1:length(RL_nonExcludedTrials)) <= mid_trialN);
%     RL_lastTrials                   = (RL_nonExcludedTrials == 1).*( (1:length(RL_nonExcludedTrials)) > mid_trialN);
    %%
    % gain/neutral/loss
    gainPair_firstTrials_idx    = find(gain_idx, mid_G_or_L_trials, 'first');
    gainPair_lastTrials_idx     = find(gain_idx, mid_G_or_L_trials, 'last');
    ntalPair_firstTrials_idx    = find(ntal_idx, mid_nt_trials, 'first');
    ntalPair_lastTrials_idx     = find(ntal_idx, mid_nt_trials, 'last');
    lossPair_firstTrials_idx    = find(loss_idx, mid_G_or_L_trials, 'first');
    lossPair_lastTrials_idx     = find(loss_idx, mid_G_or_L_trials, 'last');
    % same for GL pair
    GLPairs_firstTrials_idx = find( (gain_idx | loss_idx), mid_G_or_L_trials*2, 'first');
    GLPairs_lastTrials_idx = find( (gain_idx | loss_idx), mid_G_or_L_trials*2, 'last');
    
    gainPair_firstTrials        = gainPairsTrials;
    gainPair_firstTrials(gainPair_lastTrials_idx) = false; % remove last trials
    gainPair_lastTrials         = gainPairsTrials;
    gainPair_lastTrials(gainPair_firstTrials_idx) = false; % remove first trials
    ntalPair_firstTrials        = neutralPairsTrials;
    ntalPair_firstTrials(ntalPair_lastTrials_idx) = false; % remove last trials
    ntalPair_lastTrials         = neutralPairsTrials;
    ntalPair_lastTrials(ntalPair_firstTrials_idx) = false; % remove first trials
    lossPair_firstTrials        = lossPairsTrials;
    lossPair_firstTrials(lossPair_lastTrials_idx) = false; % remove last trials
    lossPair_lastTrials         = lossPairsTrials;
    lossPair_lastTrials(lossPair_firstTrials_idx) = false; % remove first trials
    % GL
    GLPairs_firstTrials         = GL_PairsTrials;
    GLPairs_firstTrials(GLPairs_lastTrials_idx) = false;
    GLPairs_lastTrials         = GL_PairsTrials;
    GLPairs_lastTrials(GLPairs_firstTrials_idx) = false;
    
    % apply for gain/loss pair trials
    G_RL_nonExcludedTrials = find(feedback(npair == 1) ~= 0);
    L_RL_nonExcludedTrials = find(feedback(npair == 3) ~= 0);

    
    %% store missed trials
    % gain pairs missed
    for iP = 1:length(trial_periods)
        trial_period_nm = trial_periods{iP};
        trial_period_onset_nm = trial_periods_for_loadedOnset{iP};
        
        if ~isempty(missedTrials_gainPairsTrials)
            % onsets
            learn.onset.missedTrials.gainPair.(trial_period_nm).main    = onset.(trial_period_onset_nm)(missedTrials_gainPairsTrials) - T0;
            % duration
            learn.duration.missedTrials.gainPair.(trial_period_nm).main = duration.(trial_period_onset_nm)(missedTrials_gainPairsTrials);
        else
            % onsets
            learn.onset.missedTrials.gainPair.(trial_period_nm).main    = [];
            % duration
            learn.duration.missedTrials.gainPair.(trial_period_nm).main = [];
        end
        
        % neutral pairs
        if ~isempty(missedTrials_neutralPairsTrials)
            % onsets
            learn.onset.missedTrials.ntalPair.(trial_period_nm).main = onset.(trial_period_onset_nm)(missedTrials_neutralPairsTrials) - T0;
            % durations
            learn.duration.missedTrials.ntalPair.dispOptions.main    = duration.(trial_period_onset_nm)(missedTrials_neutralPairsTrials);
        else
            % onsets
            learn.onset.missedTrials.ntalPair.(trial_period_nm).main = [];
            % durations
            learn.duration.missedTrials.ntalPair.dispOptions.main    = [];
        end
        
        % loss pairs
        if ~isempty(missedTrials_lossPairsTrials)
            % onsets
            learn.onset.missedTrials.lossPair.(trial_period_nm).main    = onset.(trial_period_onset_nm)(missedTrials_lossPairsTrials) - T0;
            % durations
            learn.duration.missedTrials.lossPair.dispOptions.main       = duration.(trial_period_onset_nm)(missedTrials_lossPairsTrials);
        else
            % onsets
            learn.onset.missedTrials.lossPair.(trial_period_nm).main    = [];
            % durations
            learn.duration.missedTrials.lossPair.dispOptions.main       = [];
        end
        
        % all pairs together
        if ~isempty(omissionTrials_idx)
            % onsets
            learn.onset.missedTrials.allPairs.(trial_period_nm).main    = onset.(trial_period_onset_nm)(omissionTrials) - T0;
            % durations
            learn.duration.missedTrials.allPairs.dispOptions.main       = duration.(trial_period_onset_nm)(omissionTrials);
        else
            % onsets
            learn.onset.missedTrials.allPairs.(trial_period_nm).main    = [];
            % durations
            learn.duration.missedTrials.allPairs.dispOptions.main       = [];
        end
        
    end % trial period loop
    
    %% exclude trials of non-interest from the rest of the onsets & durations
    % onsets
    onset.cross_ITI(omissionTrials_idx)             = NaN;
    onset.displayOptions(omissionTrials_idx)        = NaN;
    onset.choice(omissionTrials_idx)                = NaN;
    onset.gainFeedback(omissionTrials_idx)          = NaN;
    onset.neutralFeedback(omissionTrials_idx)       = NaN;
    onset.lossFeedback(omissionTrials_idx)          = NaN;
    onset.displayChoice(omissionTrials_idx)         = NaN;
    % durations
    duration.cross_ITI(omissionTrials_idx)          = NaN;
    duration.displayOptions(omissionTrials_idx)     = NaN;
    duration.choice(omissionTrials_idx)             = NaN;
    duration.gainFeedback(omissionTrials_idx)       = NaN;
    duration.neutralFeedback(omissionTrials_idx)    = NaN;
    duration.lossFeedback(omissionTrials_idx)       = NaN;
    duration.displayChoice(omissionTrials_idx)      = NaN;
    
    %% normalize onsets of interest with T0 and store them into learn
    % ITI
    learn.onset.ITI = onset.cross_ITI - T0;
    
    for iP = 1:length(trial_periods)
        trial_period_nm = trial_periods{iP};
        trial_period_onset_nm = trial_periods_for_loadedOnset{iP};
        
        % gain pairs
        learn.onset.gainPair.(trial_period_nm).main         = onset.(trial_period_onset_nm)(gainPairsTrials)   - T0;
        % gain feedback
        learn.onset.gainPair.(trial_period_nm).GainTrial    = onset.(trial_period_onset_nm)(gainPair_gainFeedbackTrials)    - T0;
        % neutral feedback
        learn.onset.gainPair.(trial_period_nm).NeutralTrial = onset.(trial_period_onset_nm)(gainPair_neutralFeedbackTrials) - T0;
        
        % neutral pairs
        learn.onset.ntalPair.(trial_period_nm).NeutralTrial      = onset.(trial_period_onset_nm)(neutralPairsTrials)        - T0;
        % neutral trial and main are the same for neutral since all neutral
        % images give neutral feedback => avoids any script particular
        % modification
        learn.onset.ntalPair.(trial_period_nm).main              = learn.onset.ntalPair.(trial_period_nm).NeutralTrial;
        
        % loss pairs
        learn.onset.lossPair.(trial_period_nm).main         = onset.(trial_period_onset_nm)(lossPairsTrials)   - T0;
        % loss feedback
        learn.onset.lossPair.(trial_period_nm).LossTrial    = onset.(trial_period_onset_nm)(lossPair_lossFeedbackTrials)   - T0;
        % neutral feedback
        learn.onset.lossPair.(trial_period_nm).NeutralTrial = onset.(trial_period_onset_nm)(lossPair_neutralFeedbackTrials)    - T0;
        
        % pool gain and loss pairs together
        learn.onset.GL_Pairs.(trial_period_nm).main         = onset.(trial_period_onset_nm)(GL_PairsTrials)   - T0;
        
        % all pairs together
        if ismember(trial_period_nm,{'dispOptions','choice','dispChoice'})
            learn.onset.allPairs.(trial_period_nm).main	 = onset.(trial_period_onset_nm)(RL_nonExcludedTrials_idx)  - T0;
        elseif strcmp(trial_period_nm,'feedback')
            learn.onset.allPairs.feedback.main       = [learn.onset.gainPair.feedback.main;...
                learn.onset.ntalPair.feedback.main;...
                learn.onset.lossPair.feedback.main]; % onset of the feedback on screen
            [learn.onset.allPairs.feedback.main, fdb_sort_idx] = sort(learn.onset.allPairs.feedback.main); % re-order by timing ascendant
        else
            error('This should not happen');
        end
        
        % split by trial number
        % gain
        % first trials
        learn.onset.gainPair_first.(trial_period_nm).main         = onset.(trial_period_onset_nm)(gainPair_firstTrials)   - T0;
        % last trials
        learn.onset.gainPair_last.(trial_period_nm).main         = onset.(trial_period_onset_nm)(gainPair_lastTrials)   - T0;
        
        % loss
        % first trials
        learn.onset.lossPair_first.(trial_period_nm).main         = onset.(trial_period_onset_nm)(lossPair_firstTrials)   - T0;
        % last trials
        learn.onset.lossPair_last.(trial_period_nm).main         = onset.(trial_period_onset_nm)(lossPair_lastTrials)   - T0;
        
        % neutral
        % first trials
        learn.onset.ntalPair_first.(trial_period_nm).main         = onset.(trial_period_onset_nm)(ntalPair_firstTrials)   - T0;
        % last trials
        learn.onset.ntalPair_last.(trial_period_nm).main         = onset.(trial_period_onset_nm)(ntalPair_lastTrials)   - T0;
        
        % gain + loss
        % first trials
        learn.onset.GLPairs_first.(trial_period_nm).main         = onset.(trial_period_onset_nm)(GLPairs_firstTrials)   - T0;
        % last trials
        learn.onset.GLPairs_last.(trial_period_nm).main         = onset.(trial_period_onset_nm)(GLPairs_lastTrials)   - T0;
        
        %% extract durations
        learn.duration.ITI = duration.cross_ITI;
        % gain pairs
        learn.duration.gainPair.(trial_period_nm).main         = duration.(trial_period_onset_nm)(gainPairsTrials);
        % gain feedback
        if ismember(trial_period_nm,{'dispOptions','choice','dispChoice'})
            learn.duration.gainPair.(trial_period_nm).GainTrial    = duration.(trial_period_onset_nm)(gainPair_gainFeedbackTrials);
        elseif strcmp(trial_period_nm,'feedback')
            learn.duration.gainPair.feedback.GainTrial       = duration.gainFeedback(gainPair_gainFeedbackTrials);
        end
        % neutral feedback
        if ismember(trial_period_nm,{'dispOptions','choice','dispChoice'})
            learn.duration.gainPair.(trial_period_nm).NeutralTrial = duration.(trial_period_onset_nm)(gainPair_neutralFeedbackTrials);
        elseif strcmp(trial_period_nm,'feedback')
            learn.duration.gainPair.feedback.NeutralTrial	 = duration.neutralFeedback(gainPair_neutralFeedbackTrials);
        end
        
        % neutral pairs
        if ismember(trial_period_nm,{'dispOptions','choice','dispChoice'})
            learn.duration.ntalPair.(trial_period_nm).NeutralTrial      = duration.(trial_period_onset_nm)(neutralPairsTrials);
        elseif strcmp(trial_period_nm,'feedback')
            learn.duration.ntalPair.feedback.NeutralTrial         = duration.neutralFeedback(neutralPairsTrials); % neutral feedback
        end
        % neutral trial and main are the same for neutral since all neutral
        % images give neutral feedback => avoids any script particular
        % modification
        learn.duration.ntalPair.(trial_period_nm).main              = learn.duration.ntalPair.(trial_period_nm).NeutralTrial;
        
        % loss pairs
        learn.duration.lossPair.(trial_period_nm).main         = duration.(trial_period_onset_nm)(lossPairsTrials);
        % loss feedback
        if ismember(trial_period_nm,{'dispOptions','choice','dispChoice'})
            learn.duration.lossPair.(trial_period_nm).LossTrial    = duration.(trial_period_onset_nm)(lossPair_lossFeedbackTrials);
        elseif strcmp(trial_period_nm,'feedback')
            learn.duration.lossPair.feedback.LossTrial       = duration.lossFeedback(lossPair_lossFeedbackTrials);
        end
        % neutral feedback
        if ismember(trial_period_nm,{'dispOptions','choice','dispChoice'})
            learn.duration.lossPair.(trial_period_nm).NeutralTrial = duration.(trial_period_onset_nm)(lossPair_neutralFeedbackTrials);
        elseif strcmp(trial_period_nm,'feedback')
            learn.duration.lossPair.feedback.NeutralTrial    = duration.neutralFeedback(lossPair_neutralFeedbackTrials);
        end
        
        % gain + loss pairs
        learn.duration.GL_Pairs.(trial_period_nm).main         = duration.(trial_period_onset_nm)(GL_PairsTrials);
        
        % all pairs together
        if ismember(trial_period_nm,{'dispOptions','choice','dispChoice'})
            learn.duration.allPairs.(trial_period_nm).main	 = duration.(trial_period_onset_nm)(RL_nonExcludedTrials_idx); % duration display of the pair of options on screen
        elseif strcmp(trial_period_nm,'feedback')
            learn.duration.allPairs.feedback.main       = [learn.duration.gainPair.feedback.main;...
                learn.duration.ntalPair.feedback.main;...
                learn.duration.lossPair.feedback.main]; % duration of the feedback on screen
            learn.duration.allPairs.feedback.main = learn.duration.allPairs.feedback.main(fdb_sort_idx); % re-order by timing ascendant
        end
        
        % split first/last trials
        % gain pairs
        % first
        learn.duration.gainPair_first.(trial_period_nm).main         = duration.(trial_period_onset_nm)(gainPair_firstTrials);
        % last
        learn.duration.gainPair_last.(trial_period_nm).main         = duration.(trial_period_onset_nm)(gainPair_lastTrials);
        
        % loss pairs
        % first
        learn.duration.lossPair_first.(trial_period_nm).main         = duration.(trial_period_onset_nm)(lossPair_firstTrials);
        % last
        learn.duration.lossPair_last.(trial_period_nm).main         = duration.(trial_period_onset_nm)(lossPair_lastTrials);
        
        % neutral pairs
        % first
        learn.duration.ntalPair_first.(trial_period_nm).main         = duration.(trial_period_onset_nm)(ntalPair_firstTrials);
        % last
        learn.duration.ntalPair_last.(trial_period_nm).main         = duration.(trial_period_onset_nm)(ntalPair_lastTrials);
        
        % gain+loss pairs
        % first
        learn.duration.GLPairs_first.(trial_period_nm).main         = duration.(trial_period_onset_nm)(GLPairs_firstTrials);
        % last
        learn.duration.GLPairs_last.(trial_period_nm).main         = duration.(trial_period_onset_nm)(GLPairs_lastTrials);
        
    end % trial period loop
    
    %% extract modulators
    %% extract RT first press = first press moment (choice) - onset of
    % options on screen
    learn.mod.RT_fp.allPairs.main.raw          = learn.onset.allPairs.choice.main               - learn.onset.allPairs.dispOptions.main;
    learn.mod.RT_fp.GL_Pairs.main.raw          = learn.onset.GL_Pairs.choice.main               - learn.onset.GL_Pairs.dispOptions.main;
    learn.mod.RT_fp.gainPair.main.raw          = learn.onset.gainPair.choice.main               - learn.onset.gainPair.dispOptions.main;
    learn.mod.RT_fp.lossPair.main.raw          = learn.onset.lossPair.choice.main               - learn.onset.lossPair.dispOptions.main;
    learn.mod.RT_fp.ntalPair.main.raw          = learn.onset.ntalPair.choice.main               - learn.onset.ntalPair.dispOptions.main;
    learn.mod.RT_fp.gainPair.Gain.raw          = learn.onset.gainPair.choice.GainTrial          - learn.onset.gainPair.dispOptions.GainTrial;
    learn.mod.RT_fp.gainPair.Neutral.raw       = learn.onset.gainPair.choice.NeutralTrial       - learn.onset.gainPair.dispOptions.NeutralTrial;
    learn.mod.RT_fp.ntalPair.Neutral.raw       = learn.onset.ntalPair.choice.NeutralTrial       - learn.onset.ntalPair.dispOptions.NeutralTrial;
    learn.mod.RT_fp.lossPair.Neutral.raw       = learn.onset.lossPair.choice.NeutralTrial       - learn.onset.lossPair.dispOptions.NeutralTrial;
    learn.mod.RT_fp.lossPair.Loss.raw          = learn.onset.lossPair.choice.LossTrial          - learn.onset.lossPair.dispOptions.LossTrial;
    learn.mod.RT_fp.gainPair_first.main.raw    = learn.onset.gainPair_first.choice.main         - learn.onset.gainPair_first.dispOptions.main;
    learn.mod.RT_fp.gainPair_last.main.raw     = learn.onset.gainPair_last.choice.main          - learn.onset.gainPair_last.dispOptions.main;
    learn.mod.RT_fp.lossPair_first.main.raw    = learn.onset.lossPair_first.choice.main         - learn.onset.lossPair_first.dispOptions.main;
    learn.mod.RT_fp.lossPair_last.main.raw     = learn.onset.lossPair_last.choice.main          - learn.onset.lossPair_last.dispOptions.main;
    learn.mod.RT_fp.ntalPair_first.main.raw    = learn.onset.ntalPair_first.choice.main         - learn.onset.ntalPair_first.dispOptions.main;
    learn.mod.RT_fp.ntalPair_last.main.raw     = learn.onset.ntalPair_last.choice.main          - learn.onset.ntalPair_last.dispOptions.main;
    learn.mod.RT_fp.GL_Pairs.main.raw          = learn.onset.GL_Pairs.choice.main               - learn.onset.GL_Pairs.dispOptions.main;
    learn.mod.RT_fp.GLPairs_first.main.raw     = learn.onset.GLPairs_first.choice.main - learn.onset.GLPairs_first.dispOptions.main;
    learn.mod.RT_fp.GLPairs_last.main.raw     = learn.onset.GLPairs_last.choice.main - learn.onset.GLPairs_last.dispOptions.main;
    % RT zscore per run
    learn.mod.RT_fp.allPairs.main.zPerRun       = nanzscore( learn.mod.RT_fp.allPairs.main.raw );
    learn.mod.RT_fp.GL_Pairs.main.zPerRun       = nanzscore( learn.mod.RT_fp.GL_Pairs.main.raw );
    learn.mod.RT_fp.gainPair.main.zPerRun       = nanzscore( learn.mod.RT_fp.gainPair.main.raw );
    learn.mod.RT_fp.lossPair.main.zPerRun       = nanzscore( learn.mod.RT_fp.lossPair.main.raw);
    learn.mod.RT_fp.ntalPair.main.zPerRun       = nanzscore( learn.mod.RT_fp.ntalPair.main.raw);
    learn.mod.RT_fp.gainPair.Gain.zPerRun       = nanzscore( learn.mod.RT_fp.gainPair.Gain.raw);
    learn.mod.RT_fp.gainPair.Neutral.zPerRun    = nanzscore( learn.mod.RT_fp.gainPair.Neutral.raw);
    learn.mod.RT_fp.ntalPair.Neutral.zPerRun    = nanzscore( learn.mod.RT_fp.ntalPair.Neutral.raw);
    learn.mod.RT_fp.lossPair.Neutral.zPerRun    = nanzscore( learn.mod.RT_fp.lossPair.Neutral.raw);
    learn.mod.RT_fp.lossPair.Loss.zPerRun       = nanzscore( learn.mod.RT_fp.lossPair.Loss.raw);
    learn.mod.RT_fp.gainPair_first.main.zPerRun = nanzscore( learn.mod.RT_fp.gainPair_first.main.raw);
    learn.mod.RT_fp.gainPair_last.main.zPerRun  = nanzscore( learn.mod.RT_fp.gainPair_last.main.raw);
    learn.mod.RT_fp.lossPair_first.main.zPerRun = nanzscore( learn.mod.RT_fp.lossPair_first.main.raw);
    learn.mod.RT_fp.lossPair_last.main.zPerRun  = nanzscore( learn.mod.RT_fp.lossPair_last.main.raw);
    learn.mod.RT_fp.ntalPair_first.main.zPerRun = nanzscore( learn.mod.RT_fp.ntalPair_first.main.raw);
    learn.mod.RT_fp.ntalPair_last.main.zPerRun  = nanzscore( learn.mod.RT_fp.ntalPair_last.main.raw);
    learn.mod.RT_fp.GLPairs_first.main.zPerRun     = nanzscore( learn.mod.RT_fp.GLPairs_first.main.raw );
    learn.mod.RT_fp.GLPairs_last.main.zPerRun     = nanzscore( learn.mod.RT_fp.GLPairs_last.main.raw);
    
    %% trial number
    % trial number (raw)
    learn.mod.trialN.raw.allPairs.main                      = trialN_raw(RL_nonExcludedTrials);
    learn.mod.trialN.raw.GL_Pairs.main                      = trialN_raw(GL_PairsTrials);
    learn.mod.trialN.raw.gainPair.main                      = trialN_raw(gainPairsTrials);
    learn.mod.trialN.raw.ntalPair.main                      = trialN_raw(neutralPairsTrials);
    learn.mod.trialN.raw.lossPair.main                      = trialN_raw(lossPairsTrials);
    learn.mod.trialN.raw.gainPair.GainTrial                 = trialN_raw(gainPair_gainFeedbackTrials);
    learn.mod.trialN.raw.gainPair.NeutralTrial              = trialN_raw(gainPair_neutralFeedbackTrials);
    learn.mod.trialN.raw.ntalPair.NeutralTrial              = learn.mod.trialN.raw.ntalPair.main;
    learn.mod.trialN.raw.lossPair.NeutralTrial              = trialN_raw(lossPair_neutralFeedbackTrials);
    learn.mod.trialN.raw.lossPair.LossTrial                 = trialN_raw(lossPair_lossFeedbackTrials);
    learn.mod.trialN.raw.gainPair_first.main                = trialN_raw(gainPair_firstTrials);
    learn.mod.trialN.raw.gainPair_last.main                 = trialN_raw(gainPair_lastTrials);
    learn.mod.trialN.raw.ntalPair_first.main                = trialN_raw(ntalPair_firstTrials);
    learn.mod.trialN.raw.ntalPair_last.main                 = trialN_raw(ntalPair_lastTrials);
    learn.mod.trialN.raw.lossPair_first.main                = trialN_raw(lossPair_firstTrials);
    learn.mod.trialN.raw.lossPair_last.main                 = trialN_raw(lossPair_lastTrials);
    learn.mod.trialN.raw.GLPairs_first.main                 = trialN_raw(GLPairs_firstTrials);
    learn.mod.trialN.raw.GLPairs_last.main                  = trialN_raw(GLPairs_lastTrials);
    % trial number (zscore/run)
    learn.mod.trialN.zPerRun.allPairs.main                  = trialN_zPerRun(RL_nonExcludedTrials);
    learn.mod.trialN.zPerRun.GL_Pairs.main                  = trialN_zPerRun(GL_PairsTrials);
    learn.mod.trialN.zPerRun.gainPair.main                  = trialN_zPerRun(gainPairsTrials);
    learn.mod.trialN.zPerRun.ntalPair.main                  = trialN_zPerRun(neutralPairsTrials);
    learn.mod.trialN.zPerRun.lossPair.main                  = trialN_zPerRun(lossPairsTrials);
    learn.mod.trialN.zPerRun.gainPair.GainTrial             = trialN_zPerRun(gainPair_gainFeedbackTrials);
    learn.mod.trialN.zPerRun.gainPair.NeutralTrial          = trialN_zPerRun(gainPair_neutralFeedbackTrials);
    learn.mod.trialN.zPerRun.ntalPair.NeutralTrial          = learn.mod.trialN.zPerRun.ntalPair.main;
    learn.mod.trialN.zPerRun.lossPair.NeutralTrial          = trialN_zPerRun(lossPair_neutralFeedbackTrials);
    learn.mod.trialN.zPerRun.lossPair.LossTrial             = trialN_zPerRun(lossPair_lossFeedbackTrials);
    learn.mod.trialN.zPerRun.gainPair_first.main            = trialN_zPerRun(gainPair_firstTrials);
    learn.mod.trialN.zPerRun.gainPair_last.main             = trialN_zPerRun(gainPair_lastTrials);
    learn.mod.trialN.zPerRun.ntalPair_first.main            = trialN_zPerRun(ntalPair_firstTrials);
    learn.mod.trialN.zPerRun.ntalPair_last.main             = trialN_zPerRun(ntalPair_lastTrials);
    learn.mod.trialN.zPerRun.lossPair_first.main            = trialN_zPerRun(lossPair_firstTrials);
    learn.mod.trialN.zPerRun.lossPair_last.main             = trialN_zPerRun(lossPair_lastTrials);
    learn.mod.trialN.zPerRun.GLPairs_first.main             = trialN_zPerRun(GLPairs_firstTrials);
    learn.mod.trialN.zPerRun.GLPairs_last.main              = trialN_zPerRun(GLPairs_lastTrials);
    
    % trial number (zscore across runs)
    trials_run_idx = (1:n_learn_trials_per_run) + n_learn_trials_per_run*(jLearnRun - 1);
    trialN_zARun = trialN_zARuns(trials_run_idx);
    learn.mod.trialN.zAcrossRuns.allPairs.main              = trialN_zARun(RL_nonExcludedTrials);
    learn.mod.trialN.zAcrossRuns.GL_Pairs.main              = trialN_zARun(GL_PairsTrials);
    learn.mod.trialN.zAcrossRuns.gainPair.main              = trialN_zARun(gainPairsTrials);
    learn.mod.trialN.zAcrossRuns.ntalPair.main              = trialN_zARun(neutralPairsTrials);
    learn.mod.trialN.zAcrossRuns.lossPair.main              = trialN_zARun(lossPairsTrials);
    learn.mod.trialN.zAcrossRuns.gainPair.GainTrial         = trialN_zARun(gainPair_gainFeedbackTrials);
    learn.mod.trialN.zAcrossRuns.gainPair.NeutralTrial      = trialN_zARun(gainPair_neutralFeedbackTrials);
    learn.mod.trialN.zAcrossRuns.ntalPair.NeutralTrial      = learn.mod.trialN.zAcrossRuns.ntalPair.main;
    learn.mod.trialN.zAcrossRuns.lossPair.NeutralTrial      = trialN_zARun(lossPair_neutralFeedbackTrials);
    learn.mod.trialN.zAcrossRuns.lossPair.LossTrial         = trialN_zARun(lossPair_lossFeedbackTrials);
    learn.mod.trialN.zAcrossRuns.gainPair_first.main        = trialN_zARun(gainPair_firstTrials);
    learn.mod.trialN.zAcrossRuns.gainPair_last.main         = trialN_zARun(gainPair_lastTrials);
    learn.mod.trialN.zAcrossRuns.ntalPair_first.main        = trialN_zARun(ntalPair_firstTrials);
    learn.mod.trialN.zAcrossRuns.ntalPair_last.main         = trialN_zARun(ntalPair_lastTrials);
    learn.mod.trialN.zAcrossRuns.lossPair_first.main        = trialN_zARun(lossPair_firstTrials);
    learn.mod.trialN.zAcrossRuns.lossPair_last.main         = trialN_zARun(lossPair_lastTrials);
    learn.mod.trialN.zAcrossRuns.GLPairs_first.main         = trialN_zARun(GLPairs_firstTrials);
    learn.mod.trialN.zAcrossRuns.GLPairs_last.main          = trialN_zARun(GLPairs_lastTrials);
    
    %% Q.values and related variables from VBA models
    % need to extract feedback for Prediction Error computation
    fbk = feedback;
    fbk(lossPair_neutralFeedbackTrials) = 0;
    fbk(lossPair_lossFeedbackTrials) = -1;
    fbk(neutralPairsTrials) = 0;
    
    % Q.values from VBA models
    for iMdl_type = 1:n_mdl_types
        mdl_nm = model_names{iMdl_type};
        model_nm_learn = [mdl_nm,'_VBA_models'];
        % identify how many models have been made
        RL_model_bis_names = ls([fMRI_Dir,'RLmodel_bis_model*_sub',subid,'.mat']);
        nModels = size(RL_model_bis_names, 1);
        
        for iModel = 6%1:nModels
            
            %% extract the data
            loadStruct_tmp = load([fMRI_Dir, 'RLmodel_bis_model',num2str(iModel),'_sub',subid,'.mat']);
            Qval_tmp = loadStruct_tmp.Qval_sub;
            uncertainty_tmp = loadStruct_tmp.uncertainty_sub;
            run_nm = ['run_',num2str(jLearnRun)];

            % extract Q.value for each cue
            Q_gainPair_gainItem = Qval_tmp.mean.GP_best.(run_nm)';
            Q_gainPair_ntalItem = Qval_tmp.mean.GP_worse.(run_nm)';
            Q_lossPair_ntalItem = Qval_tmp.mean.LP_best.(run_nm)';
            Q_lossPair_lossItem = Qval_tmp.mean.LP_worse.(run_nm)';
            % extract p(choice = worse option) across pairs and
            % trials
            pChoice_worseItem = uncertainty_tmp.pChoice_Worse.(run_nm)';
            pChoice_bestItem = 1 - pChoice_worseItem;

            % extract Q.value across pairs for each trial
            [Q_allPairs_goodCue_wExclTrials,...
                Q_allPairs_badCue_wExclTrials,...
                QchItem_all_wExclTrials, QunchItem_all_wExclTrials,...
                pChoice_leftCue_wExclTrials,...
                pChoice_leftCue_centered_wExclTrials,...
                pChosen_wExclTrials,...
                pUnchosen_wExclTrials,...
                dQ_LR_wExclTrials] = deal(NaN(1,n_trials));
            for iTrial = 1:n_trials

                % extract Q.value according to pair type
                switch npair(iTrial)
                    case 1 % gain pair
                        Q_allPairs_goodCue_wExclTrials(iTrial)  = Q_gainPair_gainItem(iTrial);
                        Q_allPairs_badCue_wExclTrials(iTrial)   = Q_gainPair_ntalItem(iTrial);
                    case 2 % neutral pair
                        Q_allPairs_goodCue_wExclTrials(iTrial)  = 0;
                        Q_allPairs_badCue_wExclTrials(iTrial)   = 0;
                    case 3 % loss pair
                        Q_allPairs_goodCue_wExclTrials(iTrial)  = Q_lossPair_ntalItem(iTrial);
                        Q_allPairs_badCue_wExclTrials(iTrial)   = Q_lossPair_lossItem(iTrial);
                end

                % extract Q.value of the chosen option
                switch choice_BestOrWorse(iTrial)
                    case 1
                        QchItem_all_wExclTrials(iTrial)     = Q_allPairs_goodCue_wExclTrials(iTrial);
                        QunchItem_all_wExclTrials(iTrial)   = Q_allPairs_badCue_wExclTrials(iTrial);
                    case 0
                        QchItem_all_wExclTrials(iTrial)     = Q_allPairs_badCue_wExclTrials(iTrial);
                        QunchItem_all_wExclTrials(iTrial)   = Q_allPairs_goodCue_wExclTrials(iTrial);
                end

                % extract p(choice) and Qleft-Qright according to choice side
                switch sideBest(iTrial)
                    case -1 % best option on the left
                        pChoice_leftCue_wExclTrials(iTrial) = pChoice_bestItem(iTrial);
                        dQ_LR_wExclTrials(iTrial) =...
                            Q_allPairs_goodCue_wExclTrials(iTrial) - Q_allPairs_badCue_wExclTrials(iTrial);
                    case 1 % best option on the right
                        pChoice_leftCue_wExclTrials(iTrial) = pChoice_worseItem(iTrial);
                        dQ_LR_wExclTrials(iTrial) =...
                            Q_allPairs_badCue_wExclTrials(iTrial) - Q_allPairs_goodCue_wExclTrials(iTrial);
                end

                % extract p(chosen)
                switch choice_LR(iTrial)
                    case -1 % left option chosen
                        pChosen_wExclTrials(iTrial)     = pChoice_leftCue_wExclTrials(iTrial);
                        pUnchosen_wExclTrials(iTrial)   = 1-pChoice_leftCue_wExclTrials(iTrial);
                    case 1
                        pChosen_wExclTrials(iTrial)     = 1-pChoice_leftCue_wExclTrials(iTrial);
                        pUnchosen_wExclTrials(iTrial)   = pChoice_leftCue_wExclTrials(iTrial);
                end
            end % trial loop
            dQ_all = Q_allPairs_goodCue_wExclTrials - Q_allPairs_badCue_wExclTrials;

            % filter for bad trials
            % Q.values
            Q_gainPair_gainItem_all = Q_gainPair_gainItem(RL_nonExcludedTrials);
            Q_gainPair_ntalItem_all = Q_gainPair_ntalItem(RL_nonExcludedTrials);
            Q_lossPair_ntalItem_all = Q_lossPair_ntalItem(RL_nonExcludedTrials);
            Q_lossPair_lossItem_all = Q_lossPair_lossItem(RL_nonExcludedTrials);
            Q_allPairs_goodItem_all = Q_allPairs_goodCue_wExclTrials(RL_nonExcludedTrials);
            Q_allPairs_badItem_all  = Q_allPairs_badCue_wExclTrials(RL_nonExcludedTrials);
            Q_GL_Pairs_goodItem_all = Q_allPairs_goodCue_wExclTrials(GL_PairsTrials);
            Q_GL_Pairs_badItem_all  = Q_allPairs_badCue_wExclTrials(GL_PairsTrials);
            Q_GL_Pairs_goodItem_first = Q_allPairs_goodCue_wExclTrials(GLPairs_firstTrials);
            Q_GL_Pairs_goodItem_last = Q_allPairs_goodCue_wExclTrials(GLPairs_lastTrials);
            Q_GL_Pairs_badItem_first  = Q_allPairs_badCue_wExclTrials(GLPairs_firstTrials);
            Q_GL_Pairs_badItem_last  = Q_allPairs_badCue_wExclTrials(GLPairs_lastTrials);
            % p(choice)
            pChoice_worseItem_all = pChoice_worseItem(RL_nonExcludedTrials);
            pChoice_bestItem_all = pChoice_bestItem(RL_nonExcludedTrials);
            pChoice_leftItem_all = pChoice_leftCue_wExclTrials(RL_nonExcludedTrials);
            pChosen_all = pChosen_wExclTrials(RL_nonExcludedTrials);
            pUnchosen_all = pUnchosen_wExclTrials(RL_nonExcludedTrials);
            pChoice_bestItem_GL_Pairs = pChoice_bestItem(GL_PairsTrials);
            pChoice_leftItem_GL_Pairs = pChoice_leftCue_wExclTrials(GL_PairsTrials);
            pChosen_GL_Pairs = pChosen_wExclTrials(GL_PairsTrials);
            pUnchosen_GL_Pairs = pUnchosen_wExclTrials(GL_PairsTrials);
            pChoice_bestItem_GL_Pairs_first = pChoice_bestItem(GLPairs_firstTrials);
            pChoice_bestItem_GL_Pairs_last = pChoice_bestItem(GLPairs_lastTrials);
            pChoice_leftItem_GL_Pairs_first = pChoice_leftCue_wExclTrials(GLPairs_firstTrials);
            pChoice_leftItem_GL_Pairs_last = pChoice_leftCue_wExclTrials(GLPairs_lastTrials);
            pChoice_leftItem_centered_all = pChoice_leftCue_wExclTrials(RL_nonExcludedTrials) - nanmean(pChoice_leftCue_wExclTrials(RL_nonExcludedTrials));
            pChoice_leftItem_centered_GL_Pairs = pChoice_leftCue_wExclTrials(GL_PairsTrials) - nanmean(pChoice_leftCue_wExclTrials(GL_PairsTrials));
            pChoice_leftItem_centered_GL_Pairs_first = pChoice_leftCue_wExclTrials(GLPairs_firstTrials) - nanmean(pChoice_leftCue_wExclTrials(GLPairs_firstTrials));
            pChoice_leftItem_centered_GL_Pairs_last = pChoice_leftCue_wExclTrials(GLPairs_lastTrials) - nanmean(pChoice_leftCue_wExclTrials(GLPairs_lastTrials));
            pChosen_GL_Pairs_first = pChosen_wExclTrials(GLPairs_firstTrials);
            pChosen_GL_Pairs_last = pChosen_wExclTrials(GLPairs_lastTrials);
            pUnchosen_GL_Pairs_first = pUnchosen_wExclTrials(GLPairs_firstTrials);
            pUnchosen_GL_Pairs_last = pUnchosen_wExclTrials(GLPairs_lastTrials);

            % additional uncertainty measures for those models
            if strcmp(mdl_nm,'Nico') % impossible to compute for 'Nico_avg_prm'
                uncertainty_sigmaOpt_all = uncertainty_tmp.sigma_options.(run_nm)(RL_nonExcludedTrials)';
                uncertainty_overlap_distrib_all = uncertainty_tmp.overlap_distrib_options.(run_nm)(RL_nonExcludedTrials)';
            end
            
            % proba for right item
            pChoice_leftItem = pChoice_leftCue_wExclTrials;
            pChoice_rightItem = 1 - pChoice_leftItem;
            pChoice_worseItem_GL_Pairs = 1 - pChoice_bestItem_GL_Pairs; % proba of choosing worse option = 1 - proba of choosing the best
            pChoice_rightItem_all   = 1 - pChoice_leftItem_all;
            pChoice_rightItem_GL_Pairs = 1 - pChoice_leftItem_GL_Pairs;
            % prediction error
            PE_all                  = fbk - QchItem_all_wExclTrials; % both fbk and QchItem_all_wExclTrials should vary between -1 and +1 in each model
            % Q.value
            QchosenItem_all = QchItem_all_wExclTrials(RL_nonExcludedTrials);
            QunchosenItem_all = QunchItem_all_wExclTrials(RL_nonExcludedTrials);
            QchosenItem_GL_PairsTrials = QchItem_all_wExclTrials(GL_PairsTrials);
            QunchosenItem_GL_PairsTrials = QunchItem_all_wExclTrials(GL_PairsTrials);
            % split first/second half
            pChoice_worseItem_GL_Pairs_first = 1 - pChoice_bestItem_GL_Pairs_first; % proba of choosing worse option = 1 - proba of choosing the best
            pChoice_worseItem_GL_Pairs_last = 1 - pChoice_bestItem_GL_Pairs_last; % proba of choosing worse option = 1 - proba of choosing the best
            pChoice_rightItem_GL_Pairs_first = 1 - pChoice_leftItem_GL_Pairs_first;
            pChoice_rightItem_GL_Pairs_last = 1 - pChoice_leftItem_GL_Pairs_last;
            QchosenItem_GL_PairsTrials_first = QchItem_all_wExclTrials(GLPairs_firstTrials);
            QchosenItem_GL_PairsTrials_last = QchItem_all_wExclTrials(GLPairs_lastTrials);
            QunchosenItem_GL_PairsTrials_first = QunchItem_all_wExclTrials(GLPairs_firstTrials);
            QunchosenItem_GL_PairsTrials_last = QunchItem_all_wExclTrials(GLPairs_lastTrials);
            
            
            %% save the data
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.allPairsTrials.GainItem            = Q_gainPair_gainItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.allPairsTrials.bestItem            = Q_gainPair_gainItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.allPairsTrials.NeutralItem         = Q_gainPair_ntalItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.allPairsTrials.worseItem           = Q_gainPair_ntalItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.allPairsTrials.NeutralItem         = Q_lossPair_ntalItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.allPairsTrials.bestItem            = Q_lossPair_ntalItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.allPairsTrials.LossItem            = Q_lossPair_lossItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.allPairsTrials.worseItem           = Q_lossPair_lossItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.allPairs.allPairsTrials.bestItem            = Q_allPairs_goodItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.allPairs.allPairsTrials.worseItem           = Q_allPairs_badItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.allPairs.allPairsTrials.chosenItem          = QchosenItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.allPairs.allPairsTrials.unchosenItem        = QunchosenItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GL_Pairs.GL_PairsTrials.bestItem            = Q_GL_Pairs_goodItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GL_Pairs.GL_PairsTrials.worseItem           = Q_GL_Pairs_badItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GL_Pairs.GL_PairsTrials.chosenItem          = QchosenItem_GL_PairsTrials;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GL_Pairs.GL_PairsTrials.unchosenItem        = QunchosenItem_GL_PairsTrials;
            % split first/second half
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GLPairs_first.GLPairs_firstTrials.bestItem            = Q_GL_Pairs_goodItem_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GLPairs_last.GLPairs_lastTrials.bestItem            = Q_GL_Pairs_goodItem_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GLPairs_first.GLPairs_firstTrials.worseItem           = Q_GL_Pairs_badItem_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GLPairs_last.GLPairs_lastTrials.worseItem           = Q_GL_Pairs_badItem_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GLPairs_first.GLPairs_firstTrials.chosenItem          = QchItem_all_wExclTrials(GLPairs_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.GLPairs_last.GLPairs_lastTrials.unchosenItem        = QchItem_all_wExclTrials(GLPairs_lastTrials);
            
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.allPairs.allPairsTrials = dQ_all(RL_nonExcludedTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.allPairs.allPairsTrials = dQ_LR_wExclTrials(RL_nonExcludedTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.allPairs.allPairsTrials = PE_all(RL_nonExcludedTrials);
            % pool only gain and loss pairs
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.GL_Pairs.GL_PairsTrials = dQ_all(GL_PairsTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.GL_Pairs.GL_PairsTrials = dQ_LR_wExclTrials(GL_PairsTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.GL_Pairs.GL_PairsTrials = PE_all(GL_PairsTrials);
            % split first/last
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.GLPairs_first.GLPairs_firstTrials = PE_all(GLPairs_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.GLPairs_last.GLPairs_lastTrials = PE_all(GLPairs_lastTrials);
            
            % proba choice best/worse
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.allPairs.allPairsTrials   = pChoice_bestItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.allPairs.allPairsTrials  = pChoice_worseItem_all;
            % left/right
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.allPairs.allPairsTrials   = pChoice_leftItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.allPairs.allPairsTrials  = pChoice_rightItem_all;
            % p(choice) -m(p(choice))
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice_centered.left.allPairs.allPairsTrials   = pChoice_leftItem_centered_all;
            % p(chosen)/p(unchosen)
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.chosen.allPairs.allPairsTrials   = pChosen_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.unchosen.allPairs.allPairsTrials  = pUnchosen_all;
            
            % proba choice best/worse
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.GL_Pairs.GL_PairsTrials   = pChoice_bestItem_GL_Pairs;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.GL_Pairs.GL_PairsTrials  = pChoice_worseItem_GL_Pairs;
            % left/right
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.GL_Pairs.GL_PairsTrials   = pChoice_leftItem_GL_Pairs;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.GL_Pairs.GL_PairsTrials  = pChoice_rightItem_GL_Pairs;
            % p(chosen)/p(unchosen)
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.chosen.GL_Pairs.GL_PairsTrials   = pChosen_GL_Pairs;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.unchosen.GL_Pairs.GL_PairsTrials  = pUnchosen_GL_Pairs;
            
            % split first/second half
            % proba choice best/worse
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.GLPairs_first.GLPairs_firstTrials   = pChoice_worseItem_GL_Pairs_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.GLPairs_last.GLPairs_lastTrials   = pChoice_bestItem_GL_Pairs_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.GLPairs_first.GLPairs_firstTrials  = pChoice_worseItem_GL_Pairs_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.GLPairs_last.GLPairs_lastTrials  = pChoice_worseItem_GL_Pairs_last;
            % left/right
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.GLPairs_first.GLPairs_firstTrials   = pChoice_leftItem_GL_Pairs_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.GLPairs_last.GLPairs_lastTrials   = pChoice_leftItem_GL_Pairs_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.GLPairs_first.GLPairs_firstTrials  = pChoice_rightItem_GL_Pairs_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.GLPairs_last.GLPairs_lastTrials  = pChoice_rightItem_GL_Pairs_last;
            % p(chosen)/p(unchosen)
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.chosen.GLPairs_first.GLPairs_firstTrials   = pChosen_GL_Pairs_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.chosen.GLPairs_last.GLPairs_lastTrials   = pChosen_GL_Pairs_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.unchosen.GLPairs_first.GLPairs_firstTrials  = pUnchosen_GL_Pairs_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.unchosen.GLPairs_last.GLPairs_lastTrials  = pUnchosen_GL_Pairs_last;
            
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.allPairs.allPairsTrials   =...
                pChosen_all.*QchosenItem_all + pUnchosen_all.*QunchosenItem_all;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.GL_Pairs.GL_PairsTrials   =...
                pChosen_GL_Pairs.*QchosenItem_GL_PairsTrials + pUnchosen_GL_Pairs.*QunchosenItem_GL_PairsTrials;
            % first/second half
             learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.GLPairs_first.GLPairs_firstTrials   =...
                pChosen_GL_Pairs_first.*QchosenItem_GL_PairsTrials_first + pUnchosen_GL_Pairs_first.*QunchosenItem_GL_PairsTrials_first;
             learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.GLPairs_last.GLPairs_lastTrials   =...
                pChosen_GL_Pairs_last.*QchosenItem_GL_PairsTrials_last + pUnchosen_GL_Pairs_last.*QunchosenItem_GL_PairsTrials_last;
            
            % classify according to pair type
            Q_gainPair_gainItem_gainPair = Q_gainPair_gainItem(gainPairsTrials);
            Q_gainPair_ntalItem_gainPair = Q_gainPair_ntalItem(gainPairsTrials);
            Q_gainPair_chItem_gainPair   = QchItem_all_wExclTrials(gainPairsTrials);
            Q_gainPair_unchItem_gainPair = QunchItem_all_wExclTrials(gainPairsTrials);
            Q_lossPair_ntalItem_lossPair = Q_lossPair_ntalItem(lossPairsTrials);
            Q_lossPair_lossItem_lossPair = Q_lossPair_lossItem(lossPairsTrials);
            Q_lossPair_chItem_lossPair   = QchItem_all_wExclTrials(lossPairsTrials);
            Q_lossPair_unchItem_lossPair = QunchItem_all_wExclTrials(lossPairsTrials);
            % classify according to pair type and split first/second half trials
            Q_gainPair_gainItem_gainPair_first  = Q_gainPair_gainItem(gainPair_firstTrials);
            Q_gainPair_gainItem_gainPair_last   = Q_gainPair_gainItem(gainPair_lastTrials);
            Q_gainPair_ntalItem_gainPair_first  = Q_gainPair_ntalItem(gainPair_firstTrials);
            Q_gainPair_ntalItem_gainPair_last   = Q_gainPair_ntalItem(gainPair_lastTrials);
            Q_gainPair_chItem_gainPair_first    = QchItem_all_wExclTrials(gainPair_firstTrials);
            Q_gainPair_chItem_gainPair_last     = QchItem_all_wExclTrials(gainPair_lastTrials);
            Q_gainPair_unchItem_gainPair_first  = QunchItem_all_wExclTrials(gainPair_firstTrials);
            Q_gainPair_unchItem_gainPair_last   = QunchItem_all_wExclTrials(gainPair_lastTrials);
            Q_lossPair_ntalItem_lossPair_first  = Q_lossPair_ntalItem(lossPair_firstTrials);
            Q_lossPair_ntalItem_lossPair_last   = Q_lossPair_ntalItem(lossPair_lastTrials);
            Q_lossPair_lossItem_lossPair_first  = Q_lossPair_lossItem(lossPair_firstTrials);
            Q_lossPair_lossItem_lossPair_last   = Q_lossPair_lossItem(lossPair_lastTrials);
            Q_lossPair_chItem_lossPair_first    = QchItem_all_wExclTrials(lossPair_firstTrials);
            Q_lossPair_chItem_lossPair_last     = QchItem_all_wExclTrials(lossPair_lastTrials);
            Q_lossPair_unchItem_lossPair_first  = QunchItem_all_wExclTrials(lossPair_firstTrials);
            Q_lossPair_unchItem_lossPair_last   = QunchItem_all_wExclTrials(lossPair_lastTrials);
            
            % best/worse
            pChoice_gainPair_gainItem_gainPair = pChoice_bestItem(gainPairsTrials);
            pChoice_gainPair_ntalItem_gainPair = pChoice_worseItem(gainPairsTrials); % proba of choosing worse option = 1 - proba of choosing the best
            pChoice_lossPair_ntalItem_lossPair = pChoice_bestItem(lossPairsTrials);
            pChoice_lossPair_lossItem_lossPair = pChoice_worseItem(lossPairsTrials);
            % best/worse first/second half trials
            pChoice_gainPair_gainItem_gainPair_first    = pChoice_bestItem(gainPair_firstTrials);
            pChoice_gainPair_gainItem_gainPair_last     = pChoice_bestItem(gainPair_lastTrials);
            pChoice_gainPair_ntalItem_gainPair_first    = pChoice_worseItem(gainPair_firstTrials);
            pChoice_gainPair_ntalItem_gainPair_last     = pChoice_worseItem(gainPair_lastTrials);
            pChoice_lossPair_ntalItem_lossPair_first    = pChoice_bestItem(lossPair_firstTrials);
            pChoice_lossPair_ntalItem_lossPair_last     = pChoice_bestItem(lossPair_lastTrials);
            pChoice_lossPair_lossItem_lossPair_first    = pChoice_worseItem(lossPair_firstTrials);
            pChoice_lossPair_lossItem_lossPair_last     = pChoice_worseItem(lossPair_lastTrials);
            
            % left/right
            pChoice_gainPair_leftItem_gainPair  = pChoice_leftItem(gainPairsTrials);
            pChoice_gainPair_rightItem_gainPair = pChoice_rightItem(gainPairsTrials); % proba of choosing worse option = 1 - proba of choosing the best
            pChoice_lossPair_leftItem_lossPair  = pChoice_leftItem(lossPairsTrials);
            pChoice_lossPair_rightItem_lossPair = pChoice_rightItem(lossPairsTrials);
            % left/right first/second half of trials
            pChoice_gainPair_leftItem_gainPair_first    = pChoice_leftItem(gainPair_firstTrials);
            pChoice_gainPair_leftItem_gainPair_last     = pChoice_leftItem(gainPair_lastTrials);
            pChoice_gainPair_rightItem_gainPair_first   = pChoice_rightItem(gainPair_firstTrials);
            pChoice_gainPair_rightItem_gainPair_last    = pChoice_rightItem(gainPair_lastTrials);
            pChoice_lossPair_leftItem_lossPair_first    = pChoice_leftItem(lossPair_firstTrials);
            pChoice_lossPair_leftItem_lossPair_last     = pChoice_leftItem(lossPair_lastTrials);
            pChoice_lossPair_rightItem_lossPair_first   = pChoice_rightItem(lossPair_firstTrials);
            pChoice_lossPair_rightItem_lossPair_last    = pChoice_rightItem(lossPair_lastTrials);
            
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.gainPairTrials.GainItem       = Q_gainPair_gainItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.gainPairTrials.bestItem       = Q_gainPair_gainItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.gainPairTrials.NeutralItem    = Q_gainPair_ntalItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.gainPairTrials.worseItem      = Q_gainPair_ntalItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.gainPairTrials.chosenItem     = Q_gainPair_chItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.gainPairTrials.unchosenItem   = Q_gainPair_unchItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.lossPairTrials.NeutralItem    = Q_lossPair_ntalItem_lossPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.lossPairTrials.bestItem       = Q_lossPair_ntalItem_lossPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.lossPairTrials.LossItem       = Q_lossPair_lossItem_lossPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.lossPairTrials.worseItem      = Q_lossPair_lossItem_lossPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.lossPairTrials.chosenItem     = Q_lossPair_chItem_lossPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.lossPairTrials.unchosenItem   = Q_lossPair_unchItem_lossPair;
            % split first/second half
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_first.gainPairTrials.GainItem      = Q_gainPair_gainItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_last.gainPairTrials.GainItem       = Q_gainPair_gainItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_first.gainPairTrials.bestItem      = Q_gainPair_gainItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_last.gainPairTrials.bestItem       = Q_gainPair_gainItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_first.gainPairTrials.NeutralItem   = Q_gainPair_ntalItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_last.gainPairTrials.NeutralItem    = Q_gainPair_ntalItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_first.gainPairTrials.worseItem     = Q_gainPair_ntalItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_last.gainPairTrials.worseItem      = Q_gainPair_ntalItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_first.gainPairTrials.chosenItem    = Q_gainPair_chItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_last.gainPairTrials.chosenItem     = Q_gainPair_chItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_first.gainPairTrials.unchosenItem  = Q_gainPair_unchItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair_last.gainPairTrials.unchosenItem   = Q_gainPair_unchItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_first.lossPairTrials.NeutralItem   = Q_lossPair_ntalItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_last.lossPairTrials.NeutralItem    = Q_lossPair_ntalItem_lossPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_first.lossPairTrials.bestItem      = Q_lossPair_ntalItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_last.lossPairTrials.bestItem       = Q_lossPair_ntalItem_lossPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_first.lossPairTrials.LossItem      = Q_lossPair_lossItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_last.lossPairTrials.LossItem       = Q_lossPair_lossItem_lossPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_first.lossPairTrials.worseItem     = Q_lossPair_lossItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_last.lossPairTrials.worseItem      = Q_lossPair_lossItem_lossPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_first.lossPairTrials.chosenItem    = Q_lossPair_chItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_last.lossPairTrials.chosenItem     = Q_lossPair_chItem_lossPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_first.lossPairTrials.unchosenItem  = Q_lossPair_unchItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair_last.lossPairTrials.unchosenItem   = Q_lossPair_unchItem_lossPair_last;
            
            % dQ best-worse
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.gainPair.gainPairTrials             = dQ_all(gainPairsTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.lossPair.lossPairTrials             = dQ_all(lossPairsTrials);
            % dQ left-right
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.gainPair.gainPairTrials             = dQ_LR_wExclTrials(gainPairsTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.lossPair.lossPairTrials             = dQ_LR_wExclTrials(lossPairsTrials);
            % split first/second half
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.gainPair_first.gainPairTrials             = dQ_all(gainPair_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.gainPair_last.gainPairTrials             = dQ_all(gainPair_lastTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.lossPair_first.lossPairTrials             = dQ_all(lossPair_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.lossPair_last.lossPairTrials             = dQ_all(lossPair_lastTrials);

            % dQ left - right
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.gainPair_first.gainPairTrials             = dQ_LR_wExclTrials(gainPair_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.gainPair_last.gainPairTrials             = dQ_LR_wExclTrials(gainPair_lastTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.lossPair_first.lossPairTrials             = dQ_LR_wExclTrials(lossPair_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.lossPair_last.lossPairTrials             = dQ_LR_wExclTrials(lossPair_lastTrials);
            
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.gainPair.gainPairTrials             = PE_all(gainPairsTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.ntalPair.ntalPairTrials             = PE_all(neutralPairsTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.lossPair.lossPairTrials             = PE_all(lossPairsTrials);
            % split first/second half
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.gainPair_first.gainPairTrials            = PE_all(gainPair_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.gainPair_last.gainPairTrials             = PE_all(gainPair_lastTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.ntalPair_first.ntalPairTrials            = PE_all(ntalPair_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.ntalPair_last.ntalPairTrials             = PE_all(ntalPair_lastTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.lossPair_first.lossPairTrials            = PE_all(lossPair_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.lossPair_last.lossPairTrials             = PE_all(lossPair_lastTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.GLPairs_first.GLPairs_firstTrials        = PE_all(GLPairs_firstTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.GLPairs_last.GLPairs_lastTrials          = PE_all(GLPairs_lastTrials);
            
            % split based on best/worse
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.gainPair.gainPairTrials   = pChoice_gainPair_gainItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.gainPair.gainPairTrials  = pChoice_gainPair_ntalItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.lossPair.lossPairTrials   = pChoice_lossPair_ntalItem_lossPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.lossPair.lossPairTrials  = pChoice_lossPair_lossItem_lossPair;
            % split based on left/right
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.gainPair.gainPairTrials   = pChoice_gainPair_leftItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.gainPair.gainPairTrials  = pChoice_gainPair_rightItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.lossPair.lossPairTrials   = pChoice_lossPair_leftItem_lossPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.lossPair.lossPairTrials  = pChoice_lossPair_rightItem_lossPair;
            % split best/worse based on first/second half
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.gainPair_first.gainPairTrials   = pChoice_gainPair_gainItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.gainPair_last.gainPairTrials    = pChoice_gainPair_gainItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.gainPair_first.gainPairTrials  = pChoice_gainPair_ntalItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.gainPair_last.gainPairTrials   = pChoice_gainPair_ntalItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.lossPair_first.lossPairTrials   = pChoice_lossPair_ntalItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.lossPair_last.lossPairTrials    = pChoice_lossPair_ntalItem_lossPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.lossPair_first.lossPairTrials  = pChoice_lossPair_lossItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.lossPair_last.lossPairTrials   = pChoice_lossPair_lossItem_lossPair_last;
            % split left/right based on first/second half
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.gainPair_first.gainPairTrials   = pChoice_gainPair_leftItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.gainPair_last.gainPairTrials    = pChoice_gainPair_leftItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.gainPair_first.gainPairTrials  = pChoice_gainPair_rightItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.gainPair_last.gainPairTrials   = pChoice_gainPair_rightItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.lossPair_first.lossPairTrials   = pChoice_lossPair_leftItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.lossPair_last.lossPairTrials    = pChoice_lossPair_leftItem_lossPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.lossPair_first.lossPairTrials  = pChoice_lossPair_rightItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.lossPair_last.lossPairTrials   = pChoice_lossPair_rightItem_lossPair_last;
            
            
            %% split by feedback type
            %% note: missing split by first/second half of trials for feedback split, you can add it if you want to perform a GLM with that
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.gainPair.gainPairTrials   =...
                pChoice_gainPair_gainItem_gainPair.*Q_gainPair_gainItem_gainPair + pChoice_gainPair_ntalItem_gainPair.*Q_gainPair_ntalItem_gainPair;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.lossPair.lossPairTrials   =...
                pChoice_lossPair_ntalItem_lossPair.*Q_lossPair_ntalItem_lossPair + pChoice_lossPair_lossItem_lossPair.*Q_lossPair_lossItem_lossPair;
            % first/second half
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.gainPair_first.gainPairTrials   =...
                pChoice_gainPair_gainItem_gainPair_first.*Q_gainPair_gainItem_gainPair_first + pChoice_gainPair_ntalItem_gainPair_first.*Q_gainPair_ntalItem_gainPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.gainPair_last.gainPairTrials   =...
                pChoice_gainPair_gainItem_gainPair_last.*Q_gainPair_gainItem_gainPair_last + pChoice_gainPair_ntalItem_gainPair_last.*Q_gainPair_ntalItem_gainPair_last;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.lossPair_first.lossPairTrials   =...
                pChoice_lossPair_ntalItem_lossPair_first.*Q_lossPair_ntalItem_lossPair_first + pChoice_lossPair_lossItem_lossPair_first.*Q_lossPair_lossItem_lossPair_first;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.lossPair_last.lossPairTrials   =...
                pChoice_lossPair_ntalItem_lossPair_last.*Q_lossPair_ntalItem_lossPair_last + pChoice_lossPair_lossItem_lossPair_last.*Q_lossPair_lossItem_lossPair_last;
            
            % classify according to feedback type
            Q_gainPair_gainItem_gainPair_gainFbk = Q_gainPair_gainItem(gainPair_gainFeedbackTrials);
            Q_gainPair_ntalItem_gainPair_gainFbk = Q_gainPair_gainItem(gainPair_gainFeedbackTrials);
            Q_gainPair_gainItem_gainPair_ntalFbk = Q_gainPair_ntalItem(gainPair_neutralFeedbackTrials);
            Q_gainPair_ntalItem_gainPair_ntalFbk = Q_gainPair_ntalItem(gainPair_neutralFeedbackTrials);
            Q_lossPair_ntalItem_lossPair_ntalFbk = Q_lossPair_ntalItem(lossPair_neutralFeedbackTrials);
            Q_lossPair_lossItem_lossPair_ntalFbk = Q_lossPair_ntalItem(lossPair_neutralFeedbackTrials);
            Q_lossPair_ntalItem_lossPair_lossFbk = Q_lossPair_lossItem(lossPair_lossFeedbackTrials);
            Q_lossPair_lossItem_lossPair_lossFbk = Q_lossPair_lossItem(lossPair_lossFeedbackTrials);
            % best/worse
            pChoice_gainPair_gainItem_gainPair_gainFbk = pChoice_bestItem(gainPair_gainFeedbackTrials);
            pChoice_gainPair_ntalItem_gainPair_gainFbk = pChoice_worseItem(gainPair_gainFeedbackTrials); % proba of choosing worse option = 1 - proba of choosing the best
            pChoice_gainPair_gainItem_gainPair_ntalFbk = pChoice_bestItem(gainPair_neutralFeedbackTrials);
            pChoice_gainPair_ntalItem_gainPair_ntalFbk = pChoice_worseItem(gainPair_neutralFeedbackTrials); % proba of choosing worse option = 1 - proba of choosing the best
            pChoice_lossPair_ntalItem_lossPair_ntalFbk = pChoice_bestItem(lossPair_neutralFeedbackTrials);
            pChoice_lossPair_lossItem_lossPair_ntalFbk = pChoice_worseItem(lossPair_neutralFeedbackTrials);
            pChoice_lossPair_ntalItem_lossPair_lossFbk = pChoice_bestItem(lossPair_lossFeedbackTrials);
            pChoice_lossPair_lossItem_lossPair_lossFbk = pChoice_worseItem(lossPair_lossFeedbackTrials);
            % left/right
            pChoice_gainPair_leftItem_gainPair_gainFbk = pChoice_leftItem(gainPair_gainFeedbackTrials);
            pChoice_gainPair_rightItem_gainPair_gainFbk = pChoice_rightItem(gainPair_gainFeedbackTrials); % proba of choosing worse option = 1 - proba of choosing the best
            pChoice_gainPair_leftItem_gainPair_ntalFbk = pChoice_leftItem(gainPair_neutralFeedbackTrials);
            pChoice_gainPair_rightItem_gainPair_ntalFbk = pChoice_rightItem(gainPair_neutralFeedbackTrials); % proba of choosing worse option = 1 - proba of choosing the best
            pChoice_lossPair_leftItem_lossPair_ntalFbk = pChoice_leftItem(lossPair_neutralFeedbackTrials);
            pChoice_lossPair_rightItem_lossPair_ntalFbk = pChoice_rightItem(lossPair_neutralFeedbackTrials);
            pChoice_lossPair_leftItem_lossPair_lossFbk = pChoice_leftItem(lossPair_lossFeedbackTrials);
            pChoice_lossPair_rightItem_lossPair_lossFbk = pChoice_rightItem(lossPair_lossFeedbackTrials);
            
            % Q.values by item type
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.GainTrial.GainItem                        = Q_gainPair_gainItem_gainPair_gainFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.GainTrial.NeutralItem                     = Q_gainPair_ntalItem_gainPair_gainFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.NeutralTrial.GainItem                     = Q_gainPair_gainItem_gainPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.NeutralTrial.NeutralItem                  = Q_gainPair_ntalItem_gainPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.NeutralTrial.NeutralItem                  = Q_lossPair_ntalItem_lossPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.NeutralTrial.LossItem                     = Q_lossPair_lossItem_lossPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.LossTrial.NeutralItem                     = Q_lossPair_ntalItem_lossPair_lossFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.LossTrial.LossItem                        = Q_lossPair_lossItem_lossPair_lossFbk;
            
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.GainTrial.bestItem                        = Q_gainPair_gainItem_gainPair_gainFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.GainTrial.worseItem                       = Q_gainPair_ntalItem_gainPair_gainFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.NeutralTrial.bestItem                     = Q_gainPair_gainItem_gainPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.NeutralTrial.worseItem                    = Q_gainPair_ntalItem_gainPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.NeutralTrial.bestItem                     = Q_lossPair_ntalItem_lossPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.NeutralTrial.worseItem                    = Q_lossPair_lossItem_lossPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.LossTrial.bestItem                        = Q_lossPair_ntalItem_lossPair_lossFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.LossTrial.worseItem                       = Q_lossPair_lossItem_lossPair_lossFbk;
            
            % Q.chosen
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.GainTrial.chosenItem                  = QchItem_all_wExclTrials(gainPair_gainFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.NeutralTrial.chosenItem               = QchItem_all_wExclTrials(gainPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.NeutralTrial.chosenItem               = QchItem_all_wExclTrials(lossPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.LossTrial.chosenItem                  = QchItem_all_wExclTrials(lossPair_lossFeedbackTrials);
            
            % Q.unchosen
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.GainTrial.unchosenItem                  = QunchItem_all_wExclTrials(gainPair_gainFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.gainPair.NeutralTrial.unchosenItem               = QunchItem_all_wExclTrials(gainPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.NeutralTrial.unchosenItem               = QunchItem_all_wExclTrials(lossPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.lossPair.LossTrial.unchosenItem                  = QunchItem_all_wExclTrials(lossPair_lossFeedbackTrials);
            
            % dQ (best-worse option)
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.gainPair.GainTrial                  = dQ_all(gainPair_gainFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.gainPair.NeutralTrial               = dQ_all(gainPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.lossPair.NeutralTrial               = dQ_all(lossPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ.lossPair.LossTrial                  = dQ_all(lossPair_lossFeedbackTrials);
            
            % dQ (left-right option)
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.gainPair.GainTrial                  = dQ_LR_wExclTrials(gainPair_gainFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.gainPair.NeutralTrial               = dQ_LR_wExclTrials(gainPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.lossPair.NeutralTrial               = dQ_LR_wExclTrials(lossPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.dQ_LR.lossPair.LossTrial                  = dQ_LR_wExclTrials(lossPair_lossFeedbackTrials);
            
            % PE
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.gainPair.GainTrial                  = PE_all(gainPair_gainFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.gainPair.NeutralTrial               = PE_all(gainPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.ntalPair.NeutralTrial               = PE_all(neutralPairsTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.lossPair.NeutralTrial               = PE_all(lossPair_neutralFeedbackTrials);
            learn.mod.(model_nm_learn).Q_model(iModel).raw.PE.lossPair.LossTrial                  = PE_all(lossPair_lossFeedbackTrials);
            
            % p(choice = best)
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.gainPair.GainTrial        = pChoice_gainPair_gainItem_gainPair_gainFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.gainPair.GainTrial       = pChoice_gainPair_ntalItem_gainPair_gainFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.gainPair.NeutralTrial     = pChoice_gainPair_gainItem_gainPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.gainPair.NeutralTrial    = pChoice_gainPair_ntalItem_gainPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.lossPair.NeutralTrial     = pChoice_lossPair_ntalItem_lossPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.lossPair.NeutralTrial    = pChoice_lossPair_lossItem_lossPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.best.lossPair.LossTrial        = pChoice_lossPair_ntalItem_lossPair_lossFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.worse.lossPair.LossTrial       = pChoice_lossPair_lossItem_lossPair_lossFbk;
            % p(choice = left/right)
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.gainPair.GainTrial        = pChoice_gainPair_leftItem_gainPair_gainFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.gainPair.GainTrial       = pChoice_gainPair_rightItem_gainPair_gainFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.gainPair.NeutralTrial     = pChoice_gainPair_leftItem_gainPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.gainPair.NeutralTrial    = pChoice_gainPair_rightItem_gainPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.lossPair.NeutralTrial     = pChoice_lossPair_leftItem_lossPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.lossPair.NeutralTrial    = pChoice_lossPair_rightItem_lossPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.left.lossPair.LossTrial        = pChoice_lossPair_leftItem_lossPair_lossFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.pChoice.right.lossPair.LossTrial       = pChoice_lossPair_rightItem_lossPair_lossFbk;
            
            % SV
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.gainPair.gainTrial   =...
                pChoice_gainPair_gainItem_gainPair_gainFbk.*Q_gainPair_gainItem_gainPair_gainFbk + pChoice_gainPair_ntalItem_gainPair_gainFbk.*Q_gainPair_ntalItem_gainPair_gainFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.gainPair.NeutralTrial   =...
                pChoice_gainPair_gainItem_gainPair_ntalFbk.*Q_gainPair_gainItem_gainPair_ntalFbk + pChoice_gainPair_ntalItem_gainPair_ntalFbk.*Q_gainPair_ntalItem_gainPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.lossPair.NeutralTrial   =...
                pChoice_lossPair_ntalItem_lossPair_ntalFbk.*Q_lossPair_ntalItem_lossPair_ntalFbk + pChoice_lossPair_lossItem_lossPair_ntalFbk.*Q_lossPair_lossItem_lossPair_ntalFbk;
            learn.mod.(model_nm_learn).Q_model(iModel).raw.SV.lossPair.lossTrial   =...
                pChoice_lossPair_ntalItem_lossPair_lossFbk.*Q_lossPair_ntalItem_lossPair_lossFbk + pChoice_lossPair_lossItem_lossPair_lossFbk.*Q_lossPair_lossItem_lossPair_lossFbk;
        end % model loop
    end % model type loop
    
    %% feedback (0/1 neutral/gain gainPair and loss/neutral for lossPair)
    % feedback (raw)
    best_fbk = feedback; %(-1/1 neutral/gain gainPair and loss/neutral for lossPair)
    best_fbk(feedback == -1) = 0; % transform -1 values into 0 to get a binary variable
    learn.mod.fbk.raw.allPairs.main                      = best_fbk(RL_nonExcludedTrials);
    learn.mod.fbk.raw.GL_Pairs.main                      = best_fbk(GL_PairsTrials);
    learn.mod.fbk.raw.GLPairs_first.main                = best_fbk(GLPairs_firstTrials);
    learn.mod.fbk.raw.GLPairs_last.main                 = best_fbk(GLPairs_lastTrials);
    learn.mod.fbk.raw.gainPair.main                      = best_fbk(gainPairsTrials);
    learn.mod.fbk.raw.ntalPair.main                      = best_fbk(neutralPairsTrials);
    learn.mod.fbk.raw.lossPair.main                      = best_fbk(lossPairsTrials);
    learn.mod.fbk.raw.gainPair.GainTrial                 = best_fbk(gainPair_gainFeedbackTrials);
    learn.mod.fbk.raw.gainPair.NeutralTrial              = best_fbk(gainPair_neutralFeedbackTrials);
    learn.mod.fbk.raw.ntalPair.NeutralTrial              = learn.mod.fbk.raw.ntalPair.main;
    learn.mod.fbk.raw.lossPair.NeutralTrial              = best_fbk(lossPair_neutralFeedbackTrials);
    learn.mod.fbk.raw.lossPair.LossTrial                 = best_fbk(lossPair_lossFeedbackTrials);
    % split first/second half trials
    learn.mod.fbk.raw.gainPair_first.main                = best_fbk(gainPair_firstTrials);
    learn.mod.fbk.raw.gainPair_last.main                 = best_fbk(gainPair_lastTrials);
    learn.mod.fbk.raw.ntalPair_first.main                = best_fbk(ntalPair_firstTrials);
    learn.mod.fbk.raw.ntalPair_last.main                 = best_fbk(ntalPair_lastTrials);
    learn.mod.fbk.raw.lossPair_first.main                = best_fbk(lossPair_firstTrials);
    learn.mod.fbk.raw.lossPair_last.main                 = best_fbk(lossPair_lastTrials);
    
    %% feedback bis (-1/0/1 loss/neutral/gain feedback)
    learn.mod.fbk_bis.raw.allPairs.main                      = gain(RL_nonExcludedTrials);
    learn.mod.fbk_bis.raw.GL_Pairs.main                      = gain(GL_PairsTrials);
    learn.mod.fbk_bis.raw.GLPairs_first.main                = gain(GLPairs_firstTrials);
    learn.mod.fbk_bis.raw.GLPairs_last.main                 = gain(GLPairs_lastTrials);
    learn.mod.fbk_bis.raw.gainPair.main                      = gain(gainPairsTrials);
    learn.mod.fbk_bis.raw.ntalPair.main                      = gain(neutralPairsTrials);
    learn.mod.fbk_bis.raw.lossPair.main                      = gain(lossPairsTrials);
    learn.mod.fbk_bis.raw.gainPair.GainTrial                 = gain(gainPair_gainFeedbackTrials);
    learn.mod.fbk_bis.raw.gainPair.NeutralTrial              = gain(gainPair_neutralFeedbackTrials);
    learn.mod.fbk_bis.raw.ntalPair.NeutralTrial              = learn.mod.fbk_bis.raw.ntalPair.main;
    learn.mod.fbk_bis.raw.lossPair.NeutralTrial              = gain(lossPair_neutralFeedbackTrials);
    learn.mod.fbk_bis.raw.lossPair.LossTrial                 = gain(lossPair_lossFeedbackTrials);
    % split first/second half trials
    learn.mod.fbk_bis.raw.gainPair_first.main                = gain(gainPair_firstTrials);
    learn.mod.fbk_bis.raw.gainPair_last.main                 = gain(gainPair_lastTrials);
    learn.mod.fbk_bis.raw.ntalPair_first.main                = gain(ntalPair_firstTrials);
    learn.mod.fbk_bis.raw.ntalPair_last.main                 = gain(ntalPair_lastTrials);
    learn.mod.fbk_bis.raw.lossPair_first.main                = gain(lossPair_firstTrials);
    learn.mod.fbk_bis.raw.lossPair_last.main                 = gain(lossPair_lastTrials);
    
    %% total gain
    learn.mod.totalGain.raw.allPairs.main                      = totalGain(RL_nonExcludedTrials);
    learn.mod.totalGain.raw.GL_Pairs.main                      = totalGain(GL_PairsTrials);
    learn.mod.totalGain.raw.GLPairs_first.main                  = totalGain(GLPairs_firstTrials);
    learn.mod.totalGain.raw.GLPairs_last.main                   = totalGain(GLPairs_lastTrials);
    learn.mod.totalGain.raw.gainPair.main                      = totalGain(gainPairsTrials);
    learn.mod.totalGain.raw.ntalPair.main                      = totalGain(neutralPairsTrials);
    learn.mod.totalGain.raw.lossPair.main                      = totalGain(lossPairsTrials);
    learn.mod.totalGain.raw.gainPair.GainTrial                 = totalGain(gainPair_gainFeedbackTrials);
    learn.mod.totalGain.raw.gainPair.NeutralTrial              = totalGain(gainPair_neutralFeedbackTrials);
    learn.mod.totalGain.raw.ntalPair.NeutralTrial              = learn.mod.totalGain.raw.ntalPair.main;
    learn.mod.totalGain.raw.lossPair.NeutralTrial              = totalGain(lossPair_neutralFeedbackTrials);
    learn.mod.totalGain.raw.lossPair.LossTrial                 = totalGain(lossPair_lossFeedbackTrials);
    % split first/second half trials
    learn.mod.totalGain.raw.gainPair_first.main                = totalGain(gainPair_firstTrials);
    learn.mod.totalGain.raw.gainPair_last.main                 = totalGain(gainPair_lastTrials);
    learn.mod.totalGain.raw.ntalPair_first.main                = totalGain(ntalPair_firstTrials);
    learn.mod.totalGain.raw.ntalPair_last.main                 = totalGain(ntalPair_lastTrials);
    learn.mod.totalGain.raw.lossPair_first.main                = totalGain(lossPair_firstTrials);
    learn.mod.totalGain.raw.lossPair_last.main                 = totalGain(lossPair_lastTrials);
    
    %% side best option
    learn.mod.sideBest.raw.allPairs.main            = sideBest(RL_nonExcludedTrials);
    learn.mod.sideBest.raw.GL_Pairs.main            = sideBest(GL_PairsTrials);
    learn.mod.sideBest.raw.GLPairs_first.main       = sideBest(GLPairs_firstTrials);
    learn.mod.sideBest.raw.GLPairs_last.main        = sideBest(GLPairs_lastTrials);
    learn.mod.sideBest.raw.gainPair.main            = sideBest(gainPairsTrials);
    learn.mod.sideBest.raw.ntalPair.main            = sideBest(neutralPairsTrials);
    learn.mod.sideBest.raw.lossPair.main            = sideBest(lossPairsTrials);
    learn.mod.sideBest.raw.gainPair.GainTrial       = sideBest(gainPair_gainFeedbackTrials);
    learn.mod.sideBest.raw.gainPair.NeutralTrial    = sideBest(gainPair_neutralFeedbackTrials);
    learn.mod.sideBest.raw.ntalPair.NeutralTrial    = learn.mod.sideBest.raw.ntalPair.main;
    learn.mod.sideBest.raw.lossPair.NeutralTrial    = sideBest(lossPair_neutralFeedbackTrials);
    learn.mod.sideBest.raw.lossPair.LossTrial       = sideBest(lossPair_lossFeedbackTrials);
    % split first/second half trials
    learn.mod.sideBest.raw.gainPair_first.main      = sideBest(gainPair_firstTrials);
    learn.mod.sideBest.raw.gainPair_last.main       = sideBest(gainPair_lastTrials);
    learn.mod.sideBest.raw.ntalPair_first.main      = sideBest(ntalPair_firstTrials);
    learn.mod.sideBest.raw.ntalPair_last.main       = sideBest(ntalPair_lastTrials);
    learn.mod.sideBest.raw.lossPair_first.main      = sideBest(lossPair_firstTrials);
    learn.mod.sideBest.raw.lossPair_last.main       = sideBest(lossPair_lastTrials);
    
    %% choice best option
    learn.mod.choiceBest.raw.allPairs.main          = choice_BestOrWorse(RL_nonExcludedTrials);
    learn.mod.choiceBest.raw.GL_Pairs.main          = choice_BestOrWorse(GL_PairsTrials);
    learn.mod.choiceBest.raw.GLPairs_first.main     = choice_BestOrWorse(GLPairs_firstTrials);
    learn.mod.choiceBest.raw.GLPairs_last.main      = choice_BestOrWorse(GLPairs_lastTrials);
    learn.mod.choiceBest.raw.gainPair.main          = choice_BestOrWorse(gainPairsTrials);
    learn.mod.choiceBest.raw.ntalPair.main          = choice_BestOrWorse(neutralPairsTrials);
    learn.mod.choiceBest.raw.lossPair.main          = choice_BestOrWorse(lossPairsTrials);
    learn.mod.choiceBest.raw.gainPair.GainTrial     = choice_BestOrWorse(gainPair_gainFeedbackTrials);
    learn.mod.choiceBest.raw.gainPair.NeutralTrial  = choice_BestOrWorse(gainPair_neutralFeedbackTrials);
    learn.mod.choiceBest.raw.ntalPair.NeutralTrial  = learn.mod.choiceBest.raw.ntalPair.main;
    learn.mod.choiceBest.raw.lossPair.NeutralTrial  = choice_BestOrWorse(lossPair_neutralFeedbackTrials);
    learn.mod.choiceBest.raw.lossPair.LossTrial     = choice_BestOrWorse(lossPair_lossFeedbackTrials);
    % split first/second half trials
    learn.mod.choiceBest.raw.gainPair_first.main    = choice_BestOrWorse(gainPair_firstTrials);
    learn.mod.choiceBest.raw.gainPair_last.main     = choice_BestOrWorse(gainPair_lastTrials);
    learn.mod.choiceBest.raw.ntalPair_first.main    = choice_BestOrWorse(ntalPair_firstTrials);
    learn.mod.choiceBest.raw.ntalPair_last.main     = choice_BestOrWorse(ntalPair_lastTrials);
    learn.mod.choiceBest.raw.lossPair_first.main    = choice_BestOrWorse(lossPair_firstTrials);
    learn.mod.choiceBest.raw.lossPair_last.main     = choice_BestOrWorse(lossPair_lastTrials);
    
    %% valence pair (-1) loss; (0) neutral; (1) gain
    valence = npair;
    valence(npair == 3) = -1;
    valence(npair == 2) = 0;
    learn.mod.valence.raw.allPairs.main         = valence(RL_nonExcludedTrials);
    learn.mod.valence.raw.GL_Pairs.main         = valence(GL_PairsTrials);
    learn.mod.valence.raw.GLPairs_first.main    = valence(GLPairs_firstTrials);
    learn.mod.valence.raw.GLPairs_last.main     = valence(GLPairs_lastTrials);
    learn.mod.valence.raw.gainPair.main         = valence(gainPairsTrials);
    learn.mod.valence.raw.ntalPair.main         = valence(neutralPairsTrials);
    learn.mod.valence.raw.lossPair.main         = valence(lossPairsTrials);
    learn.mod.valence.raw.gainPair.GainTrial    = valence(gainPair_gainFeedbackTrials);
    learn.mod.valence.raw.gainPair.NeutralTrial = valence(gainPair_neutralFeedbackTrials);
    learn.mod.valence.raw.ntalPair.NeutralTrial = learn.mod.valence.raw.ntalPair.main;
    learn.mod.valence.raw.lossPair.NeutralTrial = valence(lossPair_neutralFeedbackTrials);
    learn.mod.valence.raw.lossPair.LossTrial    = valence(lossPair_lossFeedbackTrials);
    % split first/second half trials
    learn.mod.valence.raw.gainPair_first.main   = valence(gainPair_firstTrials);
    learn.mod.valence.raw.gainPair_last.main    = valence(gainPair_lastTrials);
    learn.mod.valence.raw.ntalPair_first.main   = valence(ntalPair_firstTrials);
    learn.mod.valence.raw.ntalPair_last.main    = valence(ntalPair_lastTrials);
    learn.mod.valence.raw.lossPair_first.main   = valence(lossPair_firstTrials);
    learn.mod.valence.raw.lossPair_last.main    = valence(lossPair_lastTrials);
    
    %% pool all data which is relevant for the GLM together in a new file
    save([fMRI_Dir 'onsets_sub',subid,'_learning_run',sessnber,'.mat'],'learn');
    
    %% clear learn structure before next run to avoid interferences
    clear('learn');
end % runs

end % function