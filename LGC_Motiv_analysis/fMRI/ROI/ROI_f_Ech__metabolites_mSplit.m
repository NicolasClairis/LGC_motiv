%% extract the level of activity of a given ROI in function of each effort level chosen
% and split subjects with a median split based on the metabolites

%% general parameters
% display figures?
dispFig = true;

%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directory
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% subject definition
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load ROI
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);
% define which ROI, and which time period is of interest to you
% define ROI
ROI_names = fieldnames(ROI_trial_b_trial);
ROI_subList = ROI_trial_b_trial.subject_id;
ROI_names(strcmp(ROI_names,'subject_id')) = [];
if length(ROI_names) > 1
    error(['There should be only 1 ROI selected, not ',num2str(length(ROI_names))])
else
    ROI_nm = ROI_names;
end
ROI_short_nm = inputdlg('ROI short name?');
ROI_short_nm = ROI_short_nm{1};
% define task
task_names = {'Ep','Em','EpEmPool'};
which_task = listdlg('PromptString','Which task?','ListString',task_names);
task_to_look = task_names{which_task};
% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};

%% load behavior 
nRuns = 4;
nTrialsPerRun = 54;
nTrials = nTrialsPerRun*nRuns;
% input/mediator/output
[input_prm, ROI_mediator, choice_hE] = deal(NaN(nTrials, NS));

%% prepare median split based on metabolites
[low_met_subs, high_met_subs, metabolite_nm, ROI_nm] = medSplit_metabolites(study_nm, subject_id);

%% perform the median split

%% figure