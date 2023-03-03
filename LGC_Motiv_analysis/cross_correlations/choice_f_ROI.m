function[] = choice_f_ROI()
% choice_f_ROI will look at percentage of choices depending on ROI level of
% activity. Then it will also split the data depending on the effort level
% proposed and try to look whether the results vary accordingly.


%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% subject selection
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% general parameters
if ~exist('nBins','var') || isempty(nBins)
    nBins = 9;
end
tasks = {'Ep','Em'};
nTasks = length(tasks);
nTrialsPerRun = 54;
nRunsPerTask = 2;
nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
[fMRI_allTrials.Ep, choice_hE_allTrials.Ep,...
    fMRI_allTrials.Em, choice_hE_allTrials.Em] = deal(NaN(nTrialsPerTask, NS));
[fMRI_bins.Ep, choice_hE_bins.Ep,...
    fMRI_bins.Em, choice_hE_bins.Em] = deal(NaN(nBins, NS));

%% extract ROI activity for all subjects
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);
% define which ROI, and which time period is of interest to you
% define ROI
[fMRI_ROI_nm, fMRI_ROI_short_nm,...
    ROI_task_to_look,...
    timePeriod_nm] = extract_ROI_betas_onsets_only_questInfos(ROI_trial_b_trial);

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
    % extract runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    jRun = 0;
    for iRun = 1:length(okRuns)
        kRun = okRuns(iRun);
        run_nm = num2str(kRun);
        jRun = jRun + 1;
        task_nm_tmp = taskNames{jRun};
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        % define which task session it is
        switch kRun
            case {1,2}
                taskRun_idx = 1;
            case {3,4}
                taskRun_idx = 2;
        end
        run_nm_bis = ['run',num2str(taskRun_idx)];
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(taskRun_idx-1);
        
        % extract fMRI ROI as well
    end % run loop
    
    %% extract the bins
    
end % subject loop

%% average data across subjects

%% figure

end % function