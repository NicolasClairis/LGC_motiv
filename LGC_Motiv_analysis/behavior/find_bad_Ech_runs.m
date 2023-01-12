function[badSubs]=find_bad_Ech_runs(study_nm, condition, subject_id, NS)
% [badSubs]=find_bad_Ech_runs(study_nm, condition, subject_id, NS)
% The goal of find_bad_Ech_runs is to identify the runs and subjects who
% have not enough power in order to compare high E chosen vs low E chosen
% slopes with high effort level as a reference.
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% condition: condition (will be asked through subject_condition.m if not
% filled)
%
% subject_id: list of subjects to look at
%
% NS: number of subjects
%
% OUTPUTS
% badSubs: structure with info about the bad subjects and which runs should
% be removed for this analysis

%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% subject selection
if ~exist('subject_id','var') || ~exist('NS','var') ||...
        isempty(subject_id) || isempty(NS)
    if ~exist('condition','var') || isempty(condition)
        condition = subject_condition;
    end
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
end
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
        % define trial index for relevant variable to extract
        switch behavioral_task_to_look
            case {'Ep','Em'}
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(taskRun_idx-1);
            case 'EpEmPool'
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(kRun-1);
        end
        %% load the data
        behaviorStruct_tmp = load([subBehaviorFolder,...
            'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
            '_task.mat']);
            
    end % run loop
end % subject loop