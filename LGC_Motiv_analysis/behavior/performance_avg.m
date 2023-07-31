function[perf_avg] = performance_avg(study_nm, subject_id, condition)
% [perf_avg] = performance_avg(study_nm, subject_id, condition)
% performance_avg will look at the average performance for each task of
% each subject
%
% INPUTS
% study_nm: study name
%
% subject_id: list of subjects
%
% condition: condition
%
% OUTPUTS
% perf_avg: structure with average performance

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
end

%% working directory
pcRoot = LGCM_root_paths;
study_path = [pcRoot, study_nm, filesep];

%% initialize variables of interest
task_names = {'Ep','Em'};
nTasks = length(task_names);
n_maxRunsPerTask = 2;
[perf_avg.Ep.avg_perSub_perRun,...
    perf_avg.Em.avg_perSub_perRun] = deal(NaN(n_maxRunsPerTask, NS));
[perf_avg.Ep.avg_perSub,...
    perf_avg.Em.avg_perSub] = deal(NaN(1, NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [study_path, 'CID',sub_nm, filesep, 'behavior', filesep];

    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    for iRun = 1:n_runs
        kRun = runs.runsToKeep(iRun);
        run_nm = num2str(kRun);
        task_run_id = runs.tasks{iRun};
        switch task_run_id
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end
        switch kRun
            case {1,2}
                jRun = 1;
            case {3,4}
                jRun = 2;
        end
        [perf_run_tmp] = extract_perf(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        perf_avg.(task_run_id).avg_perSub_perRun(jRun, iS) = mean(perf_run_tmp,'omitnan');
    end % run loop

    for iTask = 1:nTasks
        task_nm = task_names{iTask};
        perf_avg.(task_nm).avg_perSub(iS) = mean(perf_avg.(task_nm).avg_perSub_perRun(:,iS),1,'omitnan');
    end
end % subject loop

for iTask = 1:nTasks
    task_nm = task_names{iTask};
    [perf_avg.(task_nm).avg,...
        perf_avg.(task_nm).sem] = mean_sem_sd(perf_avg.(task_nm).avg_perSub,2);
end
end % function