function [avg_RT] = avg_RT_perSub(study_nm, sub_nm, condition)
% [avg_RT] = avg_RT_perSub(study_nm, sub_nm, condition)
% avg_RT_perSub will extract average reaction times (RT) for the subject 
% defined in input. The average is extracted across all tasks and
% separately for each task as well.
%
% INPUTS
% study_nm: study name
%
% sub_nm: subject name
%
% condition: will define which runs should be considered
%
% OUTPUTS
% avg_RT: structure with average reaction times (RT)

%% working directories
% computerRoot = LGCM_root_paths;
computerRoot = ['E:',filesep];
dataRoot = [computerRoot, filesep, study_nm, filesep];

%% general parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
n_runsPerTask = 2;
for iT = 1:nTasks
    task_nm1 = tasks{iT};
    [RT.(task_nm1).perRun] = deal(NaN(n_runsPerTask,1));
end % task loop

subBehaviorFolder = [dataRoot, filesep, 'CID',sub_nm, filesep, 'behavior',filesep];

%% loop through runs
[runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
okRuns = runs.runsToKeep;
for iRun = 1:n_runs
    kRun = okRuns(iRun);
    run_nm = num2str(kRun);
    switch kRun
        case {1,2}
            jRun = 1;
        case {3,4}
            jRun = 2;
    end
    task_short_nm = runs.tasks{iRun};
    task_fullName = task_fullName_extraction(task_short_nm);
    [RT_tmp] = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
    RT.(task_short_nm).perRun(jRun) = mean(RT_tmp,2,'omitnan');
end % run loop

%% average data
avg_RT.Ep = mean(RT.Ep.perRun, 1, 'omitnan');
avg_RT.Em = mean(RT.Em.perRun, 1, 'omitnan');
avg_RT.EpEm = mean([avg_RT.Ep, avg_RT.Em],2,'omitnan');

end % function