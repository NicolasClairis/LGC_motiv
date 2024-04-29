function[accuracy, n_errors] = extract_Em_accuracy_and_errors_aRuns(study_nm, subject_id, condition)
% [accuracy, n_errors] = extract_Em_accuracy_and_errors_aRuns(study_nm, subject_id, condition)
% extract_Em_accuracy_and_errors_aRuns will extract the average accuracy
% (defined as (n.correct - n.errors)/(n.correct + n.errors)) and the
% average number of errors performed by each subject across the 2 mental
% effort runs.
%
% INPUTS (most will be asked or defined automatically if not defined in
% inputs)
%
% study_nm: study name
%
% subject_id: list of subjects
%
% condition: condition
%
% OUTPUTS
% accuracy: 1*NS vector with average accuracy
%
% n.errors: 1*NS vector with average number of errors

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end

%% working directory
% pcRoot = LGCM_root_paths;
pcRoot = 'E:\';
study_path = [pcRoot, study_nm, filesep];

%% initialize variable of interest
n_maxRunsPerTask = 2;
[accuracy_perR,...
    n_errors_perR] = deal(NaN(n_maxRunsPerTask, NS));
[accuracy, n_errors] = deal(NaN(1,NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [study_path, 'CID',sub_nm, filesep, 'behavior', filesep];
    
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    for iRun = 1:n_runs
        kRun = runs.runsToKeep(iRun);
        run_nm = num2str(kRun);
        task_run_id = runs.tasks{iRun};
        if strcmp(task_run_id,'Em')
            switch kRun
                case {1,2}
                    jRun = 1;
                case {3,4}
                    jRun = 2;
            end
            % extract number of correct and number of errors for current
            % run
            [~, n_correct_structTmp, n_errors_structTmp] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
            n_correct_tmp = n_correct_structTmp.allTrials;
            n_errors_tmp = n_errors_structTmp.allTrials;
            n_accuracy_tmp = (n_correct_tmp - n_errors_tmp)./(n_correct_tmp + n_errors_tmp);
            accuracy_perR(jRun, iS) = mean(n_accuracy_tmp,2,'omitnan');
            n_errors_perR(jRun, iS) = mean(n_errors_tmp,2,'omitnan');
        end
    end % run loop
    
    % average across runs
    accuracy(iS) = mean(accuracy_perR(:,iS),1,'omitnan');
    n_errors(iS) = mean(n_errors_perR(:,iS),1,'omitnan');
end % subject loop

end % function