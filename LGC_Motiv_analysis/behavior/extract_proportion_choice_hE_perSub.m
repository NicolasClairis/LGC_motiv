function[choice_hE] = extract_proportion_choice_hE_perSub(study_nm, subject_id, condition, task_nm)
% [choice_hE] = extract_proportion_choice_hE_perSub(study_nm, subject_id, condition, task_nm)
% extract_proportion_choice_hE_perSub will do a loop across subjects and
% will extract the proportion of high effort choices per subject across all
% trials and per effort level.
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% subject_id: list of subjects to look
%
% condition: which runs to consider
%
% task_nm: task name ('Ep'/'Em') of the task to consider
%
% OUTPUTS
% choice_hE: structure with relevant information for proportion of high
% effort choices

%% working directory
computerRoot = 'E:';
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% main parameters
[task_fullName] = task_fullName_extraction(task_nm);
nRunsPerTask = 2;
n_hE_levels = 3;
NS = length(subject_id);

%% initialize variables of interest
[choice_hE.(task_nm).perRun_allSubs] = deal(NaN(n_hE_levels, NS, nRunsPerTask));
[choice_hE.(task_nm).allSubs] = deal(NaN(n_hE_levels, NS));
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
    runsStruct = runs_definition(study_nm, sub_nm, condition);
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    for iRun = 1:length(okRuns)
        kRun = okRuns(iRun);
        run_nm = num2str(kRun);
        switch kRun
            case {1,2}
                taskRun_idx = 1;
            case {3,4}
                taskRun_idx = 2;
        end
        task_nm_tmp = taskNames{iRun};
        if strcmp(task_nm_tmp, task_nm)
            choice_highE_tmp = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            hE_level_tmp = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            
            for iE = 1:n_hE_levels
                E_idx = hE_level_tmp == iE;
                choice_highE_perE_tmp = choice_highE_tmp(E_idx);
                choice_proportion_tmp = sum(choice_highE_perE_tmp == 1)./sum(ismember(choice_highE_perE_tmp,[0,1]));
                choice_hE.(task_nm).perRun_allSubs(iE,iS,taskRun_idx) = choice_proportion_tmp;
            end % effort level
        end % task filter
    end % run loop
    
    %% average across runs
    for iE = 1:n_hE_levels
        choice_hE.(task_nm).allSubs(iE,iS) = mean(choice_hE.(task_nm).perRun_allSubs(iE,iS,:),3,'omitnan');
    end % effort level
end % subject loop

end % function