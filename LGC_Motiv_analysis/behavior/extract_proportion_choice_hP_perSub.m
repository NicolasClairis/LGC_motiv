function[choice_hEhP] = extract_proportion_choice_hP_perSub(study_nm, subject_id, condition)
% [choice_hEhP] = extract_proportion_choice_hP_perSub(study_nm, subject_id, condition)
% extract_proportion_choice_hP_perSub will do a loop across subjects and
% will extract the proportion of high effort/high punishment choices per subject across all
% trials and per punishment level + split by task
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% subject_id: list of subjects to look
%
% condition: which runs to consider
%
% OUTPUTS
% choice_hEhP: structure with relevant information for proportion of high
% effort/high punishment choices

%% working directory
computerRoot = 'E:';
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% main parameters
nRunsPerTask = 2;
nRuns = 4;
n_hP_levels = 3;
NS = length(subject_id);

%% initialize variables of interest
[choice_hEhP.Ep.perRun_allSubs,...
    choice_hEhP.Ep.perRun_allSubs] = deal(NaN(n_hP_levels, NS, nRunsPerTask));
choice_hEhP.all.perRun_allSubs = NaN(n_hP_levels, NS, nRuns);
[choice_hEhP.Ep.allSubs,...
    choice_hEhP.Em.allSubs,...
    choice_hEhP.all.allSubs] = deal(NaN(n_hP_levels, NS));
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
    runsStruct = runs_definition(study_nm, sub_nm, condition);
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    for iRun = 1:length(okRuns)
        jRun = okRuns(iRun);
        run_nm = num2str(jRun);
        switch jRun
            case {1,2}
                kRun = 1;
            case {3,4}
                kRun = 2;
        end
        task_nm_tmp = taskNames{iRun};
        [task_fullName] = task_fullName_extraction(task_nm_tmp); % convert "Ep/Em" into "physical/mental"
        choice_highE_tmp = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        hP_level_tmp = extract_hP_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);

        for iP = 1:n_hP_levels
            P_idx = hP_level_tmp == iP;
            choice_hE_perP_tmp = choice_highE_tmp(P_idx);
            choice_proportion_tmp = sum(choice_hE_perP_tmp == 1)./sum(ismember(choice_hE_perP_tmp,[0,1]));
            choice_hEhP.(task_nm_tmp).perRun_allSubs(iP,iS,kRun) = choice_proportion_tmp;
            choice_hEhP.all.perRun_allSubs(iP,iS,jRun) = choice_proportion_tmp;
        end % effort level
    end % run loop
    
    %% average per subject across runs
    for iP = 1:n_hP_levels
        choice_hEhP.Ep.allSubs(iP,iS) = mean(choice_hEhP.Ep.perRun_allSubs(iP,iS,:),3,'omitnan');
        choice_hEhP.Em.allSubs(iP,iS) = mean(choice_hEhP.Em.perRun_allSubs(iP,iS,:),3,'omitnan');
        choice_hEhP.all.allSubs(iP,iS) = mean(choice_hEhP.all.perRun_allSubs(iP,iS,:),3,'omitnan');
    end % effort level
end % subject loop

end % function