function[badSubs, allSubs]=find_bad_Ech_runs(study_nm, condition, subject_id, NS)
% [badSubs, allSubs]=find_bad_Ech_runs(study_nm, condition, subject_id, NS)
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
%
% allSubs: structure with info about all subjects and runs

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
n_hE_lvl = 3;
n_runsPerTask = 2;
[choice_hE_perhE, choice_lE_perhE,...
    choice_hE_perhE_yn, choice_lE_perhE_yn,] = deal(NaN(n_hE_lvl, n_runsPerTask, NS));
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
        %% load the data
        behaviorStruct_tmp = getfield(load([subBehaviorFolder,...
            'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
            '_task_behavioral_tmp.mat']),'summary');
        % extract relevant info
        hE_lvl_tmp = max(behaviorStruct_tmp.choiceOptions.E.left, behaviorStruct_tmp.choiceOptions.E.right);
        Ech_tmp = behaviorStruct_tmp.E_chosen;
        for iE = 1:n_hE_lvl
            curr_E_idx = hE_lvl_tmp == iE;
            option_chosen_tmp = Ech_tmp(curr_E_idx);
            % extract number of times that choice was high or low effort
            % for the current effort level
            choice_hE_perhE(iE, taskRun_idx, iS) = sum(option_chosen_tmp == iE);
            choice_lE_perhE(iE, taskRun_idx, iS) = sum(option_chosen_tmp == 0);
            
            % extract binary variable indicating whether at least one
            % choice of the high or low effort was made for the current
            % effort level
            if sum(option_chosen_tmp == iE) > 0
                choice_hE_perhE_yn(iE, taskRun_idx, iS) = 1;
            else
                choice_hE_perhE_yn(iE, taskRun_idx, iS) = 0;
            end
            % same but for low effort option selection
            if sum(option_chosen_tmp == 0) > 0
                choice_lE_perhE_yn(iE, taskRun_idx, iS) = 1;
            else
                choice_lE_perhE_yn(iE, taskRun_idx, iS) = 0;
            end
        end % effort level loop
        
        % store information
        allSubs.(['CID',sub_nm]).(['run',run_nm]).choice_highE = choice_hE_perhE(:,taskRun_idx, iS);
        allSubs.(['CID',sub_nm]).(['run',run_nm]).choice_lowE = choice_lE_perhE(:,taskRun_idx, iS);
        % bad subjects = those who have not at least 2 elements for the
        % curve of choice = f(E level) in both high effort and low effort
        % choice
        if sum(choice_hE_perhE_yn(:,taskRun_idx, iS)) < 2 || sum(choice_lE_perhE_yn(:,taskRun_idx, iS)) < 2
            badSubs.(['CID',sub_nm]).(['run',run_nm]).choice_highE = choice_hE_perhE(:,taskRun_idx, iS);
            badSubs.(['CID',sub_nm]).(['run',run_nm]).choice_lowE = choice_lE_perhE(:,taskRun_idx, iS);
        end
        % bad subjects bis = those who have not some element for all the
        % points of the curve of choice = f(E level) in both high effort 
        % and low effort choice
        if sum(choice_hE_perhE_yn(:,taskRun_idx, iS)) < 3 || sum(choice_lE_perhE_yn(:,taskRun_idx, iS)) < 3
            badSubs.bis.(['CID',sub_nm]).(['run',run_nm]).choice_highE = choice_hE_perhE(:,taskRun_idx, iS);
            badSubs.bis.(['CID',sub_nm]).(['run',run_nm]).choice_lowE = choice_lE_perhE(:,taskRun_idx, iS);
        end
    end % run loop
end % subject loop

end % function