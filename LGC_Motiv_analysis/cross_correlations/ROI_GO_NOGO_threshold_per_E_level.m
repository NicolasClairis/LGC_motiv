
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

% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};

%% main parameters of interest
task_names = {'Ep','Em'};
nTasks = length(task_names);
nTrialsPerRun = 54;
nRunsPerTask = 2;
n_hE_lvl = 3;
for iT = 1:nTasks
    task_nm = task_names{iT};
    % average per run per subject
    [ROI_meanDelta_yn_perSubPerRun.(task_nm).E1,...
            ROI_meanDelta_yn_perSubPerRun.(task_nm).E2,...
            ROI_meanDelta_yn_perSubPerRun.(task_nm).E3,...
            ROI_min_choice_yes_perSubPerRun.(task_nm).E1,...
            ROI_min_choice_yes_perSubPerRun.(task_nm).E2,...
            ROI_min_choice_yes_perSubPerRun.(task_nm).E3,...
            ROI_max_choice_yes_perSubPerRun.(task_nm).E1,...
            ROI_max_choice_yes_perSubPerRun.(task_nm).E2,...
            ROI_max_choice_yes_perSubPerRun.(task_nm).E3,...
            ROI_min_choice_no_perSubPerRun.(task_nm).E1,...
            ROI_min_choice_no_perSubPerRun.(task_nm).E2,...
            ROI_min_choice_no_perSubPerRun.(task_nm).E3,...
            ROI_max_choice_no_perSubPerRun.(task_nm).E1,...
            ROI_max_choice_no_perSubPerRun.(task_nm).E2,...
            ROI_max_choice_no_perSubPerRun.(task_nm).E3] = deal(NaN(nRunsPerTask,NS));
        
    % average per subject
    [ROI_meanDelta_yn.(task_nm).E1,...
        ROI_meanDelta_yn.(task_nm).E2,...
        ROI_meanDelta_yn.(task_nm).E3,...
        ROI_min_choice_yes.(task_nm).E1,...
        ROI_min_choice_yes.(task_nm).E2,...
        ROI_min_choice_yes.(task_nm).E3,...
        ROI_max_choice_yes.(task_nm).E1,...
        ROI_max_choice_yes.(task_nm).E2,...
        ROI_max_choice_yes.(task_nm).E3,...
        ROI_min_choice_no.(task_nm).E1,...
        ROI_min_choice_no.(task_nm).E2,...
        ROI_min_choice_no.(task_nm).E3,...
        ROI_max_choice_no.(task_nm).E1,...
        ROI_max_choice_no.(task_nm).E2,...
        ROI_max_choice_no.(task_nm).E3] = deal(NaN(1,NS));
end

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    if ~ismember(sub_nm,{'054','055','061'})
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
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(taskRun_idx-1);
        
        %% check if run is ok
        [choices] = count_choices_perElvl(subBehaviorFolder,...
            sub_nm, run_nm, task_fullName);
        if choices.allTrialsFullFilled == 1
            %% extract index of each effort level
            hE_level = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            %% extract choice for each trial (high E or low E)
            choice_hE = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            
            %% extract relevant information
            for iE = 1:n_hE_lvl
                E_nm = ['E',num2str(iE)];
                hE_idx = hE_level == iE;
                choice_yes = choice_hE == 1;
                choice_no = choice_hE == 0;
                hE_choice_yes_idx = hE_idx.*choice_yes == 1;
                hE_choice_no_idx = hE_idx.*choice_no == 1;
                
                % ROI BOLD for choice yes with current effort level
                ROI_hE_choice_yes_tmp = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(['run',num2str(taskRun_idx)]).(timePeriod_nm)(hE_choice_yes_idx,iS);
                ROI_hE_choice_no_tmp = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(['run',num2str(taskRun_idx)]).(timePeriod_nm)(hE_choice_no_idx,iS);
                ROI_meanDelta_yn_perSubPerRun.(task_nm_tmp).(E_nm)(taskRun_idx,iS) = mean(ROI_hE_choice_yes_tmp,1,'omitnan') - mean(ROI_hE_choice_no_tmp,1,'omitnan');
                ROI_min_choice_yes_perSubPerRun.(task_nm_tmp).(E_nm)(taskRun_idx,iS) = min(ROI_hE_choice_yes_tmp,[],1,'omitnan');
                ROI_max_choice_yes_perSubPerRun.(task_nm_tmp).(E_nm)(taskRun_idx,iS) = max(ROI_hE_choice_yes_tmp,[],1,'omitnan');
                ROI_min_choice_no_perSubPerRun.(task_nm_tmp).(E_nm)(taskRun_idx,iS) = min(ROI_hE_choice_no_tmp,[],1,'omitnan');
                ROI_max_choice_no_perSubPerRun.(task_nm_tmp).(E_nm)(taskRun_idx,iS) = max(ROI_hE_choice_no_tmp,[],1,'omitnan');
            end % effort level loop
        end % filter choices
    end % run loop
    
    % average per task across runs
    for iT = 1:nTasks
        task_nm = task_names{iT};
        for iE = 1:n_hE_lvl
            E_nm = ['E',num2str(iE)];
            ROI_meanDelta_yn.(task_nm).(E_nm)(iS) = mean(ROI_meanDelta_yn_perSubPerRun.(task_nm).(E_nm)(:,iS),1,'omitnan');
            ROI_min_choice_yes.(task_nm).(E_nm)(iS) = mean(ROI_min_choice_yes_perSubPerRun.(task_nm).(E_nm)(:,iS),1,'omitnan');
            ROI_max_choice_yes.(task_nm).(E_nm)(iS) = mean(ROI_max_choice_yes_perSubPerRun.(task_nm).(E_nm)(:,iS),1,'omitnan');
            ROI_min_choice_no.(task_nm).(E_nm)(iS) = mean(ROI_min_choice_no_perSubPerRun.(task_nm).(E_nm)(:,iS),1,'omitnan');
            ROI_max_choice_no.(task_nm).(E_nm)(iS) = mean(ROI_max_choice_no_perSubPerRun.(task_nm).(E_nm)(:,iS),1,'omitnan');
        end % effort loop
    end % task loop
    end
end % subject loop

%% average across individuals
for iT = 1:nTasks
    task_nm = task_names{iT};
    % mean delta
    m_ROI_meanDelta_yn_allSubs.(task_nm) = [mean(ROI_meanDelta_yn.(task_nm).E1,2,'omitnan'),...
        mean(ROI_meanDelta_yn.(task_nm).E2,2,'omitnan'),...
        mean(ROI_meanDelta_yn.(task_nm).E3,2,'omitnan')];
    sem_ROI_meanDelta_yn_allSubs.(task_nm) = [sem(ROI_meanDelta_yn.(task_nm).E1,2),...
        sem(ROI_meanDelta_yn.(task_nm).E2,2),...
        sem(ROI_meanDelta_yn.(task_nm).E3,2)];
    % min yes
    m_ROI_min_choice_yes_allSubs.(task_nm) = [mean(ROI_min_choice_yes.(task_nm).E1,2,'omitnan'),...
        mean(ROI_min_choice_yes.(task_nm).E2,2,'omitnan'),...
        mean(ROI_min_choice_yes.(task_nm).E3,2,'omitnan')];
    sem_ROI_min_choice_yes_allSubs.(task_nm) = [sem(ROI_min_choice_yes.(task_nm).E1,2),...
        sem(ROI_min_choice_yes.(task_nm).E2,2),...
        sem(ROI_min_choice_yes.(task_nm).E3,2)];
    % min no
    m_ROI_min_choice_no_allSubs.(task_nm) = [mean(ROI_min_choice_no.(task_nm).E1,2,'omitnan'),...
        mean(ROI_min_choice_no.(task_nm).E2,2,'omitnan'),...
        mean(ROI_min_choice_no.(task_nm).E3,2,'omitnan')];
    sem_ROI_min_choice_no_allSubs.(task_nm) = [sem(ROI_min_choice_no.(task_nm).E1,2),...
        sem(ROI_min_choice_no.(task_nm).E2,2),...
        sem(ROI_min_choice_no.(task_nm).E3,2)];
    % max yes
    m_ROI_max_choice_yes_allSubs.(task_nm) = [mean(ROI_max_choice_yes.(task_nm).E1,2,'omitnan'),...
        mean(ROI_max_choice_yes.(task_nm).E2,2,'omitnan'),...
        mean(ROI_max_choice_yes.(task_nm).E3,2,'omitnan')];
    sem_ROI_max_choice_yes_allSubs.(task_nm) = [sem(ROI_max_choice_yes.(task_nm).E1,2),...
        sem(ROI_max_choice_yes.(task_nm).E2,2),...
        sem(ROI_max_choice_yes.(task_nm).E3,2)];
    % max no
    m_ROI_max_choice_no_allSubs.(task_nm) = [mean(ROI_max_choice_no.(task_nm).E1,2,'omitnan'),...
        mean(ROI_max_choice_no.(task_nm).E2,2,'omitnan'),...
        mean(ROI_max_choice_no.(task_nm).E3,2,'omitnan')];
    sem_ROI_max_choice_no_allSubs.(task_nm) = [sem(ROI_max_choice_no.(task_nm).E1,2),...
        sem(ROI_max_choice_no.(task_nm).E2,2),...
        sem(ROI_max_choice_no.(task_nm).E3,2)];
end % task loop

%% figure
pSize = 30;
lWidth = 3;
xlim_vals = [0.75 3.25];

%% threshold
fig;
% Ep
subplot(1,2,1); hold on;
er_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_meanDelta_yn_allSubs.Ep,...
    sem_ROI_meanDelta_yn_allSubs.Ep);
er_hdl.LineWidth = lWidth;
er_hdl.Color = 'k';
xticks(1:n_hE_lvl);
xticklabels({'E1','E2','E3'});
xlabel('Physical effort level');
xlim(xlim_vals);
ylabel([ROI_short_nm,' YES - NO']);
legend_size(pSize);

% Em
subplot(1,2,2); hold on;
er_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_meanDelta_yn_allSubs.Em,...
    sem_ROI_meanDelta_yn_allSubs.Em);
er_hdl.LineWidth = lWidth;
er_hdl.Color = 'k';
refL = line(xlim_vals, [0 0]);
refL.LineWidth = 1;
refL.Color = 'k';
xticks(1:n_hE_lvl);
xticklabels({'E1','E2','E3'});
xlabel('Mental effort level');
xlim(xlim_vals);
ylabel([ROI_short_nm,' YES - NO']);
legend_size(pSize);

%% Look at Min vs Max
fig;
% Ep
subplot(1,2,1); hold on;
% yes choice
er_min_yes_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_min_choice_yes_allSubs.Ep,...
    sem_ROI_min_choice_yes_allSubs.Ep);
er_max_yes_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_max_choice_yes_allSubs.Ep,...
    sem_ROI_max_choice_yes_allSubs.Ep);
er_min_yes_hdl.LineWidth = lWidth;
er_min_yes_hdl.Color = 'g';
er_max_yes_hdl.LineWidth = lWidth;
er_max_yes_hdl.Color = 'g';
% no choice
er_min_no_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_min_choice_no_allSubs.Ep,...
    sem_ROI_min_choice_no_allSubs.Ep);
er_max_no_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_max_choice_no_allSubs.Ep,...
    sem_ROI_max_choice_no_allSubs.Ep);
er_min_no_hdl.LineWidth = lWidth;
er_min_no_hdl.Color = 'r';
er_max_no_hdl.LineWidth = lWidth;
er_max_no_hdl.Color = 'r';
refL = line(xlim_vals, [0 0]);
refL.LineWidth = 1;
refL.Color = 'k';
xticks(1:n_hE_lvl);
xticklabels({'E1','E2','E3'});
xlabel('Physical effort level');
xlim(xlim_vals);
ylabel(ROI_short_nm);
legend_size(pSize);

% Em
subplot(1,2,2); hold on;
% yes choice
er_min_yes_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_min_choice_yes_allSubs.Em,...
    sem_ROI_min_choice_yes_allSubs.Em);
er_max_yes_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_max_choice_yes_allSubs.Em,...
    sem_ROI_max_choice_yes_allSubs.Em);
er_min_yes_hdl.LineWidth = lWidth;
er_min_yes_hdl.Color = 'g';
er_max_yes_hdl.LineWidth = lWidth;
er_max_yes_hdl.Color = 'g';
% no choice
er_min_no_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_min_choice_no_allSubs.Em,...
    sem_ROI_min_choice_no_allSubs.Em);
er_max_no_hdl = errorbar(1:n_hE_lvl,...
    m_ROI_max_choice_no_allSubs.Em,...
    sem_ROI_max_choice_no_allSubs.Em);
er_min_no_hdl.LineWidth = lWidth;
er_min_no_hdl.Color = 'r';
er_max_no_hdl.LineWidth = lWidth;
er_max_no_hdl.Color = 'r';
refL = line(xlim_vals, [0 0]);
refL.LineWidth = 1;
refL.Color = 'k';
xticks(1:n_hE_lvl);
xticklabels({'E1','E2','E3'});
xlabel('Physical effort level');
xlim(xlim_vals);
ylabel(ROI_short_nm);
legend_size(pSize);