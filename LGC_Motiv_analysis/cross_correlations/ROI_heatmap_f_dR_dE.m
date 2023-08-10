function[] = ROI_heatmap_f_dR_dE(study_nm, condition)
% ROI_heatmap_f_dR_dE will display a heatmap corresponding to the average
% level of activity of the ROI entered in input depending on the incentives
% and effort at stake for each effort type.
%

%% working directory
[computerRoot] = LGCM_root_paths();

%% extract list of subjects
if ~exist('study_nm','var') || ~isempty(study_nm)
    study_nm = 'study1';
end
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];
if ~exist('condition','var') || ~isempty(condition)
    condition = subject_condition;
end
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);
[condition_for_fMRI_extraction] = condition_for_fMRI(condition);

%% select and extract the ROI data
[ROI_trial_b_trial, ROI_subList,...
    fMRI_ROI_nm, ROI_short_nm,...
    ~, timePeriod_nm] = extract_ROI_betas_onsets_only_bis(computerRoot,...
    study_nm, subject_id, condition_for_fMRI_extraction);

%% prepare output of interest
task_names = {'Ep','Em'};
nTasks = length(task_names);
n_dR = 3;
n_dP = 3;
n_dInc = n_dR + n_dP;
n_dE = 3;
[ROI_avg_perSub.Ep, choice_avg_perSub.Ep, RT_avg_perSub.Ep,...
    ROI_avg_perSub.Em, choice_avg_perSub.Em, RT_avg_perSub.Em] = deal(NaN(n_dInc, n_dE, NS));
nRunsPerTask = 2;
[ROI_avg_perSubperRun.Ep, choice_avg_perSubperRun.Ep, RT_avg_perSubperRun.Ep,...
    ROI_avg_perSubperRun.Em, choice_avg_perSubperRun.Em, RT_avg_perSubperRun.Em] = deal(NaN(n_dInc, n_dE, NS, nRunsPerTask));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder,...
        'CID',sub_nm, filesep, 'behavior',filesep];
    
    % extract runs
    [runsStruct, n_runs] = runs_definition(study_nm, sub_nm, condition);
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    %% loop through runs
    for iRun = 1:n_runs
        kRun = okRuns(iRun);
        task_nm_tmp = taskNames{iRun};
        run_nm = num2str(kRun);
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        switch kRun
            case {1,2}
                jRun = 1;
            case {3,4}
                jRun = 2;
        end
        run_nm_bis = ['run',num2str(jRun)];
        
        %% extract choices for the current session
        choice_hE_allTrials_tmp = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        R_level_tmp = extract_hR_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        P_level_tmp = extract_hP_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        E_level_tmp = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        RT_allTrials_tmp = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        %% extract fMRI ROI
        fMRI_allTrials_tmp = ROI_trial_b_trial.(fMRI_ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
        
        %% average per incentive and effort level
        for iE = 1:n_dE
            jIncentiveLine = 0;
            for iP = n_dP:(-1):1 % revert order to have symetry with rewards
                jIncentiveLine = jIncentiveLine + 1;
                PE_trials_idx = (E_level_tmp == iE).*(P_level_tmp == iP) == 1;
                ROI_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, jRun) = mean(fMRI_allTrials_tmp(PE_trials_idx),'omitnan');
                choice_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, jRun) = mean(choice_hE_allTrials_tmp(PE_trials_idx),'omitnan');
                RT_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, jRun) = mean(RT_allTrials_tmp(PE_trials_idx),'omitnan');
            end % punishment
            
            for iR = 1:n_dR
                jIncentiveLine = jIncentiveLine + 1;
                RE_trials_idx = (E_level_tmp == iE).*(R_level_tmp == iR) == 1;
                ROI_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, jRun) = mean(fMRI_allTrials_tmp(RE_trials_idx),'omitnan');
                choice_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, jRun) = mean(choice_hE_allTrials_tmp(RE_trials_idx),'omitnan');
                RT_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, jRun) = mean(RT_allTrials_tmp(RE_trials_idx),'omitnan');
            end % reward
        end % effort level
    end % run loop
    
    %% average runs together
    for iTask = 1:nTasks
        task_nm = task_names{iTask};
        ROI_avg_perSub.(task_nm)(:,:,iS) = mean(ROI_avg_perSubperRun.(task_nm)(:, :, iS, :), 4,'omitnan');
        choice_avg_perSub.(task_nm)(:,:,iS) = mean(choice_avg_perSubperRun.(task_nm)(:, :, iS, :), 4,'omitnan');
        RT_avg_perSub.(task_nm)(:,:,iS) = mean(RT_avg_perSubperRun.(task_nm)(:, :, iS, :), 4,'omitnan');
    end % task loop
    
    % indicate where we are at
    disp(['subject ',num2str(iS),'/',num2str(NS),' done']);
end % subject loop

%% average subjects together
for iTask = 1:nTasks
    task_nm = task_names{iTask};
    ROI_avg.(task_nm) = mean(ROI_avg_perSub.(task_nm), 3,'omitnan');
    choice_avg.(task_nm) = mean(choice_avg_perSub.(task_nm), 3,'omitnan');
    RT_avg.(task_nm) = mean(RT_avg_perSub.(task_nm), 3,'omitnan');
end % task loop

%% figures
nPlotsPerLine = 3;
[pSize] = general_fig_prm;
ROI_range = [min(min(min(ROI_avg.Ep, ROI_avg.Em))),...
    max(max(max(ROI_avg.Ep, ROI_avg.Em)))];
choice_range = [0 100];
RT_range = [1.8 2.4];
fig;
% define which colormap you want to use (see full list here if you are not
% happy with the selection:
% https://ch.mathworks.com/help/matlab/ref/colormap.html)
% color_range_ROI_choices = 'hot';
% color_range_ROI_choices = 'turbo';
color_range_ROI_choices = redblue(45);

for iTask = 1:nTasks
    task_nm = task_names{iTask};
    % ROI
    ROI_hdl = subplot(nTasks, nPlotsPerLine, 1 + nPlotsPerLine*(iTask - 1));
    imagesc(ROI_avg.(task_nm),ROI_range);
%     imagesc(ROI_avg.(task_nm));
    colormap(ROI_hdl, color_range_ROI_choices);
    xticks(1:n_dE);
    xticklabels({'E1','E2','E3'});
    yticks(1:n_dInc)
    yticklabels({'P3','P2','P1','R1','R2','R3'});
    title(ROI_short_nm);
    legend_size(pSize);
    colorbar;
    
    % choices
    choice_hdl = subplot(nTasks, nPlotsPerLine, 2 + nPlotsPerLine*(iTask - 1));
    imagesc(choice_avg.(task_nm).*100, choice_range);
    colormap(choice_hdl, color_range_ROI_choices);
    xticks(1:n_dE);
    xticklabels({'E1','E2','E3'});
    yticks(1:n_dInc)
    yticklabels({'P3','P2','P1','R1','R2','R3'});
    title('Choice (%)');
    legend_size(pSize);
    colorbar;
    
    % RT
    RT_hdl = subplot(nTasks, nPlotsPerLine, 3 + nPlotsPerLine*(iTask - 1));
    imagesc(RT_avg.(task_nm), RT_range);
    colormap(RT_hdl, color_range_ROI_choices);
    % revert axis for RT to have warmer colours for faster answers
    RT_hdl.Colormap = flipud(RT_hdl.Colormap);
    xticks(1:n_dE);
    xticklabels({'E1','E2','E3'});
    yticks(1:n_dInc)
    yticklabels({'P3','P2','P1','R1','R2','R3'});
    title('RT (s)');
    legend_size(pSize);
    colorbar;
end % task loop
end % function