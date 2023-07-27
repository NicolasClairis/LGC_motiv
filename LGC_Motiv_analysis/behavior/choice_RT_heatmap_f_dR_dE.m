function[] = choice_RT_heatmap_f_dR_dE(study_nm, condition)
% choice_RT_heatmap_f_dR_dE will display a heatmap corresponding to the average
% proportion of effortful choices and reaction times (RT) depending on the 
% incentives and effort at stake for each effort type.
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

%% prepare output of interest
task_names = {'Ep','Em'};
nTasks = length(task_names);
n_dR = 3;
n_dP = 3;
n_dInc = n_dR + n_dP;
n_dE = 3;
[choice_avg_perSub.Ep,...
    RT_avg_perSub.Ep, RT_fit_avg_perSub.Ep,...
    choice_avg_perSub.Em,...
    RT_avg_perSub.Em, RT_fit_avg_perSub.Em] = deal(NaN(n_dInc, n_dE, NS));
nRunsPerTask = 2;
[choice_avg_perSubperRun.Ep,...
    RT_avg_perSubperRun.Ep, RT_fit_avg_perSubperRun.Ep,...
    choice_avg_perSubperRun.Em,...
    RT_avg_perSubperRun.Em, RT_fit_avg_perSubperRun.Em] = deal(NaN(n_dInc, n_dE, NS, nRunsPerTask));
RT_fit_GLM = 2; % GLM to use for RT fit extraction

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder,...
        'CID',sub_nm, filesep, 'behavior',filesep];
    
    % extract RT fit for this particular subject
    [~, ~, ~, ~,...
        RT_fit_perRun_tmp] = RT_GLM(0, computerRoot, study_nm, {sub_nm}, RT_fit_GLM);
    
    % extract runs
    [runsStruct, n_runs] = runs_definition(study_nm, sub_nm, condition);
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    %% loop through runs
    for iRun = 1:n_runs
        kRun = okRuns(iRun);
        task_nm_tmp = taskNames{iRun};
        run_nm = num2str(kRun);
        run_fullNm = ['run',num2str(run_nm)];
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        
        %% extract choices for the current session
        choice_hE_allTrials_tmp = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        R_level_tmp = extract_hR_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        P_level_tmp = extract_hP_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        E_level_tmp = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        RT_allTrials_tmp = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        RT_fit_allTrials_tmp = RT_fit_perRun_tmp.(run_fullNm);
        
        %% average per incentive and effort level
        for iE = 1:n_dE
            jIncentiveLine = 0;
            for iP = n_dP:(-1):1 % revert order to have symetry with rewards
                jIncentiveLine = jIncentiveLine + 1;
                PE_trials_idx = (E_level_tmp == iE).*(P_level_tmp == iP) == 1;
                choice_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, kRun) = mean(choice_hE_allTrials_tmp(PE_trials_idx),'omitnan');
                RT_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, kRun) = mean(RT_allTrials_tmp(PE_trials_idx),'omitnan');
                RT_fit_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, kRun) = mean(RT_fit_allTrials_tmp(PE_trials_idx),'omitnan');
            end % punishment
            
            for iR = 1:n_dR
                jIncentiveLine = jIncentiveLine + 1;
                RE_trials_idx = (E_level_tmp == iE).*(R_level_tmp == iR) == 1;
                choice_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, kRun) = mean(choice_hE_allTrials_tmp(RE_trials_idx),'omitnan');
                RT_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, kRun) = mean(RT_allTrials_tmp(RE_trials_idx),'omitnan');
                RT_fit_avg_perSubperRun.(task_nm_tmp)(jIncentiveLine, iE, iS, kRun) = mean(RT_fit_allTrials_tmp(RE_trials_idx),'omitnan');
            end % reward
        end % effort level
    end % run loop
    
    %% average runs together
    for iTask = 1:nTasks
        task_nm = task_names{iTask};
        choice_avg_perSub.(task_nm)(:,:,iS) = mean(choice_avg_perSubperRun.(task_nm)(:, :, iS, :), 4,'omitnan');
        RT_avg_perSub.(task_nm)(:,:,iS) = mean(RT_avg_perSubperRun.(task_nm)(:, :, iS, :), 4,'omitnan');
        RT_fit_avg_perSub.(task_nm)(:,:,iS) = mean(RT_fit_avg_perSubperRun.(task_nm)(:, :, iS, :), 4,'omitnan');
    end % task loop
    
    % indicate where we are at
    disp(['subject ',num2str(iS),'/',num2str(NS),' done']);
end % subject loop

%% average subjects together
for iTask = 1:nTasks
    task_nm = task_names{iTask};
    choice_avg.(task_nm) = mean(choice_avg_perSub.(task_nm), 3,'omitnan');
    RT_avg.(task_nm) = mean(RT_avg_perSub.(task_nm), 3,'omitnan');
    RT_fit_avg.(task_nm) = mean(RT_fit_avg_perSub.(task_nm), 3,'omitnan');
end % task loop

%% figures
nLines = 2;
[pSize] = general_fig_prm;
choice_range = [0 1];
RT_range = [1 1.5];
choice_fig = fig;
RT_fig = fig;

for iTask = 1:nTasks
    task_nm = task_names{iTask};
    
    %% choice figure
    figure(choice_fig);
    iChoiceLine = 1;
    choice_plot_hdl = subplot(nLines, nTasks, iTask + nTasks*(iChoiceLine - 1));
    imagesc(choice_avg.(task_nm), choice_range);
    % define which colormap you want to use (see full list here if you are not
    % happy with the selection:
    % https://ch.mathworks.com/help/matlab/ref/colormap.html)
    % colormap hot;
    % colormap jet;
    colormap(choice_plot_hdl,redblue(45));
    xticks(1:n_dE);
    xticklabels({'E1','E2','E3'});
    yticks(1:n_dInc)
    yticklabels({'P3','P2','P1','R1','R2','R3'});
    title('Choice (%)');
    legend_size(pSize);
    colorbar;
    
    % should add bayesian fit here
    iChoiceLine = 2;
    warning('please find a way to include Arthur''s bayesian model results here');
    
    %% RT figures
    figure(RT_fig);
    iRTline = 1;
    RT_plot_hdl = subplot(nLines, nTasks, iTask +  + nTasks*(iRTline - 1));
    imagesc(RT_avg.(task_nm), RT_range);
    % define which colormap you want to use (see full list here if you are not
    % happy with the selection:
    % https://ch.mathworks.com/help/matlab/ref/colormap.html)
    % colormap hot;
    % colormap jet;
    % revert color axis for RT (so that faster is shown with hoter colors)
    colormap(RT_plot_hdl,flipud(redblue(45)));
    xticks(1:n_dE);
    xticklabels({'E1','E2','E3'});
    yticks(1:n_dInc)
    yticklabels({'P3','P2','P1','R1','R2','R3'});
    title('RT (s)');
    legend_size(pSize);
    colorbar;
    
    % RT fit
    iRTline = 2;
    RT_fit_plot_hdl = subplot(nLines, nTasks, iTask +  + nTasks*(iRTline - 1));
    imagesc(RT_fit_avg.(task_nm), RT_range);
    % define which colormap you want to use (see full list here if you are not
    % happy with the selection:
    % https://ch.mathworks.com/help/matlab/ref/colormap.html)
    % colormap hot;
    % colormap jet;
    % revert color axis for RT (so that faster is shown with hoter colors)
    colormap(RT_fit_plot_hdl,flipud(redblue(45)));
    xticks(1:n_dE);
    xticklabels({'E1','E2','E3'});
    yticks(1:n_dInc)
    yticklabels({'P3','P2','P1','R1','R2','R3'});
    title('RT (s)');
    legend_size(pSize);
    colorbar;
end % task loop
end % function