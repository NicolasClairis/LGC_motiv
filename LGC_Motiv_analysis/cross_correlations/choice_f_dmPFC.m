% use GLM 102 and 103 where the BOLD is extracted either across
% all trials (GLM103) or separately per effort level (GLM102) allowing to
% look whether the average level of activity predicts the average
% percentage of choice per run, per task, or globally.

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% main parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
runsPerTask = 2;
E_tasks = {'E','Ep','Em'};
n_E_tasks = length(E_tasks);

%% extract ROI level
% across all trials
GLM = 103;
[con_vec_all_GLM103,...
    con_avg_GLM103, con_sem_GLM103, con_sd_GLM103,...
    con_names_GLM103,...
    ROI_coords_GLM103, ttest_ROI_GLM103] = ROI_extraction_group(study_nm, GLM,...
    subject_id, condition, 0);

% per effort level
GLM = 102;
[con_vec_all_GLM102,...
    con_avg_GLM102, con_sem_GLM102, con_sd_GLM102,...
    con_names_GLM102,...
    ROI_coords_GLM102, ttest_ROI_GLM102] = ROI_extraction_group(study_nm, GLM,...
    subject_id, condition, 0);

%% extract choice proportion per subject/task/run
E_trialType = {'E','E1','E2','E3'};
for iE = 1:length(E_trialType)
    E_nm = E_trialType{iE};
    [choice.E.(E_nm).all,...
        choice.Ep.(E_nm).all, choice.Em.(E_nm).all,...
        choice.E.(E_nm).run1, choice.E.(E_nm).run2,...
        choice.E.(E_nm).run3, choice.E.(E_nm).run4,...
        choice.Ep.(E_nm).run1, choice.Ep.(E_nm).run2,...
        choice.Em.(E_nm).run1, choice.Em.(E_nm).run2,...
        fMRI_BOLD.E.(E_nm).all,...
        fMRI_BOLD.Ep.(E_nm).all, fMRI_BOLD.Em.(E_nm).all] = deal(NaN(1,NS));
end % effort type loop

for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
    runsStruct = runs_definition(study_nm, sub_nm, condition);
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    for iRun = 1:length(okRuns)
        kRun = okRuns(iRun);
        task_nm_tmp = taskNames{iRun};
        run_nm = num2str(kRun);
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
        
        % global extraction
        choice_hE_tmp = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        choice.E.E.(['run',run_nm])(iS) = sum(choice_hE_tmp==1)./sum(ismember(choice_hE_tmp,[0,1]));
        choice.(task_nm_tmp).E.(run_nm_bis)(iS) = sum(choice_hE_tmp==1)./sum(ismember(choice_hE_tmp,[0,1]));
        
        % extraction per effort level
        hE_lvl_tmp = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        for iE = 1:3
            E_nm = ['E',num2str(iE)];
            hE_idx = hE_lvl_tmp == iE;
            choice_hE_perE_tmp = choice_hE_tmp(hE_idx);
            choice.E.(E_nm).(['run',run_nm])(iS) = sum(choice_hE_perE_tmp==1)./sum(ismember(choice_hE_perE_tmp,[0,1]));
            choice.(task_nm_tmp).(E_nm).(run_nm_bis)(iS) = sum(choice_hE_perE_tmp==1)./sum(ismember(choice_hE_perE_tmp,[0,1]));
        end % effort level
    end % run loop
    
    %% average per subject across runs
    % average choices across runs
    choice.E.E.all(iS) = mean([choice.E.E.run1(iS), choice.E.E.run2(iS),...
        choice.E.E.run3(iS), choice.E.E.run4(iS)],'omitnan');
    choice.Ep.E.all(iS) = mean([choice.Ep.E.run1(iS),...
        choice.Ep.E.run2(iS)],'omitnan');
    choice.Em.E.all(iS) = mean([choice.Em.E.run1(iS),...
        choice.Em.E.run2(iS)],'omitnan');
    % average fMRI across runs
    E_con_idx_GLM103 = strcmp(con_names_GLM103,'Ep+Em ONSET choice RP E');
    fMRI_BOLD.E.E.all(iS) = con_vec_all_GLM103(E_con_idx_GLM103,iS);
    Ep_con_idx_GLM103 = strcmp(con_names_GLM103,'Ep ONSET choice RP E');
    fMRI_BOLD.Ep.E.all(iS) = con_vec_all_GLM103(Ep_con_idx_GLM103,iS);
    Em_con_idx_GLM103 = strcmp(con_names_GLM103,'Em ONSET choice RP E');
    fMRI_BOLD.Em.E.all(iS) = con_vec_all_GLM103(Em_con_idx_GLM103,iS);
    
    % average for each effort level
    for iE = 1:3
        E_nm = ['E',num2str(iE)];
        % average choices across runs
        choice.E.(E_nm).all(iS) = mean([choice.E.(E_nm).run1(iS),...
            choice.E.(E_nm).run2(iS),...
            choice.E.(E_nm).run3(iS),...
            choice.E.(E_nm).run4(iS)],'omitnan');
        choice.Ep.(E_nm).all(iS) = mean([choice.Ep.(E_nm).run1(iS),...
            choice.Ep.(E_nm).run2(iS)],'omitnan');
        choice.Em.(E_nm).all(iS) = mean([choice.Em.(E_nm).run1(iS),...
            choice.Em.(E_nm).run2(iS)],'omitnan');
        % average fMRI across runs
        E_con_idx_GLM102 = strcmp(con_names_GLM102,['Ep+Em ONSET choice RP ',E_nm]);
        fMRI_BOLD.E.(E_nm).all(iS) = con_vec_all_GLM102(E_con_idx_GLM102,iS);
        Ep_con_idx_GLM102 = strcmp(con_names_GLM102,['Ep ONSET choice RP ',E_nm]);
        fMRI_BOLD.Ep.(E_nm).all(iS) = con_vec_all_GLM102(Ep_con_idx_GLM102,iS);
        Em_con_idx_GLM102 = strcmp(con_names_GLM102,['Em ONSET choice RP ',E_nm]);
        fMRI_BOLD.Em.(E_nm).all(iS) = con_vec_all_GLM102(Em_con_idx_GLM102,iS);
    end % effort level
end % subject loop

%% test across subjects
% choices across all trials = f(dmPFC)
for iE_task = 1:n_E_tasks
    E_task_nm = E_tasks{iE_task};
    for iE_lvl = 1:length(E_trialType)
        E_lvl_nm = E_trialType{iE_lvl};
        [betas.(E_task_nm).(E_lvl_nm).all, ~,stats_tmp] = glmfit(fMRI_BOLD.(E_task_nm).(E_lvl_nm).all,...
            choice.(E_task_nm).(E_lvl_nm).all,'normal');
        pval.(E_task_nm).(E_lvl_nm).all = stats_tmp.p;
        fMRI_BOLD_sorted.(E_task_nm).(E_lvl_nm).all = sort(fMRI_BOLD.(E_task_nm).(E_lvl_nm).all);
        choice_fit.(E_task_nm).(E_lvl_nm).all = glmval(betas.(E_task_nm).(E_lvl_nm).all,...
            fMRI_BOLD_sorted.(E_task_nm).(E_lvl_nm).all, 'identity');
    end % effort level E/E1/E2/E3
end % effort task E/Ep/Em


%% figures
grey = [143 143 143]./255;
lWidth = 3;
pSize = 30;

for iE_lvl = 1:length(E_trialType)
    E_lvl_nm = E_trialType{iE_lvl};
    switch E_lvl_nm
        case 'E'
            E_lvl_nm_bis = '';
        otherwise
            E_lvl_nm_bis = E_lvl_nm;
    end
    
    % choice = f(dmPFC)
    fig;
    % all efforts
    subplot(2,2,1:2);
    hold on;
    scat_hdl = scatter(fMRI_BOLD.E.(E_lvl_nm).all,...
        choice.E.(E_lvl_nm).all);
    fit_hdl = plot(fMRI_BOLD_sorted.E.(E_lvl_nm).all,...
        choice_fit.E.(E_lvl_nm).all);
    scat_hdl.MarkerEdgeColor = 'k';
    fit_hdl.LineStyle = '--';
    fit_hdl.Color = grey;
    fit_hdl.LineWidth = lWidth;
    xlabel(['dmPFC BOLD during ',E_lvl_nm_bis,' choice']);
    ylabel([E_lvl_nm_bis,' choices (%)']);
    legend_size(pSize);
    
    % Ep
    subplot(2,2,3);
    hold on;
    scat_hdl = scatter(fMRI_BOLD.Ep.(E_lvl_nm).all,...
        choice.Ep.(E_lvl_nm).all);
    fit_hdl = plot(fMRI_BOLD_sorted.Ep.(E_lvl_nm).all,...
        choice_fit.Ep.(E_lvl_nm).all);
    scat_hdl.MarkerEdgeColor = 'k';
    fit_hdl.LineStyle = '--';
    fit_hdl.Color = grey;
    fit_hdl.LineWidth = lWidth;
    xlabel(['dmPFC BOLD during ',E_lvl_nm_bis,' choice']);
    ylabel(['Ep ',E_lvl_nm_bis,' choices (%)']);
    legend_size(pSize);
    
    % Em
    subplot(2,2,4);
    hold on;
    scat_hdl = scatter(fMRI_BOLD.Em.(E_lvl_nm).all,...
        choice.Em.(E_lvl_nm).all);
    fit_hdl = plot(fMRI_BOLD_sorted.Em.(E_lvl_nm).all,...
        choice_fit.Em.(E_lvl_nm).all);
    scat_hdl.MarkerEdgeColor = 'k';
    fit_hdl.LineStyle = '--';
    fit_hdl.Color = grey;
    fit_hdl.LineWidth = lWidth;
    xlabel(['dmPFC BOLD during ',E_lvl_nm_bis,' choice']);
    ylabel(['Em ',E_lvl_nm_bis,' choices (%)']);
    legend_size(pSize);
    
end % effort