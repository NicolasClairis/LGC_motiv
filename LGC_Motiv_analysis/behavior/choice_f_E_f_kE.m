function[] = choice_f_E_f_kE(study_nm)
%[] = choice_f_E_f_kE(study_nm)
% choice_f_E_f_kE will look at the average proportion of choices per effort
% level and also by splitting on the kEp and kEm parameters from the model
%
%

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% working directory
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% main parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
fig_choice_f_E_f_kE_violin = fig;
fig_choice_f_E_f_kE_mean = fig;
n_hE_levels = 3;
x_hE = 1:n_hE_levels;
nRunsPerTask = 2;
low_col =[127 191 123]./255;
high_col =[175 141 195]./255;
pSize = 30;
lWidth = 3;

%% extract data
for iT = 1:nTasks
    task_nm = tasks{iT};
    switch task_nm
        case 'Ep'
            task_fullName = 'physical';
        case 'Em'
            task_fullName = 'mental';
    end
    
    [choice_hE.(task_nm).perRun_allSubs] = deal(NaN(n_hE_levels, NS, nRunsPerTask));
    [choice_hE.(task_nm).allSubs] = deal(NaN(n_hE_levels, NS));
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
        
        % average across runs
        for iE = 1:n_hE_levels
            choice_hE.(task_nm).allSubs(iE,iS) = mean(choice_hE.(task_nm).perRun_allSubs(iE,iS,:),3,'omitnan');
        end % effort level
    end % subject loop
    
    %% extract parameter
    prm_nm = ['k',task_nm];
    [low_prm_subs, high_prm_subs] = medSplit_prm(study_nm, subject_id, prm_nm);
    
    %% split people according to parameter
    choice_hE.(task_nm).(['low_',prm_nm]) = choice_hE.(task_nm).allSubs(:,low_prm_subs);
    choice_hE.(task_nm).(['high_',prm_nm]) = choice_hE.(task_nm).allSubs(:,high_prm_subs);
    %% average
    [m_choice_hE.(task_nm).(['low_',prm_nm]),...
        sem_choice_hE.(task_nm).(['low_',prm_nm])] = mean_sem_sd(choice_hE.(task_nm).(['low_',prm_nm]),2);
    [m_choice_hE.(task_nm).(['high_',prm_nm]),...
        sem_choice_hE.(task_nm).(['high_',prm_nm])] = mean_sem_sd(choice_hE.(task_nm).(['high_',prm_nm]),2);
    
    %% figure with violin plot
    figure(fig_choice_f_E_f_kE_violin);
    subplot(1,2,iT);
    hold on;
    choice_violin_mtrx = [choice_hE.(task_nm).(['low_',prm_nm])(1,:)',...
        choice_hE.(task_nm).(['high_',prm_nm])(1,:)',...
        choice_hE.(task_nm).(['low_',prm_nm])(2,:)',...
        choice_hE.(task_nm).(['high_',prm_nm])(2,:)',...
        choice_hE.(task_nm).(['low_',prm_nm])(3,:)',...
        choice_hE.(task_nm).(['high_',prm_nm])(3,:)'];
    % replace 0 by eps to avoid bugs
    choice_violin_mtrx(choice_violin_mtrx == 0) = eps;
    violinplot(choice_violin_mtrx,...
        {'E1_l','E1_h','E2_l','E2_h','E3_l','E3_h'},...
        'ViolinColor',[low_col;high_col;low_col;high_col;low_col;high_col]);
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);
    
    %% figure with mean
    figure(fig_choice_f_E_f_kE_mean);
    subplot(1,2,iT);
    hold on;
    low_hdl = errorbar(x_hE-0.02,...
        m_choice_hE.(task_nm).(['low_',prm_nm]),...
        sem_choice_hE.(task_nm).(['low_',prm_nm]));
    low_hdl.LineWidth = lWidth;
    low_hdl.Color = low_col;
    low_hdl.LineStyle = '-';
    high_hdl = errorbar(x_hE+0.02,...
        m_choice_hE.(task_nm).(['high_',prm_nm]),...
        sem_choice_hE.(task_nm).(['high_',prm_nm]));
    high_hdl.LineWidth = lWidth;
    high_hdl.Color = high_col;
    high_hdl.LineStyle = '-';
    legend([high_hdl, low_hdl],{['high ',prm_nm],['low ',prm_nm]});
    legend('boxoff');
    legend('Location','NorthEast');
    ylim([0 1]);
    ylabel([task_nm,' choices (%)']);
    legend_size(pSize);
end % task loop