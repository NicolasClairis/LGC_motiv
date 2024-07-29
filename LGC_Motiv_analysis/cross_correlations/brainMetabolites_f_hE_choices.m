% function[] = metabolites_f_hE_choices()
% metabolites_f_hE_choices will split participants depending on the
% proportion of high effort choices made for the high effort option (effort
% level 3) within each task. Then, the script will look at the different
% levels of metabolite depending on this.

%% working directory
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% subject identification
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% loop through subjects
tasks = {'Ep','Em'};
nTasks = length(tasks);
nTrialsPerRun = 54;
[totalChoice_hE.Ep, totalChoice_hE.Em] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];

    %% extract runs
    [runsStruct] = runs_definition(study_nm, sub_nm, condition);
    nRuns = length(runsStruct.tasks);
    runs_Ep = strcmp(runsStruct.tasks,'Ep');
    runs_Em = strcmp(runsStruct.tasks,'Em');

    %% loop through tasks
    for iT = 1:nTasks
        task_nm = tasks{iT};
        switch task_nm
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end
        nRunsPerTask = sum(strcmp(runsStruct.tasks,task_nm));
        totalChoice_hE_perRun_tmp = NaN(2,1);

        %% loop through runs
        jRun = 0;
        for iRun = 1:nRuns
            runToInclude = 0;
            switch task_nm
                case 'Ep'
                    if runs_Ep(iRun) == 1
                        jRun = jRun + 1;
                        runToInclude = 1;
                    end
                case 'Em'
                    if runs_Em(iRun) == 1
                        jRun = jRun + 1;
                        runToInclude = 1;
                    end
            end
            if runToInclude == 1
                run_nm = num2str(runsStruct.runsToKeep(iRun));
                [choice_highE_tmp] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                [hE_level] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                totalChoice_hE_perRun_tmp(jRun) = (sum(choice_highE_tmp(hE_level == 3), 'omitnan')./sum(hE_level == 3)).*100;
            end
        end % run loop

        totalChoice_hE.(task_nm)(iS) = mean(totalChoice_hE_perRun_tmp,1,'omitnan');
    end % task loop
end % subject loop

%% loop through tasks
MRS_ROIs = {'dmPFC','aIns'};
n_MRS_roi = length(MRS_ROIs);
for iT = 1:nTasks
    task_nm = tasks{iT};
    %% perform median split based on proportion of choices
    med_choice_hE.(task_nm) = median(totalChoice_hE.(task_nm),'omitnan');
    low_hE_mSplit_subs = totalChoice_hE.(task_nm) <= med_choice_hE.(task_nm);
    high_hE_mSplit_subs = totalChoice_hE.(task_nm) > med_choice_hE.(task_nm);

    %% alternative split people based on a 50% threshold
    low_hE_50Split_subs = totalChoice_hE.(task_nm) <= 50;
    high_hE_50Split_subs = totalChoice_hE.(task_nm) > 50;

    %% load all metabolites
    [allMetabolitesStruct] = metabolite_load(subject_id);
    %% extract metabolites for each group
    for iMRS_roi = 1:n_MRS_roi
        MRS_ROI_nm = MRS_ROIs{iMRS_roi};
        metabolite_names = fieldnames(allMetabolitesStruct.(MRS_ROI_nm));
        for iMb = 1:length(metabolite_names)
            mb_nm = metabolite_names{iMb};
            mb_nm_bis = [MRS_ROI_nm,'_',mb_nm];

            % zscore all metabolites across subjects to have everyone in a
            % common scale + more normal values
            metabolite_values = nanzscore(allMetabolitesStruct.(MRS_ROI_nm).(mb_nm));

            % median split
            [m_metabolites.(task_nm).mSplit.(mb_nm_bis).hE,...
                sem_metabolites.(task_nm).mSplit.(mb_nm_bis).hE,...
                sd_metabolites.(task_nm).mSplit.(mb_nm_bis).hE] = mean_sem_sd(metabolite_values(high_hE_mSplit_subs),2);
            [m_metabolites.(task_nm).mSplit.(mb_nm_bis).lE,...
                sem_metabolites.(task_nm).mSplit.(mb_nm_bis).lE,...
                sd_metabolites.(task_nm).mSplit.(mb_nm_bis).lE] = mean_sem_sd(metabolite_values(low_hE_mSplit_subs),2);
            [~,pval.(task_nm).mSplit.(mb_nm_bis)] = ttest2(metabolite_values(low_hE_mSplit_subs),...
                metabolite_values(high_hE_mSplit_subs));
            if pval.(task_nm).mSplit.(mb_nm_bis) < 0.05
                pval.signif.(task_nm).mSplit.(mb_nm_bis) = pval.(task_nm).mSplit.(mb_nm_bis);
            end % significant p.value

            % split on 50% choices
            [m_metabolites.(task_nm).midSplit.(mb_nm_bis).hE,...
                sem_metabolites.(task_nm).midSplit.(mb_nm_bis).hE,...
                sd_metabolites.(task_nm).midSplit.(mb_nm_bis).hE] = mean_sem_sd(metabolite_values(high_hE_50Split_subs),2);
            [m_metabolites.(task_nm).midSplit.(mb_nm_bis).lE,...
                sem_metabolites.(task_nm).midSplit.(mb_nm_bis).lE,...
                sd_metabolites.(task_nm).midSplit.(mb_nm_bis).lE] = mean_sem_sd(metabolite_values(low_hE_50Split_subs),2);
            [~,pval.(task_nm).midSplit.(mb_nm_bis)] = ttest2(metabolite_values(low_hE_50Split_subs),...
                metabolite_values(high_hE_50Split_subs));
            if pval.(task_nm).midSplit.(mb_nm_bis) < 0.05
                pval.signif.(task_nm).midSplit.(mb_nm_bis) = pval.(task_nm).midSplit.(mb_nm_bis);
            end % significant p.value
        end % metabolite loop
    end % MRS ROI
    
    %% figure
    allMb_names = fieldnames(m_metabolites.(task_nm).mSplit);
    nTotalMb = length(allMb_names);
    pSize = 20;
    lWidth = 3;
    orange = [230 85 13]./255;
    lightBlue = [49 130 189]./255;

    % median split
    fig;
    hold on;
    for iMb = 1:nTotalMb
        metab_nm = allMb_names{iMb};
        low_hdl = errorbar(iMb-0.02,...
            m_metabolites.(task_nm).mSplit.(metab_nm).lE,...
            sem_metabolites.(task_nm).mSplit.(metab_nm).lE);
        high_hdl = errorbar(iMb+0.02,...
            m_metabolites.(task_nm).mSplit.(metab_nm).hE,...
            sem_metabolites.(task_nm).mSplit.(metab_nm).hE);
        % shape
        low_hdl.LineWidth = lWidth;
        low_hdl.Color = orange;
        high_hdl.LineWidth = lWidth;
        high_hdl.Color = lightBlue;

        % add a line on top if the p.value is significant
        if pval.(task_nm).mSplit.(metab_nm) < 0.05
            yLine = max( (m_metabolites.(task_nm).mSplit.(metab_nm).lE+sem_metabolites.(task_nm).mSplit.(metab_nm).lE),...
                (m_metabolites.(task_nm).mSplit.(metab_nm).hE+sem_metabolites.(task_nm).mSplit.(metab_nm).hE)) +...
                +sem_metabolites.(task_nm).mSplit.(metab_nm).lE*3;
            signifLine = line([iMb-0.4, iMb+0.4],...
                [yLine, yLine]);
            signifLine.LineWidth = lWidth*2;
            signifLine.Color = 'k';
            signifLine.LineStyle = '-';
        end
    end
    xticks(1:nTotalMb);
    xticklabels(strrep(allMb_names,'_',' '));
    ylabel(['Metabolite concentration (zscored) - ',task_nm]);
    legend_size(pSize);

    % 50% split
    fig;
    hold on;
    for iMb = 1:nTotalMb
        metab_nm = allMb_names{iMb};
        low_hdl = errorbar(iMb-0.02,...
            m_metabolites.(task_nm).midSplit.(metab_nm).lE,...
            sem_metabolites.(task_nm).midSplit.(metab_nm).lE);
        high_hdl = errorbar(iMb+0.02,...
            m_metabolites.(task_nm).midSplit.(metab_nm).hE,...
            sem_metabolites.(task_nm).midSplit.(metab_nm).hE);
        % shape
        low_hdl.LineWidth = lWidth;
        low_hdl.Color = lightBlue;
        high_hdl.LineWidth = lWidth;
        high_hdl.Color = orange;

        % add a line on top if the p.value is significant
        if pval.(task_nm).midSplit.(metab_nm) < 0.05
            yLine = max( (m_metabolites.(task_nm).midSplit.(metab_nm).lE+sem_metabolites.(task_nm).midSplit.(metab_nm).lE),...
                (m_metabolites.(task_nm).midSplit.(metab_nm).hE+sem_metabolites.(task_nm).midSplit.(metab_nm).hE)) +...
                +sem_metabolites.(task_nm).midSplit.(metab_nm).lE*3;
            signifLine = line([iMb-0.4, iMb+0.4],...
                [yLine, yLine]);
            signifLine.LineWidth = lWidth*2;
            signifLine.Color = 'k';
            signifLine.LineStyle = '-';
        end
    end
    xticks(1:nTotalMb);
    xticklabels(strrep(allMb_names,'_',' '));
    ylabel(['Metabolite concentration (zscored) - ',task_nm]);
    legend_size(pSize);
end % task loop
% end % function