function[betas, pval] = mRT_f_parameters()
% [betas, pval] = mRT_f_parameters()
% mRT_f_parameters will look whether mean (or median or standard deviation)
% of the reaction times correlate with sensitivity to effort across 
% participants.
%
% OUTPUTS
% betas: structure with betas for the slope between each parameter and
% mean/median RT
%
% pval: structure with p.value corresponding to the relationship between
% each parameter and mean/median RT

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% subject selection
study_nm = 'study1';
condition = subject_condition();
gender = 'all';
[subject_id, NS] = LGCM_subject_selection('study1',condition,gender);

%% select RT type
RT_types = {'raw','log'};
RT_type_idx = listdlg('PromptString','Select RT format','ListString',RT_types,'SelectionMode','single');
RT_type = RT_types{RT_type_idx};
switch RT_type
    case 'raw'
        raw_z_log = @(x) x;
    case 'log'
        raw_z_log = @(x) log(x);
end
% careful not to zscore as that would remove any inter-individual
% difference by normalizing everybody

%% by default display figure
if ~exist('figDisp','var') || isempty(figDisp) || ~islogical(figDisp)
    figDisp = true;
end

%% define main variables
nTrialsPerRun = 54;
tasks = {'Ep','Em'};
nTasks = length(tasks);
[mean_RT.Ep, mean_RT.Em,...
    median_RT.Ep, median_RT.Em,...
     sd_RT.Ep, sd_RT.Em] = deal(NaN(1,NS));

%% extract behavioral parameters
[prm] = prm_extraction(study_nm, subject_id, [], []);
parameters = fieldnames(prm);
parameters(strcmp(parameters,'CID')) = [];
nPrm = length(parameters);

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];
    
    %% extract average choice RT across all tasks
    % extract information of runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    nTotalRuns = length(runsStruct.tasks);
    for iTask = 1:nTasks
        task_id = tasks{iTask};
        switch task_id
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end
        runs.(task_id) = strcmp(runsStruct.tasks,task_id);
        nRuns.(task_id) = sum(runs.(task_id));
        nSubjectTotalTrials.(task_id) = nTrialsPerRun*nRuns.(task_id);
        % prepare data to be extracted across runs
        [RT_perTrial_tmp.(task_id),...
            inc_perTrial_tmp.(task_id),...
            eff_perTrial_tmp.(task_id),...
            choice_hE_perTrial_tmp.(task_id)] = deal(NaN(nSubjectTotalTrials.(task_id), 1));
        
        jRun = 0;
        for iRun = 1:nTotalRuns
            runToInclude = 0;
            if runs.(task_id)(iRun) == 1
                jRun = jRun + 1;
                runToInclude = 1;
            end
            run_nm = num2str(iRun);
            
            if runToInclude == 1
                runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);
                %% load data
                behaviorStruct_tmp = load([subBehaviorFolder,...
                    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
                    '_task.mat']);
                % load relevant data
                choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
                switch task_id
                    case 'Ep'
                        choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                    case 'Em'
                        choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
                end
                % remove confidence information from choice
                choice_LR_tmp(choice_LR_tmp == 2) = 1;
                choice_LR_tmp(choice_LR_tmp == -2) = -1;
                % extract RT
                [RT_tmp] = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                RT_perTrial_tmp.(task_id)(runTrials_idx) = RT_tmp;
                % extract choice made
                defaultSide_tmp = choiceOptions_tmp.default_LR;
                choice_hE_tmp = choice_LR_tmp == -defaultSide_tmp; % high effort chosen = non-default chosen
                choice_hE_tmp = double(choice_hE_tmp);
                choice_hE_tmp(choice_LR_tmp == 0) = NaN; % place NaN when no choice was made
                choice_hE_perTrial_tmp.(task_id)(runTrials_idx) = choice_hE_tmp;
                % extract level of monetary incentive proposed for high
                % effort option for each trial
                RP_var = strcmp(choiceOptions_tmp.R_or_P, 'R');
                RP_var = double(RP_var);
                RP_var(strcmp(choiceOptions_tmp.R_or_P, 'P')) = -1;
                money_level_left_tmp    = choiceOptions_tmp.monetary_amount.left.*RP_var;
                money_level_right_tmp   = choiceOptions_tmp.monetary_amount.right.*RP_var;
                inc_perTrial_tmp.(task_id)(runTrials_idx) = money_level_left_tmp.*(defaultSide_tmp == 1) +...
                    money_level_right_tmp.*(defaultSide_tmp == -1);
                % extract level of high effort option for each trial
                hE_eff_proposed_tmp = choiceOptions_tmp.E.left.*(defaultSide_tmp == 1) +...
                    choiceOptions_tmp.E.right.*(defaultSide_tmp == -1);
                eff_perTrial_tmp.(task_id)(runTrials_idx) = hE_eff_proposed_tmp;
                
            end % run to include?
        end % run loop
        
        % average the RT within each subject
        mean_RT.(task_id)(iS) = mean(raw_z_log(RT_perTrial_tmp.(task_id)),1,'omitnan');
        median_RT.(task_id)(iS) = median(raw_z_log(RT_perTrial_tmp.(task_id)),1,'omitnan');
        sd_RT.(task_id)(iS) = std(raw_z_log(RT_perTrial_tmp.(task_id)),[],1,'omitnan');
    end % physical/mental loop
end % subject loop

%% test link between parameters and RT
for iPrm = 1:nPrm
    prm_nm = parameters{iPrm};
    % filter subjects for which Arthur's model hasn't been performed yet
    good_subs = ~isnan(prm.(prm_nm));
    for iTask = 1:nTasks
        task_id = tasks{iTask};
        
        % MEAN RT test
        [beta_meanRT.(task_id).(prm_nm), ~,stats_meanRT.(task_id).(prm_nm)] = glmfit(prm.(prm_nm)(good_subs), mean_RT.(task_id)(good_subs), 'normal');
        meanRT_f_prm.(task_id).(prm_nm) = glmval(beta_meanRT.(task_id).(prm_nm), prm.(prm_nm)(good_subs), 'identity');
        pval.mean.(task_id).(prm_nm) = stats_meanRT.(task_id).(prm_nm).p(2);
        betas.mean.(task_id).(prm_nm) = beta_meanRT.(task_id).(prm_nm)(2);
        
        % MEDIAN RT test
        [beta_medianRT.(task_id).(prm_nm), ~,stats_medianRT.(task_id).(prm_nm)] = glmfit(prm.(prm_nm)(good_subs), median_RT.(task_id)(good_subs), 'normal');
        medianRT_f_prm.(task_id).(prm_nm) = glmval(beta_medianRT.(task_id).(prm_nm), prm.(prm_nm)(good_subs), 'identity');
        pval.median.(task_id).(prm_nm) = stats_medianRT.(task_id).(prm_nm).p(2);
        betas.median.(task_id).(prm_nm) = beta_medianRT.(task_id).(prm_nm)(2);

        % SD RT test
        [beta_sd_RT.(task_id).(prm_nm), ~,stats_sd_RT.(task_id).(prm_nm)] = glmfit(prm.(prm_nm)(good_subs), sd_RT.(task_id)(good_subs), 'normal');
        sd_RT_f_prm.(task_id).(prm_nm) = glmval(beta_sd_RT.(task_id).(prm_nm), prm.(prm_nm)(good_subs), 'identity');
        pval.sd.(task_id).(prm_nm) = stats_sd_RT.(task_id).(prm_nm).p(2);
        betas.sd.(task_id).(prm_nm) = beta_sd_RT.(task_id).(prm_nm)(2);

        %% figure display
        if figDisp == true
            pSize = 50;
            lWidth = 3;
            
            % look at mean, median and std RT = f(parameter)
            fig;
            jPlot = 0;
            
            % mean RT
            jPlot = jPlot + 1;
            subplot(1,3,jPlot);
            hold on;
            scatter(prm.(prm_nm)(good_subs), mean_RT.(task_id)(good_subs));
            plot(prm.(prm_nm)(good_subs), meanRT_f_prm.(task_id).(prm_nm),...
                'LineWidth',lWidth);
            xlabel(prm_nm);
            switch RT_type
                case 'raw'
                    ylabel([task_id,' - mean RT (s)']);
                case 'log'
                    ylabel([task_id,' - mean(log(RT))']);
                otherwise
                    error(['legend not ready yet for ',RT_type,' RT']);
            end
            legend_size(pSize);
            
            % median RT
            jPlot = jPlot + 1;
            subplot(1,3,jPlot);
            hold on;
            scatter(prm.(prm_nm)(good_subs), median_RT.(task_id)(good_subs));
            plot(prm.(prm_nm)(good_subs), medianRT_f_prm.(task_id).(prm_nm),...
                'LineWidth',lWidth);
            xlabel(prm_nm);
            switch RT_type
                case 'raw'
                    ylabel([task_id,' - median RT (s)']);
                case 'log'
                    ylabel([task_id,' - median(log(RT))']);
                otherwise
                    error(['legend not ready yet for ',RT_type,' RT']);
            end
            legend_size(pSize);

            % SD RT
            jPlot = jPlot + 1;
            subplot(1,3,jPlot);
            hold on;
            scatter(prm.(prm_nm)(good_subs), sd_RT.(task_id)(good_subs));
            plot(prm.(prm_nm)(good_subs), sd_RT_f_prm.(task_id).(prm_nm),...
                'LineWidth',lWidth);
            xlabel(prm_nm);
            switch RT_type
                case 'raw'
                    ylabel([task_id,' - sd RT (s)']);
                case 'log'
                    ylabel([task_id,' - sd(log(RT))']);
                otherwise
                    error(['legend not ready yet for ',RT_type,' RT']);
            end
            legend_size(pSize);
        end % figure display
        
    end % physical/mental loop
end % parameter loop

end % function