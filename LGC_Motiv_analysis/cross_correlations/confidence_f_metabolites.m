%% script to check whether there is a global difference in confidence 
% ratings and in inferred confidence based on the level of brain
% metabolites.


%% working directory
[computerRoot] = LGCM_root_paths();
study_nm = 'study1';
studyFolder = fullfile(computerRoot, study_nm);

%% define all subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection('study1', condition);

%% define metabolite and ROI you want to focus on and extract subjects accordingly
[low_met_subs, high_met_subs, metabolite_nm] = medSplit_metabolites(study_nm, subject_id);

%% define behavioral model to use
[mdlType, mdlN] = behavioral_model_selection();
mdl_nm = ['mdl_',mdlN];
%% main parameters
nTrialsPerRun = 54;
nRuns = 4;
nTotalTrials = nTrialsPerRun*nRuns;
[conf_rated, conf_inferred] = deal(NaN(nTotalTrials,NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    % results folder
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];
    % extract runs for current subject
    [subjectRuns] = runs_definition(study_nm, sub_nm, 'behavior');

    % extract confidence inferred and confidence rated
    switch mdlType
        case 'simple'
            [~, dataInferred] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                0, 'levels', 6, 6);
        case 'bayesian'
            gitResultsFolder = [fullfile('C:','Users','clairis','Desktop',...
                'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];
    end
    for iRun = 1:nRuns
        run_nm = ['run',num2str(iRun)];
        if ismember(iRun, subjectRuns.runsToKeep)
            run_task_id = subjectRuns.tasks{iRun};
            switch run_task_id
                case 'Ep'
                    task_fullName = 'physical';
                case 'Em'
                    task_fullName = 'mental';
            end
            run_trial_idx = (1:nTrialsPerRun) + nTrialsPerRun*(nRuns - 1);

            % extract inferred confidence
            switch mdlType
                case 'simple'
                    conf_inferred(run_trial_idx,iS) = dataInferred.confidenceFitted.(mdl_nm).(run_nm);
                case 'bayesian'
                    [~, ~, confidence_run] = extract_bayesian_mdl(gitResultsFolder, subBehaviorFolder,...
                        sub_nm, num2str(iRun), task_fullName, mdl_nm);
                    conf_inferred(run_trial_idx,iS) = confidence_run;
            end

            % extract rated confidence
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                '_task.mat']);
            switch run_task_id
                case 'Ep'
                    choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                case 'Em'
                    choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
            end
            conf_tmp = choice_LR_tmp;
            conf_tmp(abs(choice_LR_tmp) == 2) = 1;
            conf_tmp(abs(choice_LR_tmp) == 1) = 0;
            conf_tmp(abs(choice_LR_tmp) == 0) = NaN;
            conf_rated(run_trial_idx, iS) = conf_tmp;
        end % filter good runs
    end % run loop
end % subject loop

%% take one value per subject
confInferredPerSub_bis = mean(conf_inferred,1,'omitnan');
confRatedPerSub_bis = mean(conf_rated,1,'omitnan');

%% average for each group
confInferred.low    = confInferredPerSub_bis(low_met_subs);
confInferred.high   = confInferredPerSub_bis(high_met_subs);
confRated.low       = confRatedPerSub_bis(low_met_subs);
confRated.high      = confRatedPerSub_bis(high_met_subs);
% compare both groups
[~,pval.confRated.low_vs_high] = ttest2(confRated.low, confRated.high);
[~,pval.confInferred.low_vs_high] = ttest2(confInferred.low, confInferred.high);
% average
[confInferred.avg.low, confInferred.sem.low] = mean_sem_sd(confInferred.low,2);
[confInferred.avg.high, confInferred.sem.high] = mean_sem_sd(confInferred.high,2);
[confRated.avg.low, confRated.sem.low] = mean_sem_sd(confRated.low,2);
[confRated.avg.high, confRated.sem.high] = mean_sem_sd(confRated.high,2);

%% display results
bWdth = 0.4;
lWdth = 3;
pSize = 30;

% inferred confidence
fig;
hold on;
bar(1-0.2, confInferred.avg.low, 'BarWidth',bWdth);
bar(1+0.2, confInferred.avg.high, 'BarWidth',bWdth);
errorbar(1-0.2, confInferred.avg.low, confInferred.sem.low, 'LineWidth',lWdth,'Color','k');
errorbar(1+0.2, confInferred.avg.high, confInferred.sem.high,'LineWidth',lWdth,'Color','k');
xticks([1-0.2, 1+0.2]);
xticklabels({['low ',metabolite_nm],['high ',metabolite_nm]});
ylabel('Conf inferred by the model');
legend_size(pSize);

% rated confidence
fig;
hold on;
bar(1-0.2, confRated.avg.low, 'BarWidth',bWdth);
bar(1+0.2, confRated.avg.high, 'BarWidth',bWdth);
errorbar(1-0.2, confRated.avg.low, confRated.sem.low, 'LineWidth',lWdth,'Color','k');
errorbar(1+0.2, confRated.avg.high, confRated.sem.high,'LineWidth',lWdth,'Color','k');
xticks([1-0.2, 1+0.2]);
xticklabels({['low ',metabolite_nm],['high ',metabolite_nm]});
ylabel('Conf rated');
legend_size(pSize);