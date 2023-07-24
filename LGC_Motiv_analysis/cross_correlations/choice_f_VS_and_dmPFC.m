function[betas, pval] = choices_f_VS_and_dmPFC()
% [betas, pval] = choices_f_VS_and_dmPFC()
% choices_f_VS_and_dmPFC will test whether the selection of a high effort
% choice can be determined by the level of activity of the ventral striatum
% (VS) and of the dorsomedial prefrontal cortex (dmPFC)
%
% INPUTS
%
% OUTPUTS
% betas: structure with betas
%
% pval: structure with corresponding p.values
%
% See also logitfitchoices_with_fMRI.m and
% logitfit_choices_group_with_fMRI.m for a more complex implementation of
% this which includes rewards and efforts

%% working directories
% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = 'E:\';
%     computerRoot = LGCM_root_paths;
end

%% study name
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];
resultFolder_a = [studyBehaviorFolder,'results',filesep];
if ~exist(resultFolder_a,'dir')
    mkdir(resultFolder_a);
end
resultFolder = [resultFolder_a,'figures',filesep];
if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end

%% subject selection
if ~exist('subject_id','var') || isempty(subject_id)
    [condition] = subject_condition;
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end
% store subject list to know which beta corresponds to which subject
betas.subList = subject_id;

%% define ROIs
[VS_ROI_infos] = load_VS_ROI();
[dmPFC_ROI_infos] = load_dmPFC_ROI();

%% load BOLD for each ROI of interest
GLM = 94;
% [VS_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
%     study_nm, subject_id, condition, GLM,...
%     VS_ROI_infos);
% [dmPFC_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
%     study_nm, subject_id, condition, GLM,...
%     dmPFC_ROI_infos);
ROI_trialPerTrial = load(fullfile(resultFolder_a,'ROI','ROI_BOLDperTrial_GLM94_65subs_VS_dmPFC.mat'),...
    'VS_trial_b_trial','dmPFC_trial_b_trial');
VS_trial_b_trial = ROI_trialPerTrial.VS_trial_b_trial;
dmPFC_trial_b_trial = ROI_trialPerTrial.dmPFC_trial_b_trial;

%% initialize variables of interest
task_names = {'Ep','Em','EpEm'};
nTasks = length(task_names);
for iT = 1:nTasks
    task_nm = task_names{iT};
    [betas.linear_mdl.(task_nm).b0,...
        betas.linear_mdl.(task_nm).bVS,...
        betas.linear_mdl.(task_nm).bdmPFC,...
        betas.sigmo_mdl.(task_nm).b0,...
        betas.sigmo_mdl.(task_nm).bVS,...
        betas.sigmo_mdl.(task_nm).bdmPFC] = deal(NaN(1,NS));
end % task loop

% prepare variables trial per trial
nTrialsPerRun = 54;
nMaxRunsPerTask = 2;
nTrialsPerTask = nTrialsPerRun.*nMaxRunsPerTask;
nMaxRuns = 4;
nTrials = nTrialsPerRun.*nMaxRuns;
[choice_hE_perTrial.Ep,...
    VS_perTrial.Ep,...
    dmPFC_perTrial.Ep,...
    choice_hE_perTrial.Em,...
    VS_perTrial.Em,...
    dmPFC_perTrial.Em] = deal(NaN(nTrialsPerTask,NS));
[choice_hE_perTrial.EpEm,...
    VS_perTrial.EpEm,...
    dmPFC_perTrial.EpEm] = deal(NaN(nTrials,NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
        'CID',sub_nm, filesep, 'behavior', filesep];
    
    % extract info about relevant runs
    [runsStruct] = runs_definition(study_nm, sub_nm, condition);
    nRuns = length(runsStruct.tasks);
    
    % loop through runs
    jRun_Ep = 0;
    jRun_Em = 0;
    jRun_EpEm = 0;
    for iRun = 1:nRuns
        if ismember(iRun, runsStruct.runsToKeep)
            jRun_EpEm = jRun_EpEm + 1;
            run_nm = num2str(runsStruct.runsToKeep(jRun_EpEm));
            % physical or mental run?
            task_nm_tmp = runsStruct.tasks{jRun_EpEm};
            switch task_nm_tmp
                case 'Ep'
                    task_fullName = 'physical';
                    jRun_Ep = jRun_Ep + 1;
                    run_nm_bis = ['run',num2str(jRun_Ep)];
                case 'Em'
                    task_fullName = 'mental';
                    jRun_Em = jRun_Em + 1;
                    run_nm_bis = ['run',num2str(jRun_Em)];
            end
            
            
            % extract choice (high/low) and VS/dmPFC BOLD activity for each trial
            [choice_highE_tmp] = extract_choice_hE(subBehaviorFolder,...
                sub_nm, run_nm, task_fullName);
            VS_tmp = VS_trial_b_trial.VS.(task_nm_tmp).(run_nm_bis).choice(:,iS);
            dmPFC_tmp = dmPFC_trial_b_trial.dmPFC.(task_nm_tmp).(run_nm_bis).choice(:,iS);
            % distribute variables in the output
            switch task_nm_tmp
                case 'Ep'
                    trials_idx_perTask.Ep = (1:nTrialsPerRun) + nTrialsPerRun.*(jRun_Ep - 1);
                case 'Em'
                    trials_idx_perTask.Em = (1:nTrialsPerRun) + nTrialsPerRun.*(jRun_Em - 1);
            end
            choice_hE_perTrial.(task_nm_tmp)(trials_idx_perTask.(task_nm_tmp),iS) = choice_highE_tmp;
            VS_perTrial.(task_nm_tmp)(trials_idx_perTask.(task_nm_tmp),iS) = VS_tmp;
            dmPFC_perTrial.(task_nm_tmp)(trials_idx_perTask.(task_nm_tmp),iS) = dmPFC_tmp;
            % all trials
            trials_idx = (1:nTrialsPerRun) + nTrialsPerRun.*(jRun_EpEm - 1);
            choice_hE_perTrial.EpEm(trials_idx,iS) = choice_highE_tmp;
            VS_perTrial.EpEm(trials_idx,iS) = VS_tmp;
            dmPFC_perTrial.EpEm (trials_idx,iS) = dmPFC_tmp;
        end % run filter
    end % run loop
    
    %% GLM
    for iT = 1:nTasks
        task_nm_tmp = task_names{iT};
        goodTrials = ~isnan(choice_hE_perTrial.(task_nm_tmp)(:,iS).*VS_perTrial.(task_nm_tmp)(:,iS).*dmPFC_perTrial.(task_nm_tmp)(:,iS));
        % linear model
        [betas_tmp,~,stats_tmp] = glmfit([VS_perTrial.(task_nm_tmp)(goodTrials,iS),...
            dmPFC_perTrial.(task_nm_tmp)(goodTrials,iS)],...
            choice_hE_perTrial.(task_nm_tmp)(goodTrials,iS),...
            'normal');
        betas.linear_mdl.(task_nm_tmp).b0(iS) = betas_tmp(1);
        betas.linear_mdl.(task_nm_tmp).bVS(iS) = betas_tmp(2);
        betas.linear_mdl.(task_nm_tmp).bdmPFC(iS) = betas_tmp(3);
        
        % sigmoid model
        [betas_tmp,~,stats_tmp] = glmfit([VS_perTrial.(task_nm_tmp)(goodTrials,iS),...
            dmPFC_perTrial.(task_nm_tmp)(goodTrials,iS)],...
            choice_hE_perTrial.(task_nm_tmp)(goodTrials,iS),...
            'binomial', 'link', 'probit');
        betas.sigmo_mdl.(task_nm_tmp).b0(iS) = betas_tmp(1);
        betas.sigmo_mdl.(task_nm_tmp).bVS(iS) = betas_tmp(2);
        betas.sigmo_mdl.(task_nm_tmp).bdmPFC(iS) = betas_tmp(3);
        
        %% extract bins
        
    end % task loop
end % subject loop

%% test if betas are significantly different from zero
for iT = 1:nTasks
    task_nm_tmp = task_names{iT};
    % linear model
    [~,pval.linear_mdl.(task_nm_tmp).b0] = ttest(betas.linear_mdl.(task_nm_tmp).b0);
    [~,pval.linear_mdl.(task_nm_tmp).bVS] = ttest(betas.linear_mdl.(task_nm_tmp).bVS);
    [~,pval.linear_mdl.(task_nm_tmp).bdmPFC] = ttest(betas.linear_mdl.(task_nm_tmp).bdmPFC);
    % sigmoid model
    [~,pval.sigmo_mdl.(task_nm_tmp).b0] = ttest(betas.sigmo_mdl.(task_nm_tmp).b0);
    [~,pval.sigmo_mdl.(task_nm_tmp).bVS] = ttest(betas.sigmo_mdl.(task_nm_tmp).bVS);
    [~,pval.sigmo_mdl.(task_nm_tmp).bdmPFC] = ttest(betas.sigmo_mdl.(task_nm_tmp).bdmPFC);
end

%% figure

end % function