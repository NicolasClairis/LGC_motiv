function[betas, choices] = logitfit_choices_with_fMRI(computerRoot, study_nm, sub_nm, condition)
% [betas, choices] = logitfit_choices_with_fMRI(computerRoot, study_nm, sub_nm, condition)
% logitfit_choices_with_fMRI will look at a GLM explaining choices like
% logitfit_choices.m but also includes fMRI information.
%
% INPUTS
% computerRoot: pathway where data is
% 
% study_nm: study name
%
% sub_nm: subject number id 'XXX'
%
% condition: condition to use
%
% OUTPUTS
% betas: structure with betas
%
% choices: structure with bins for choices

%% working directories
% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = 'E:\';
%     computerRoot = LGCM_root_paths;
end

%% working directories
subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
    'CID',sub_nm, filesep, 'behavior', filesep];

%% condition
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end

%% define ROIs
[VS_ROI_infos] = load_VS_ROI();
[dmPFC_ROI_infos] = load_dmPFC_ROI();
% [aINS_ROI_infos] = load_aINS_ROI();

%% load BOLD for each ROI of interest
GLM = 94;
[VS_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, {sub_nm}, condition, GLM,...
    VS_ROI_infos);
[dmPFC_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, {sub_nm}, condition, GLM,...
    dmPFC_ROI_infos);
% [aINS_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
%     study_nm, {sub_nm}, condition, GLM,...
%     aINS_ROI_infos);

%% initialize variables of interest
nTrialsPerRun = 54;
% nTrialsPerRPConditionPerRun = nTrialsPerRun/2;
if strcmp(study_nm,'study1') && strcmp(sub_nm,'040')
    nRunsPerTask = 1;
else
    nRunsPerTask = 2;
end
nTrials = nTrialsPerRun*nRunsPerTask;
% nTrialsPerRPCond = nTrialsPerRPConditionPerRun*nRunsPerTask;

[choice_hE.Ep,...
    deltaE_level.Ep,...
    deltaR_money.Ep,...
    deltaP_money.Ep,...
    sumPrevAUC.Ep,...
    BOLD.VS.Ep,...
    BOLD.dmPFC.Ep,...
    BOLD.aINS.Ep,...
    choice_hE.Em,...
    deltaE_level.Em,...
    deltaR_money.Em,...
    deltaP_money.Em,...
    prevEfficacy.Em,...
    BOLD.VS.Em,...
    BOLD.dmPFC.Em,...
    BOLD.aINS.Em] = deal(NaN(nTrials,1));

%% extract runs
[runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
nRuns = length(runsStruct.tasks);
runs_Ep = strcmp(runsStruct.tasks,'Ep');
runs_Em = strcmp(runsStruct.tasks,'Em');

%% loop through physical and mental
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
            task_fullName = 'physical';
        case 2
            task_id = 'Em';
            task_fullName = 'mental';
    end
    
    jRun = 0;
    for iRun = 1:nRuns
        runToInclude = 0;
        switch task_id
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
        run_nm = num2str(runsStruct.runsToKeep(iRun));
        run_nm_bis = ['run',num2str(jRun)];
        
        if runToInclude == 1
            runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);
            % extract relevant data
            choice_hE.(task_id)(runTrials_idx) = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            [deltaE_level.(task_id)(runTrials_idx)] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            [deltaR_money.(task_id)(runTrials_idx)] = extract_deltaR_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            [deltaP_money.(task_id)(runTrials_idx)] = extract_deltaP_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            switch task_id
                case 'Ep'
                    sumPrevAUC.Ep(runTrials_idx) = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                case 'Em'
                    [~, ~, ~, ~,~,~,...
                        ~,~,...
                        efficacy_bis_with2first,...
                        ~] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
                    efficacy_tmp = efficacy_bis_with2first.allTrials;
                    prevEfficacy.Em(runTrials_idx) = [0,efficacy_tmp(2:end)];
            end % task filter to adapt temporal cost
            % extract fMRI data
            BOLD.VS.(task_id)(runTrials_idx) = VS_trial_b_trial.VS.(task_id).(run_nm_bis).choice(:,1);
            BOLD.dmPFC.(task_id)(runTrials_idx) = dmPFC_trial_b_trial.dmPFC.(task_id).(run_nm_bis).choice(:,1);
%             BOLD.aINS.(task_id)(runTrials_idx) = aINS_trial_b_trial.aINS.(task_id).(run_nm_bis).choice(:,1);
        end % run to include?
    end % run loop
    %% check number of runs
    if sum(runs_Ep) < 2 || sum(runs_Em) < 2
       warning('one run had to be removed. please adapt the script accordingly where needed.'); 
    end
    %% perform the fit
    switch task_id
        case 'Ep'
            %             xModel = [deltaR_money.Ep,...
            %                 (deltaR_money.Ep.*BOLD.VS.Ep),...
            %                 deltaP_money.Ep,...
            %                 (deltaP_money.Ep.*BOLD.aINS.Ep),...
            %                 deltaE_level.Ep,...
            %                 (deltaE_level.Ep.*BOLD.dmPFC.Ep),...
            %                 (deltaE_level.Ep.*sumPrevAUC.Ep),...
            %                 (deltaE_level.Ep.*sumPrevAUC.Ep.*BOLD.dmPFC.Ep)];
            xModel = [deltaR_money.Ep,...
                (deltaR_money.Ep.*BOLD.VS.Ep),...
                deltaP_money.Ep,...
                deltaE_level.Ep,...
                (deltaE_level.Ep.*BOLD.dmPFC.Ep),...
                (deltaE_level.Ep.*sumPrevAUC.Ep)];
        case 'Em'
            %             xModel = [deltaR_money.Em,...
            %                 (deltaR_money.Em.*BOLD.VS.Em),...
            %                 deltaP_money.Em,...
            %                 (deltaP_money.Em.*BOLD.aINS.Em),...
            %                 deltaE_level.Em,...
            %                 (deltaE_level.Em.*BOLD.dmPFC.Em),...
            %                 (deltaE_level.Em.*prevEfficacy.Em),...
            %                 (deltaE_level.Em.*prevEfficacy.Em.*BOLD.dmPFC.Em)];
            xModel = [deltaR_money.Em,...
                (deltaR_money.Em.*BOLD.VS.Em),...
                deltaP_money.Em,...
                deltaE_level.Em,...
                (deltaE_level.Em.*BOLD.dmPFC.Em),...
                (deltaE_level.Em.*prevEfficacy.Em)];
    end % adapt fit to the task
    
    % check whether there are NaNs
    nTrials = size(xModel,1);
    goodTrials = false(1,nTrials);
    for iTrial = 1:nTrials
        if ~isnan(sum(xModel(iTrial,:)))
            goodTrials(iTrial) = true;
        else
            goodTrials(iTrial) = false;
        end
    end
    
    %% only perform GLM if enough non-NaN data
    if sum(goodTrials) > 3
        % perform the model and extract the betas
        betaModel = glmfit(xModel(goodTrials,:), choice_hE.(task_id)(goodTrials),...
            'binomial','link','logit');
        % extract betas
%         betas.(task_id).bias = betaModel(1);
%         betas.(task_id).bR = betaModel(2);
%         betas.(task_id).b_VS_R = betaModel(3);
%         betas.(task_id).bP = betaModel(4);
%         betas.(task_id).b_aINS_P = betaModel(5);
%         betas.(task_id).bE = betaModel(6);
%         betas.(task_id).b_dmPFC_E = betaModel(7);
%         betas.(task_id).bF = betaModel(8);
%         betas.(task_id).b_dmPFC_F = betaModel(9);
         betas.(task_id).bias = betaModel(1);
        betas.(task_id).bR = betaModel(2);
        betas.(task_id).b_VS_R = betaModel(3);
        betas.(task_id).bP = betaModel(4);
        betas.(task_id).bE = betaModel(5);
        betas.(task_id).b_dmPFC_E = betaModel(6);
        betas.(task_id).bF = betaModel(7);
        % remove other variables
        [betas.(task_id).b_aINS_P, betas.(task_id).b_dmPFC_F] = deal(NaN);
        
        % extract fitted choices
        modelFit.(task_id) = glmval(betaModel, xModel, 'logit');
        confidence.(task_id) = (modelFit.(task_id) - 0.5).^2;
    else % not possible to apply the model => put NaN everywhere
        [betas.(task_id).bias,...
            betas.(task_id).bR,...
            betas.(task_id).b_VS_R,...
            betas.(task_id).bP,...
            betas.(task_id).b_aINS_P,...
            betas.(task_id).bE,...
            betas.(task_id).b_dmPFC_E,...
            betas.(task_id).bF,...
            betas.(task_id).b_dmPFC_F] = deal(NaN);
        [modelFit, confidence] = deal(NaN(1,nTrials));
    end
end % physical/mental

% extract output
choices.actualChoices = choice_hE;
choices.choicesFitted = modelFit;
choices.confidenceFitted = confidence;
end % function