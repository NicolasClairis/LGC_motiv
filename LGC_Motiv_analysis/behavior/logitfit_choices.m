function[betas, choices] = logitfit_choices(computerRoot, study_nm, sub_nm,...
    condition, figDisp, dispMoneyOrLevels, n_NV_bins, n_trialN_bins)
% [betas, choices] = logitfit_choices(computerRoot, study_nm, sub_nm,...
%       condition, figDisp, dispMoneyOrLevels, n_NV_bins, n_trialN_bins)
% logitfit_choices will perform a logistic regression on the choices
% performed by the participants
%
% INPUTS
% computerRoot: pathway where data is
% 
% study_nm: study name
%
% sub_nm: subject number id 'XXX'
%
% condition: condition indicating which runs to include or not
%
% figDisp: display individual figure (1) or not (0)
%
% dispMoneyOrLevels: display actual money ('money') or reward levels
% ('levels')
%
% n_NV_bins: number of bins for net value
%
% n_trialN_bins: number of bins for fatigue
%
% OUTPUTS
% betas: structure with betas
%
% choices: structure with bins for choices

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = 'E:\';
%     computerRoot = LGCM_root_paths;
end

%% working directories
subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
    'CID',sub_nm, filesep, 'behavior', filesep];

%% by default, display individual figure
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
    disp(['figDisp was not defined in the inputs so that by default ',...
        'figures are displayed for each individual.']);
end

%% by default, display monetary levels instead of actual monetary amounts
if ~exist('dispMoneyOrLevels','var') || isempty(dispMoneyOrLevels)
    dispMoneyOrLevels = 'levels';
end

%% extract runs
[runsStruct] = runs_definition(study_nm, sub_nm, condition);
nRuns = length(runsStruct.tasks);
runs_Ep = strcmp(runsStruct.tasks,'Ep');
runs_Em = strcmp(runsStruct.tasks,'Em');

%% define R/P/E values
money_levels = [-3, -2, -1, 1, 2, 3];
nMoneyLevels = length(money_levels);
[actualMoney_values.Ep,...
    actualMoney_values.Em]   = deal(NaN(1,nMoneyLevels));
E_levels = [1, 2, 3];
nELevels = length(E_levels);

if ~exist('n_NV_bins','var') || isempty(n_NV_bins)
    n_NV_bins = 6;
end

if ~exist('n_trialN_bins','var') || isempty(n_trialN_bins)
    n_trialN_bins = 6;
end

%% initialize variables of interest
nTrialsPerRun = 54;
% nTrialsPerRPConditionPerRun = nTrialsPerRun/2;
if strcmp(study_nm,'study1') && strcmp(sub_nm,'040')
    nRunsPerTask = 1;
else
    nRunsPerTask = 2;
end
nTrialsPerTask = nTrialsPerRun*nRunsPerTask;
% nTrialsPerRPCond = nTrialsPerRPConditionPerRun*nRunsPerTask;

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
    
    % initialize variables to store variables across all trials
    [choice_nonDef.(task_id),...
        trialN.(task_id),...
        R_nonDef_valuePerTrial.(task_id),...
        R_default_valuePerTrial.(task_id),...
        P_default_valuePerTrial.(task_id),...
        P_nonDef_valuePerTrial.(task_id),...
        money_default_valuePerTrial.(task_id),...
        money_nonDef_valuePerTrial.(task_id),...
        money_nonDef_levelPerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id),...
        E_default_levelPerTrial.(task_id)] = deal(NaN(nTrialsPerTask,1));
    nMdl = 4;
    for iRun = 1:nRuns
        kRun = runsStruct.runsToKeep(iRun);
        runToInclude = 0;
        switch task_id
            case 'Ep'
                if runs_Ep(iRun) == 1
                    runToInclude = 1;
                end
            case 'Em'
                if runs_Em(iRun) == 1
                    runToInclude = 1;
                end
        end
        run_nm = num2str(kRun);
        
        if runToInclude == 1
            switch kRun
                case {1,2}
                    jRun = 1;
                case {3,4}
                    jRun = 2;
            end
            run_nm_bis = ['run',num2str(jRun)];
            runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);
%             runTrialsRP_idx = (1:nTrialsPerRPCond) + nTrialsPerRPCond*(jRun-1);
            
            %% load data
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(kRun),'_',task_fullName,...
                '_task.mat']);
            choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
            
            %% load relevant data
            defaultSide_tmp = choiceOptions_tmp.default_LR;
            default_left_tmp = defaultSide_tmp == -1;
            default_right_tmp = defaultSide_tmp == 1;
            RP_var = strcmp(choiceOptions_tmp.R_or_P,'R');
            switch task_id
                case 'Ep'
                    choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                case 'Em'
                    choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
            end
            % remove confidence information from choice
            choice_LR_tmp(choice_LR_tmp == -2) = -1;
            choice_LR_tmp(choice_LR_tmp == 2) = 1;
            % variable with left or right choice
            choice_left_tmp = choice_LR_tmp == -1;
            choice_right_tmp = choice_LR_tmp == 1;
            % extract money/reward/punishment data
            money_value_nonDef_tmp = (choiceOptions_tmp.monetary_amount.left.*default_right_tmp +...
                choiceOptions_tmp.monetary_amount.right.*default_left_tmp).*((RP_var == 1) - (RP_var == 0));
            money_value_default_tmp = (choiceOptions_tmp.monetary_amount.left.*default_left_tmp +...
                choiceOptions_tmp.monetary_amount.right.*default_right_tmp).*((RP_var == 1) - (RP_var == 0));
            money_chosen_value_tmp = (choiceOptions_tmp.monetary_amount.left.*choice_left_tmp +...
                choiceOptions_tmp.monetary_amount.right.*choice_right_tmp).*((RP_var == 1) - (RP_var == 0));
            R_value_nonDef_tmp = (choiceOptions_tmp.monetary_amount.left.*default_right_tmp +...
                choiceOptions_tmp.monetary_amount.right.*default_left_tmp).*(RP_var == 1);
            R_value_default_tmp = (choiceOptions_tmp.monetary_amount.left.*default_left_tmp +...
                choiceOptions_tmp.monetary_amount.right.*default_right_tmp).*(RP_var == 1);
            R_chosen_value_tmp = (choiceOptions_tmp.monetary_amount.left.*choice_left_tmp +...
                choiceOptions_tmp.monetary_amount.right.*choice_right_tmp).*(RP_var == 1);
            P_value_nonDef_tmp = -(choiceOptions_tmp.monetary_amount.left.*default_right_tmp +...
                choiceOptions_tmp.monetary_amount.right.*default_left_tmp).*(RP_var == 0);
            P_value_default_tmp = -(choiceOptions_tmp.monetary_amount.left.*default_left_tmp +...
                choiceOptions_tmp.monetary_amount.right.*default_right_tmp).*(RP_var == 0);
            P_chosen_value_tmp = -(choiceOptions_tmp.monetary_amount.left.*choice_left_tmp +...
                choiceOptions_tmp.monetary_amount.right.*choice_right_tmp).*(RP_var == 0);
            % extract money levels
            money_level_nonDef_tmp = (choiceOptions_tmp.R.left.*default_right_tmp +...
                choiceOptions_tmp.R.right.*default_left_tmp).*((RP_var == 1) - (RP_var == 0));
            % extract effort levels
            E_nonDef_tmp = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            E_default_tmp = choiceOptions_tmp.E.left.*default_left_tmp +...
                choiceOptions_tmp.E.right.*default_right_tmp;
            E_chosen_tmp = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            
            
            R_money_tmp = getfield(load([subBehaviorFolder,...
                'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
                '_task_messyAllStuff.mat'],'R_money'),'R_money');
            jMoney = 0;
            for iMoney = money_levels
                jMoney = jMoney + 1;
                if iMoney < 0
                    actualMoney_values.(task_id)(jMoney) = -R_money_tmp.(['P_',num2str(-iMoney)]);
                elseif iMoney > 0
                    actualMoney_values.(task_id)(jMoney) = R_money_tmp.(['R_',num2str(iMoney)]);
                end
            end
            
            % extract relevant data
            choice_nonDef.(task_id)(runTrials_idx) = choice_LR_tmp == -defaultSide_tmp;
            trialN.(task_id)(runTrials_idx) = 1:nTrialsPerRun;
            R_nonDef_valuePerTrial.(task_id)(runTrials_idx) = R_value_nonDef_tmp;
            R_default_valuePerTrial.(task_id)(runTrials_idx) = R_value_default_tmp;
            P_nonDef_valuePerTrial.(task_id)(runTrials_idx) = P_value_nonDef_tmp;
            P_default_valuePerTrial.(task_id)(runTrials_idx) = P_value_default_tmp;
            money_default_valuePerTrial.(task_id)(runTrials_idx) = money_value_default_tmp;
            money_nonDef_valuePerTrial.(task_id)(runTrials_idx) = money_value_nonDef_tmp;
            money_nonDef_levelPerTrial.(task_id)(runTrials_idx) = money_level_nonDef_tmp;
            E_nonDef_levelPerTrial.(task_id)(runTrials_idx) = E_nonDef_tmp;
            E_default_levelPerTrial.(task_id)(runTrials_idx) = E_default_tmp;
            % extract variables/run
            trialN_perRun.(task_id).(run_nm_bis) = 1:nTrialsPerRun;
            money_chosen_valuePerTrial.(task_id).(run_nm_bis) = money_chosen_value_tmp;
            R_chosen_valuePerTrial.(task_id).(run_nm_bis) = R_chosen_value_tmp;
            P_chosen_valuePerTrial.(task_id).(run_nm_bis) = P_chosen_value_tmp;
            E_chosen_levelPerTrial.(task_id).(run_nm_bis) = E_chosen_tmp;
            money_varOption_valuePerTrialPerRun.(task_id).(run_nm_bis) = money_value_nonDef_tmp;
            R_varOption_valuePerTrialPerRun.(task_id).(run_nm_bis) = R_value_nonDef_tmp;
            P_varOption_valuePerTrialPerRun.(task_id).(run_nm_bis) = P_value_nonDef_tmp;
            E_varOption_levelPerTrialPerRun.(task_id).(run_nm_bis) = E_nonDef_tmp;
        end % run to include?
    end % run loop
    
    %% extract trials to be included in the GLM (ie remove excluded runs with NaNs from the analysis)
    okTrials = ~isnan(choice_nonDef.(task_id));
    
    %% perform the fit
    % model 1: SV = kMoney*Money-kE*E
    xModel1 = [money_nonDef_valuePerTrial.(task_id)-money_default_valuePerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id)];
    % perform the model and extract the betas
    betamdl_1 = glmfit(xModel1(okTrials,:), choice_nonDef.(task_id)(okTrials),...
        'binomial','link','logit','Constant','off');
    betas.(task_id).mdl_1.kMoney = betamdl_1(1);
    betas.(task_id).mdl_1.kEffort = betamdl_1(2);
    % extract fitted choices
    modelFit.mdl_1.(task_id) = glmval(betamdl_1, xModel1, 'logit','Constant','off');
    % extract corresponding confidence and net value inferred by the model
    confidence.mdl_1.(task_id).allTrials = (modelFit.mdl_1.(task_id) - 0.5).^2;
    deltaNV.mdl_1.(task_id) = betamdl_1(1).*xModel1(:,1) +...
        betamdl_1(2).*xModel1(:,2);
    
    % model 2: SV=kR*R+kP*P-kE*E (split R/P factor)
    xModel2 = [R_nonDef_valuePerTrial.(task_id)-R_default_valuePerTrial.(task_id),...
        P_nonDef_valuePerTrial.(task_id)-P_default_valuePerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id)];
    % perform the model and extract the betas
    betamdl_2 = glmfit(xModel2(okTrials,:), choice_nonDef.(task_id)(okTrials),...
        'binomial','link','logit','Constant','off');
    betas.(task_id).mdl_2.kR = betamdl_2(1);
    betas.(task_id).mdl_2.kP = betamdl_2(2);
    betas.(task_id).mdl_2.kEffort = betamdl_2(3);
    % extract fitted choices
    modelFit.mdl_2.(task_id) = glmval(betamdl_2, xModel2, 'logit','Constant','off');
    % extract corresponding confidence and net value inferred by the model
    confidence.mdl_2.(task_id).allTrials = (modelFit.mdl_2.(task_id) - 0.5).^2;
    deltaNV.mdl_2.(task_id) = betamdl_2(1).*xModel2(:,1) +...
        betamdl_2(2).*xModel2(:,2) +...
        betamdl_2(3).*xModel2(:,3);
    
    % model 3: SV = kMoney*Money-kE*E+kF*E*(trial number-1) model with fatigue
    xModel3 = [money_nonDef_valuePerTrial.(task_id)-money_default_valuePerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id),...
        (E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id)).*trialN.(task_id)];
    % perform the model and extract the betas
    betamdl_3 = glmfit(xModel3(okTrials,:), choice_nonDef.(task_id)(okTrials),...
        'binomial','link','logit','Constant','off');
    betas.(task_id).mdl_3.kMoney = betamdl_3(1);
    betas.(task_id).mdl_3.kEffort = betamdl_3(2);
    betas.(task_id).mdl_3.kFatigue = betamdl_3(3);
    % extract fitted choices
    modelFit.mdl_3.(task_id) = glmval(betamdl_3, xModel3, 'logit','Constant','off');
    % extract corresponding confidence and net value inferred by the model
    confidence.mdl_3.(task_id).allTrials = (modelFit.mdl_3.(task_id) - 0.5).^2;
    deltaNV.mdl_3.(task_id) = betamdl_3(1).*xModel3(:,1) +...
        betamdl_3(2).*xModel3(:,2) +...
        betamdl_3(3).*xModel3(:,3);
    
    % model 4: SV = kR*R-kE*E-kP*P+kF*E*(trial number-1) model with fatigue
    xModel4 = [R_nonDef_valuePerTrial.(task_id)-R_default_valuePerTrial.(task_id),...
        P_nonDef_valuePerTrial.(task_id)-P_default_valuePerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id),...
        (E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id)).*trialN.(task_id)];
    % perform the model and extract the betas
    betamdl_4 = glmfit(xModel4(okTrials,:), choice_nonDef.(task_id)(okTrials),...
        'binomial','link','logit','Constant','off');
    betas.(task_id).mdl_4.kR = betamdl_4(1);
    betas.(task_id).mdl_4.kP = betamdl_4(2);
    betas.(task_id).mdl_4.kEffort = betamdl_4(3);
    betas.(task_id).mdl_4.kFatigue = betamdl_4(4);
    % extract fitted choices
    modelFit.mdl_4.(task_id) = glmval(betamdl_4, xModel4, 'logit','Constant','off');
    % extract corresponding confidence and net value inferred by the model
    confidence.mdl_4.(task_id).allTrials = (modelFit.mdl_4.(task_id) - 0.5).^2;
    deltaNV.mdl_4.(task_id) = betamdl_4(1).*xModel4(:,1) +...
        betamdl_4(2).*xModel4(:,2) +...
        betamdl_4(3).*xModel4(:,3) +...
        betamdl_4(4).*xModel4(:,4);
    
    % extract data per run for each model
    for iRun = 1:nRuns
        kRun = runsStruct.runsToKeep(iRun);
        runToInclude = 0;
        switch task_id
            case 'Ep'
                if runs_Ep(iRun) == 1
                    runToInclude = 1;
                end
            case 'Em'
                if runs_Em(iRun) == 1
                    runToInclude = 1;
                end
        end
        run_fullNm = ['run',num2str(kRun)]; % run number with one single index
        
        if runToInclude == 1
            switch kRun
                case {1,2}
                    jRun = 1;
                case {3,4}
                    jRun = 2;
            end
            run_nm_bis = ['run',num2str(jRun)]; % run number where index is independent for each task
            runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);
            % model 1
            confidence.mdl_1.(run_fullNm) = confidence.mdl_1.(task_id).allTrials(runTrials_idx);
            NV_chosen.(task_id).mdl_1.(run_nm_bis) = betas.(task_id).mdl_1.kMoney.*money_chosen_valuePerTrial.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_1.kEffort.*E_chosen_levelPerTrial.(task_id).(run_nm_bis);
            NV_varOption.(task_id).mdl_1.(run_nm_bis) = betas.(task_id).mdl_1.kMoney.*money_varOption_valuePerTrialPerRun.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_1.kEffort.*E_varOption_levelPerTrialPerRun.(task_id).(run_nm_bis);
            confidence.mdl_1.(task_id).(run_nm_bis) = confidence.mdl_1.(task_id).allTrials(runTrials_idx);
            pChoice_hE.(task_id).mdl_1.(run_nm_bis) = modelFit.mdl_1.(task_id)(runTrials_idx);
            
            % model 2
            confidence.mdl_2.(run_fullNm) = confidence.mdl_2.(task_id).allTrials(runTrials_idx);
            NV_chosen.(task_id).mdl_2.(run_nm_bis) = betas.(task_id).mdl_2.kR.*R_chosen_valuePerTrial.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_2.kP.*P_chosen_valuePerTrial.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_2.kEffort.*E_chosen_levelPerTrial.(task_id).(run_nm_bis);
            NV_varOption.(task_id).mdl_2.(run_nm_bis) = betas.(task_id).mdl_2.kR.*R_varOption_valuePerTrialPerRun.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_2.kP.*P_varOption_valuePerTrialPerRun.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_2.kEffort.*E_varOption_levelPerTrialPerRun.(task_id).(run_nm_bis);
            confidence.mdl_2.(task_id).(run_nm_bis) = confidence.mdl_2.(task_id).allTrials(runTrials_idx);
            pChoice_hE.(task_id).mdl_2.(run_nm_bis) = modelFit.mdl_2.(task_id)(runTrials_idx);
            
            % model 3
            confidence.mdl_3.(run_fullNm) = confidence.mdl_3.(task_id).allTrials(runTrials_idx);
            NV_chosen.(task_id).mdl_3.(run_nm_bis) = betas.(task_id).mdl_3.kMoney.*money_chosen_valuePerTrial.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_3.kEffort.*E_chosen_levelPerTrial.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_3.kFatigue.*E_chosen_levelPerTrial.(task_id).(run_nm_bis).*trialN_perRun.(task_id).(run_nm_bis);
            NV_varOption.(task_id).mdl_3.(run_nm_bis) = betas.(task_id).mdl_3.kMoney.*money_varOption_valuePerTrialPerRun.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_3.kEffort.*E_varOption_levelPerTrialPerRun.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_3.kFatigue.*E_varOption_levelPerTrialPerRun.(task_id).(run_nm_bis).*trialN_perRun.(task_id).(run_nm_bis);
            confidence.mdl_3.(task_id).(run_nm_bis) = confidence.mdl_3.(task_id).allTrials(runTrials_idx);
            pChoice_hE.(task_id).mdl_3.(run_nm_bis) = modelFit.mdl_3.(task_id)(runTrials_idx);
            
            % model 4
            confidence.mdl_4.(run_fullNm) = confidence.mdl_4.(task_id).allTrials(runTrials_idx);
            NV_chosen.(task_id).mdl_4.(run_nm_bis) = betas.(task_id).mdl_4.kR.*R_chosen_valuePerTrial.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_4.kP.*P_chosen_valuePerTrial.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_4.kEffort.*E_chosen_levelPerTrial.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_4.kFatigue.*E_chosen_levelPerTrial.(task_id).(run_nm_bis).*trialN_perRun.(task_id).(run_nm_bis);
            NV_varOption.(task_id).mdl_4.(run_nm_bis) = betas.(task_id).mdl_4.kR.*R_varOption_valuePerTrialPerRun.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_4.kP.*P_varOption_valuePerTrialPerRun.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_4.kEffort.*E_varOption_levelPerTrialPerRun.(task_id).(run_nm_bis) +...
                betas.(task_id).mdl_4.kFatigue.*E_varOption_levelPerTrialPerRun.(task_id).(run_nm_bis).*trialN_perRun.(task_id).(run_nm_bis);
            confidence.mdl_4.(task_id).(run_nm_bis) = confidence.mdl_4.(task_id).allTrials(runTrials_idx);
            pChoice_hE.(task_id).mdl_4.(run_nm_bis) = modelFit.mdl_4.(task_id)(runTrials_idx);
        end % run to include
    end % run loop
    
    %% extract choice for each money and effort level
    % initialize variables for the bins
    [choiceNonDef.perMoneyLevel.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl_1.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl_2.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl_3.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl_4.(task_id),...
        confidence.perMoneyLevel.mdl_1.(task_id),...
        confidence.perMoneyLevel.mdl_2.(task_id),...
        confidence.perMoneyLevel.mdl_3.(task_id),...
        confidence.perMoneyLevel.mdl_4.(task_id)] = deal(NaN(1,nMoneyLevels));
    [choiceNonDef.perEffortLevel.(task_id),...
        choiceNonDef.perEffortLevel.(task_id),...
        choiceNonDef.perEffortLevel.(task_id),...
        choiceNonDef.perEffortLevel.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl_1.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl_2.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl_3.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl_4.(task_id),...
        confidence.perEffortLevel.mdl_1.(task_id),...
        confidence.perEffortLevel.mdl_2.(task_id),...
        confidence.perEffortLevel.mdl_3.(task_id),...
        confidence.perEffortLevel.mdl_4.(task_id)] = deal(NaN(1,nELevels));
    [choiceNonDef.perNVLevel.mdl_1.(task_id),...
        choiceFitNonDef.perNVLevel.mdl_1.(task_id),...
        deltaNV_bins.mdl_1.(task_id),...
        confidence.perNVLevel.mdl_1.(task_id),...
        choiceNonDef.perNVLevel.mdl_2.(task_id),...
        choiceFitNonDef.perNVLevel.mdl_2.(task_id),...
        deltaNV_bins.mdl_2.(task_id),...
        confidence.perNVLevel.mdl_2.(task_id),...
        choiceNonDef.perNVLevel.mdl_3.(task_id),...
        choiceFitNonDef.perNVLevel.mdl_3.(task_id),...
        deltaNV_bins.mdl_3.(task_id),...
        confidence.perNVLevel.mdl_3.(task_id),...
        choiceNonDef.perNVLevel.mdl_4.(task_id),...
        choiceFitNonDef.perNVLevel.mdl_4.(task_id),...
        deltaNV_bins.mdl_4.(task_id),...
        confidence.perNVLevel.mdl_4.(task_id)] = deal(NaN(1,n_NV_bins));
    [choiceNonDef.perTrialN.(task_id),...
        trialN_bins.(task_id),...
        choiceFitNonDef.perTrialN.mdl_1.(task_id),...
        choiceFitNonDef.perTrialN.mdl_2.(task_id),...
        choiceFitNonDef.perTrialN.mdl_3.(task_id),...
        choiceFitNonDef.perTrialN.mdl_4.(task_id),...
        confidence.perTrialN.mdl_1.(task_id),...
        confidence.perTrialN.mdl_2.(task_id),...
        confidence.perTrialN.mdl_3.(task_id),...
        confidence.perTrialN.mdl_4.(task_id)] = deal(NaN(1,n_trialN_bins));
    % money level
    jM = 0;
    for iMoney = money_levels
        jM = jM + 1;
        choiceNonDef.perMoneyLevel.(task_id)(jM) = mean(choice_nonDef.(task_id)( money_nonDef_levelPerTrial.(task_id) == iMoney),1,'omitnan');
        % loop through models
        for iMdl = 1:nMdl
            mdl_nm = ['mdl_',num2str(iMdl)];
            choiceFitNonDef.perMoneyLevel.(mdl_nm).(task_id)(jM) = mean(modelFit.(mdl_nm).(task_id)( money_nonDef_levelPerTrial.(task_id) == iMoney),1,'omitnan');
            confidence.perMoneyLevel.(mdl_nm).(task_id)(jM) = mean(confidence.(mdl_nm).(task_id).allTrials( money_nonDef_levelPerTrial.(task_id) == iMoney),1,'omitnan');
        end % model loop
    end
    % effort level
    jE = 0;
    for iEffort = E_levels
        jE = jE + 1;
        choiceNonDef.perEffortLevel.(task_id)(jE) = mean(choice_nonDef.(task_id)( E_nonDef_levelPerTrial.(task_id) == iEffort),1,'omitnan');
        % loop through models
        for iMdl = 1:nMdl
            mdl_nm = ['mdl_',num2str(iMdl)];
            choiceFitNonDef.perEffortLevel.(mdl_nm).(task_id)(jE) = mean(modelFit.(mdl_nm).(task_id)( E_nonDef_levelPerTrial.(task_id) == iEffort),1,'omitnan');
            confidence.perEffortLevel.(mdl_nm).(task_id)(jE) = mean(confidence.(mdl_nm).(task_id).allTrials( E_nonDef_levelPerTrial.(task_id) == iEffort),1,'omitnan');
        end % model loop
    end % effort level
    
    % net value
    % loop through models
    for iMdl = 1:nMdl
        mdl_nm = ['mdl_',num2str(iMdl)];
        [choiceNonDef.perNVLevel.(mdl_nm).(task_id),...
            deltaNV_bins.(mdl_nm).(task_id)] = do_bin2(choice_nonDef.(task_id), deltaNV.(mdl_nm).(task_id), n_NV_bins, 0);
        [choiceFitNonDef.perNVLevel.(mdl_nm).(task_id),...
            deltaNV_bins.(mdl_nm).(task_id)] = do_bin2(modelFit.(mdl_nm).(task_id), deltaNV.(mdl_nm).(task_id), n_NV_bins, 0);
        [confidence.perNVLevel.(mdl_nm).(task_id),...
            deltaNV_bins.(mdl_nm).(task_id)] = do_bin2(confidence.(mdl_nm).(task_id).allTrials, deltaNV.(mdl_nm).(task_id), n_NV_bins, 0);
    end
    
    % fatigue
    choiceNonDef.perTrialN.(task_id) = do_bin2(choice_nonDef.(task_id), trialN.(task_id), n_trialN_bins, 0);
    % model loop
    for iMdl = 1:nMdl
        mdl_nm = ['mdl_',num2str(iMdl)];
        [choiceFitNonDef.perTrialN.(mdl_nm).(task_id),...
            trialN_bins.(task_id)] = do_bin2(modelFit.(mdl_nm).(task_id),...
            trialN.(task_id),...
            n_trialN_bins, 0);
        [confidence.perTrialN.(mdl_nm).(task_id),...
            trialN_bins.(task_id)] = do_bin2(confidence.(mdl_nm).(task_id).allTrials,...
            trialN.(task_id),...
            n_trialN_bins, 0);
    end % model loop
    %% figures
    if figDisp == 1
        pSize = 30;
        lWidth = 3;
        lWidth_borders = 1;
        %% loop through models
        for iMdl = 1:nMdl
            %% display choice = f(net value)
            fig;
            pointMdl = scatter(deltaNV_bins.(['mdl_',num2str(iMdl)]).(task_id),...
                choiceNonDef.perNVLevel.(['mdl_',num2str(iMdl)]).(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            line(xlim(),[1 1],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            lHdlMdl = plot(deltaNV_bins.(['mdl_',num2str(iMdl)]).(task_id),...
                choiceFitNonDef.perNVLevel.(['mdl_',num2str(iMdl)]).(task_id));
            lHdlMdl.LineStyle = '--';
            lHdlMdl.LineWidth = lWidth;
            lHdlMdl.Color = [0 0 0];
            ylim([-0.2 1.2]);
            xlabel([task_fullName,' net value (non-default - default) - model ',num2str(iMdl)']);
            ylabel('Choice non-default option (%)');
            legend_size(pSize);
            
            %% choice non-default = f(money levels)
            switch dispMoneyOrLevels
                case 'money'
                    money_or_levels = actualMoney_values.(task_id);
                case 'levels'
                    money_or_levels = money_levels;
            end
            fig;
            pointMdl = scatter(money_or_levels,...
                choiceNonDef.perMoneyLevel.(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            line(xlim(),[1 1],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            lHdlMdl = plot(money_or_levels,...
                choiceFitNonDef.perMoneyLevel.(['mdl_',num2str(iMdl)]).(task_id));
            lHdlMdl.LineStyle = '--';
            lHdlMdl.LineWidth = lWidth;
            lHdlMdl.Color = [0 0 0];
            ylim([-0.2 1.2]);
            switch dispMoneyOrLevels
                case 'money'
                    xlabel([task_fullName,' Money (€) - model ',num2str(iMdl)]);
                case 'levels'
                    xlabel([task_fullName,' Money level - model ',num2str(iMdl)]);
            end
            ylabel('Choice non-default option (%)');
            legend_size(pSize);
            
            %% choice non-default = f(effort levels)
            fig;
            pointMdl = scatter(E_levels,...
                choiceNonDef.perEffortLevel.(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            line(xlim(),[1 1],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            lHdlMdl = plot(E_levels,...
                choiceFitNonDef.perEffortLevel.(['mdl_',num2str(iMdl)]).(task_id));
            lHdlMdl.LineStyle = '--';
            lHdlMdl.LineWidth = lWidth;
            lHdlMdl.Color = [0 0 0];
            ylim([-0.2 1.2]);
            xlabel([task_fullName,' effort level - model ',num2str(iMdl)']);
            ylabel('Choice non-default option (%)');
            legend_size(pSize);
            
            %% choice non-default = f(trial number)
            fig;
            pointMdl = scatter(trialN_bins.(task_id),...
                choiceNonDef.perTrialN.(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            line(xlim(),[1 1],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            lHdlMdl = plot(trialN_bins.(task_id),...
                choiceFitNonDef.perTrialN.(['mdl_',num2str(iMdl)]).(task_id));
            lHdlMdl.LineStyle = '--';
            lHdlMdl.LineWidth = lWidth;
            lHdlMdl.Color = [0 0 0];
            ylim([-0.2 1.2]);
            xlabel([task_fullName,' trial number - model ',num2str(iMdl)']);
            ylabel('Choice non-default option (%)');
            legend_size(pSize);
        end % model loop
    end % figure display
    
end % physical/mental

%% extract output
choices.deltaNV = deltaNV;
choices.choiceNonDef = choiceNonDef;
choices.choiceFitNonDef = choiceFitNonDef;
choices.choicesFitted = modelFit;
choices.confidenceFitted_highE = confidence;
choices.actualMoney_values = actualMoney_values;
choices.NV_bins = deltaNV_bins;
choices.trialN_bins = trialN_bins;
choices.NV_chosen = NV_chosen;
choices.NV_varOption = NV_varOption;
choices.pChoice_hE = pChoice_hE;

end % function