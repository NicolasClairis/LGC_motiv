function[betas, choices] = logitfit_choices(computerRoot, study_nm, sub_nm,...
    figDisp, dispMoneyOrLevels, n_NV_bins, n_trialN_bins)
% [betas, choices] = logitfit_choices(computerRoot, study_nm, sub_nm,...
%       figDisp, dispMoneyOrLevels, n_NV_bins, n_trialN_bins)
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
    computerRoot = LGCM_root_paths;
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
[runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
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
nRunsPerTask = 2;
nTrials = nTrialsPerRun*nRunsPerTask;
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
        E_default_levelPerTrial.(task_id)] = deal(NaN(nTrials,1));
    nMdl = 4;
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
        
        if runToInclude == 1
            runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);
%             runTrialsRP_idx = (1:nTrialsPerRPCond) + nTrialsPerRPCond*(jRun-1);
            
            %% load data
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                '_task.mat']);
            choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
            
            %% load relevant data
            defaultSide_tmp = choiceOptions_tmp.default_LR;
            RP_var = strcmp(choiceOptions_tmp.R_or_P,'R');
            money_nonDef_tmp = (choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == -1)).*((RP_var == 1) - (RP_var == 0));
            money_default_tmp = (choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == 1)).*((RP_var == 1) - (RP_var == 0));
            R_nonDef_tmp = choiceOptions_tmp.monetary_amount.left.*(RP_var == 1).*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var == 1).*(defaultSide_tmp == -1);
            R_default_tmp = choiceOptions_tmp.monetary_amount.left.*(RP_var == 1).*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var == 1).*(defaultSide_tmp == 1);
            P_nonDef_tmp = -(choiceOptions_tmp.monetary_amount.left.*(RP_var == 0).*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var == 0).*(defaultSide_tmp == -1));
            P_default_tmp = -(choiceOptions_tmp.monetary_amount.left.*(RP_var == 0).*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var == 0).*(defaultSide_tmp == 1));
            E_nonDef_tmp = choiceOptions_tmp.E.left.*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.E.right.*(defaultSide_tmp == -1);
            E_default_tmp = choiceOptions_tmp.E.left.*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.E.right.*(defaultSide_tmp == 1);
            % extract money levels
            money_level_nonDef_tmp = (choiceOptions_tmp.R.left.*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.R.right.*(defaultSide_tmp == -1)).*((RP_var == 1) - (RP_var == 0));
            switch task_id
                case 'Ep'
                    choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                case 'Em'
                    choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
            end
            % remove confidence information
            choice_LR_tmp(choice_LR_tmp == -2) = -1;
            choice_LR_tmp(choice_LR_tmp == 2) = 1;
            
            R_money_tmp = getfield(load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
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
            R_nonDef_valuePerTrial.(task_id)(runTrials_idx) = R_nonDef_tmp;
            R_default_valuePerTrial.(task_id)(runTrials_idx) = R_default_tmp;
            P_nonDef_valuePerTrial.(task_id)(runTrials_idx) = P_nonDef_tmp;
            P_default_valuePerTrial.(task_id)(runTrials_idx) = P_default_tmp;
            money_default_valuePerTrial.(task_id)(runTrials_idx) = money_default_tmp;
            money_nonDef_valuePerTrial.(task_id)(runTrials_idx) = money_nonDef_tmp;
            money_nonDef_levelPerTrial.(task_id)(runTrials_idx) = money_level_nonDef_tmp;
            E_nonDef_levelPerTrial.(task_id)(runTrials_idx) = E_nonDef_tmp;
            E_default_levelPerTrial.(task_id)(runTrials_idx) = E_default_tmp;
        end % run to include?
    end % run loop
    
    %% perform the fit
    % model 1: SV = kMoney*Money-kE*E
    xModel1 = [money_nonDef_valuePerTrial.(task_id)-money_default_valuePerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id)];
    betaMdl1 = glmfit(xModel1, choice_nonDef.(task_id),...
        'binomial','link','logit');
    betas.(task_id).mdl1.kb0 = betaMdl1(1);
    betas.(task_id).mdl1.kMoney = betaMdl1(2);
    betas.(task_id).mdl1.kEffort = betaMdl1(3);
    modelFit.mdl1.(task_id) = glmval(betaMdl1, xModel1, 'logit');
    deltaNV.mdl1.(task_id) = betaMdl1(1) +...
        betaMdl1(2).*xModel1(:,1) +...
        betaMdl1(3).*xModel1(:,2);
    
    % model 2: SV=kR*R+kP*P-kE*E (split R/P factor)
    xModel2 = [R_nonDef_valuePerTrial.(task_id)-R_default_valuePerTrial.(task_id),...
        P_nonDef_valuePerTrial.(task_id)-P_default_valuePerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id)];
    betaMdl2 = glmfit(xModel2, choice_nonDef.(task_id),...
        'binomial','link','logit');
    betas.(task_id).mdl2.kb0 = betaMdl2(1);
    betas.(task_id).mdl2.kR = betaMdl2(2);
    betas.(task_id).mdl2.kP = betaMdl2(3);
    betas.(task_id).mdl2.kEffort = betaMdl2(4);
    modelFit.mdl2.(task_id) = glmval(betaMdl2, xModel2, 'logit');
    deltaNV.mdl2.(task_id) = betaMdl2(1) +...
        betaMdl2(2).*xModel2(:,1) +...
        betaMdl2(3).*xModel2(:,2) +...
        betaMdl2(4).*xModel2(:,3);
    
    % model 3: SV = kMoney*Money-kE*E+kF*E*(trial number-1) model with fatigue
    xModel3 = [money_nonDef_valuePerTrial.(task_id)-money_default_valuePerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id),...
        (E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id)).*(trialN.(task_id) - 1)];
    betaMdl3 = glmfit(xModel3, choice_nonDef.(task_id),...
        'binomial','link','logit');
    betas.(task_id).mdl3.kb0 = betaMdl3(1);
    betas.(task_id).mdl3.kMoney = betaMdl3(2);
    betas.(task_id).mdl3.kEffort = betaMdl3(3);
    betas.(task_id).mdl3.kFatigue = betaMdl3(4);
    modelFit.mdl3.(task_id) = glmval(betaMdl3, xModel3, 'logit');
    deltaNV.mdl3.(task_id) = betaMdl3(1) +...
        betaMdl3(2).*xModel3(:,1) +...
        betaMdl3(3).*xModel3(:,2) +...
        betaMdl3(4).*xModel3(:,3);
    
    % model 4: SV = kR*R-kE*E-kP*P+kF*E*(trial number-1) model with fatigue
    xModel4 = [R_nonDef_valuePerTrial.(task_id)-R_default_valuePerTrial.(task_id),...
        P_nonDef_valuePerTrial.(task_id)-P_default_valuePerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id),...
        (E_nonDef_levelPerTrial.(task_id)-E_default_levelPerTrial.(task_id)).*(trialN.(task_id) - 1)];
    betaMdl4 = glmfit(xModel4, choice_nonDef.(task_id),...
        'binomial','link','logit');
    betas.(task_id).mdl4.kb0 = betaMdl4(1);
    betas.(task_id).mdl4.kR = betaMdl4(2);
    betas.(task_id).mdl4.kP = betaMdl4(3);
    betas.(task_id).mdl4.kEffort = betaMdl4(4);
    betas.(task_id).mdl4.kFatigue = betaMdl4(5);
    modelFit.mdl4.(task_id) = glmval(betaMdl4, xModel4, 'logit');
    deltaNV.mdl4.(task_id) = betaMdl4(1) +...
        betaMdl4(2).*xModel4(:,1) +...
        betaMdl4(3).*xModel4(:,2) +...
        betaMdl4(4).*xModel4(:,3) +...
        betaMdl4(5).*xModel4(:,4);
    
    %% extract choice for each money and effort level
    % initialize variables for the bins
    [choiceNonDef.perMoneyLevel.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl1.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl2.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl3.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl4.(task_id)] = deal(NaN(1,nMoneyLevels));
    [choiceNonDef.perEffortLevel.(task_id),...
        choiceNonDef.perEffortLevel.(task_id),...
        choiceNonDef.perEffortLevel.(task_id),...
        choiceNonDef.perEffortLevel.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl1.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl2.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl3.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl4.(task_id)] = deal(NaN(1,nELevels));
    [choiceNonDef.perNVLevel.mdl1.(task_id),...
        choiceFitNonDef.perNVLevel.mdl1.(task_id),...
        deltaNV_bins.mdl1.(task_id),...
        choiceNonDef.perNVLevel.mdl2.(task_id),...
        choiceFitNonDef.perNVLevel.mdl2.(task_id),...
        deltaNV_bins.mdl2.(task_id),...
        choiceNonDef.perNVLevel.mdl3.(task_id),...
        choiceFitNonDef.perNVLevel.mdl3.(task_id),...
        deltaNV_bins.mdl3.(task_id),...
        choiceNonDef.perNVLevel.mdl4.(task_id),...
        choiceFitNonDef.perNVLevel.mdl4.(task_id),...
        deltaNV_bins.mdl4.(task_id)] = deal(NaN(1,n_NV_bins));
    [choiceNonDef.perTrialN.mdl1.(task_id),...
        choiceFitNonDef.perTrialN.mdl1.(task_id),...
        trialN_bins.(task_id),...
        choiceNonDef.perTrialN.mdl2.(task_id),...
        choiceFitNonDef.perTrialN.mdl2.(task_id),...
        choiceNonDef.perTrialN.mdl3.(task_id),...
        choiceFitNonDef.perTrialN.mdl3.(task_id),...
        choiceNonDef.perTrialN.mdl4.(task_id),...
        choiceFitNonDef.perTrialN.mdl4.(task_id)] = deal(NaN(1,n_trialN_bins));
    % money level
    jM = 0;
    for iMoney = money_levels
        jM = jM + 1;
        choiceNonDef.perMoneyLevel.(task_id)(jM) = nanmean(choice_nonDef.(task_id)( money_nonDef_levelPerTrial.(task_id) == iMoney));
        % loop through models
        for iMdl = 1:nMdl
            mdl_nm = ['mdl',num2str(iMdl)];
            choiceFitNonDef.perMoneyLevel.(mdl_nm).(task_id)(jM) = nanmean(modelFit.(mdl_nm).(task_id)( money_nonDef_levelPerTrial.(task_id) == iMoney));
        end % model loop
    end
    % effort level
    jE = 0;
    for iEffort = E_levels
        jE = jE + 1;
        choiceNonDef.perEffortLevel.(task_id)(jE) = nanmean(choice_nonDef.(task_id)( E_nonDef_levelPerTrial.(task_id) == iEffort));
        % loop through models
        for iMdl = 1:nMdl
            mdl_nm = ['mdl',num2str(iMdl)];
            choiceFitNonDef.perEffortLevel.(mdl_nm).(task_id)(jE) = nanmean(modelFit.(mdl_nm).(task_id)( E_nonDef_levelPerTrial.(task_id) == iEffort));
        end % model loop
    end % effort level
    
    % net value
    % loop through models
    for iMdl = 1:nMdl
        mdl_nm = ['mdl',num2str(iMdl)];
        [choiceNonDef.perNVLevel.(mdl_nm).(task_id),...
            deltaNV_bins.(mdl_nm).(task_id)] = do_bin2(choice_nonDef.(task_id), deltaNV.(mdl_nm).(task_id), n_NV_bins, 0);
        [choiceFitNonDef.perNVLevel.(mdl_nm).(task_id),...
            deltaNV_bins.(mdl_nm).(task_id)] = do_bin2(modelFit.(mdl_nm).(task_id), deltaNV.(mdl_nm).(task_id), n_NV_bins, 0);
    end
    
    % fatigue
    choiceNonDef.perTrialN.(task_id) = do_bin2(choice_nonDef.(task_id), trialN.(task_id), n_trialN_bins, 0);
    % model loop
    for iMdl = 1:nMdl
        mdl_nm = ['mdl',num2str(iMdl)];
        [choiceFitNonDef.perTrialN.(mdl_nm).(task_id),...
            trialN_bins.(task_id)] = do_bin2(modelFit.(mdl_nm).(task_id),...
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
            pointMdl = scatter(deltaNV_bins.(['mdl',num2str(iMdl)]).(task_id),...
                choiceNonDef.perNVLevel.(['mdl',num2str(iMdl)]).(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            line(xlim(),[1 1],'LineWidth',lWidth_borders,'Color',[0 0 0]);
            lHdlMdl = plot(deltaNV_bins.(['mdl',num2str(iMdl)]).(task_id),...
                choiceFitNonDef.perNVLevel.(['mdl',num2str(iMdl)]).(task_id));
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
                choiceFitNonDef.perMoneyLevel.(['mdl',num2str(iMdl)]).(task_id));
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
                choiceFitNonDef.perEffortLevel.(['mdl',num2str(iMdl)]).(task_id));
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
                choiceFitNonDef.perTrialN.(['mdl',num2str(iMdl)]).(task_id));
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
choices.actualMoney_values = actualMoney_values;
choices.NV_bins = deltaNV_bins;
choices.trialN_bins = trialN_bins;

end % function