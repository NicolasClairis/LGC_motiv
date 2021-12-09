function[betas, choices] = logitfit_choices(computerRoot, study_nm, sub_nm, figDisp)
% [betas, choices] = logitfit_choices(computerRoot, study_nm, sub_nm, figDisp)
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

%% extract runs
[runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
nRuns = length(runsStruct.tasks);
runs_Ep = strcmp(runsStruct.tasks,'Ep');
runs_Em = strcmp(runsStruct.tasks,'Em');

%% define R/P/E values
P_levels = [-3, -2, -1];
nPlevels = length(P_levels);
R_levels = [1, 2, 3];
nRlevels = length(R_levels);
money_levels = [-3, -2, -1, 1, 2, 3];
nMoneyLevels = length(money_levels);
E_levels = [1, 2, 3];
nELevels = length(E_levels);

%% initialize variables of interest
nTrialsPerRun = 54;
nTrialsPerRPConditionPerRun = nTrialsPerRun/2;
nRunsPerTask = 2;
nTrials = nTrialsPerRun*nRunsPerTask;
nTrialsPerRPCond = nTrialsPerRPConditionPerRun*nRunsPerTask;

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
    
    [choice_nonDef.(task_id),...
        R_nonDef_valuePerTrial.(task_id),...
        R_default_valuePerTrial.(task_id),...
        P_default_valuePerTrial.(task_id),...
        P_nonDef_valuePerTrial.(task_id),...
        money_default_valuePerTrial.(task_id),...
        money_nonDef_valuePerTrial.(task_id),...
        money_nonDef_levelPerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id),...
        E_default_levelPerTrial.(task_id)] = deal(NaN(nTrials,1));
    [choice_nonDef_R.(task_id), choice_nonDef_P.(task_id),...
        R_nonDef_valuePerTrial_R.(task_id), P_nonDef_valuePerTrial_P.(task_id),...
        R_default_valuePerTrial_R.(task_id), P_default_valuePerTrial_P.(task_id),...
        E_nonDef_levelPerTrial_R.(task_id), E_nonDef_levelPerTrial_P.(task_id),...
        E_default_levelPerTrial_R.(task_id), E_default_levelPerTrial_P.(task_id)] = deal(NaN(nTrialsPerRPCond,1));
    [choiceNonDef.perMoneylevel.(task_id),...
        choiceFitNonDef.perMoneylevel.(task_id)] = deal(NaN(1,nMoneyLevels));
    [choiceNonDef.perEffortlevel.(task_id),...
        choiceFitNonDef.perEffortlevel.(task_id)] = deal(NaN(1,nELevels));
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
            runTrialsRP_idx = (1:nTrialsPerRPCond) + nTrialsPerRPCond*(jRun-1);
            
            %% load data
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                '_task.mat']);
            choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
            
            %% load relevant data
            switch task_id
                case 'Ep'
                    choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                case 'Em'
                    choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
            end
            % remove confidence information
            choice_LR_tmp(choice_LR_tmp == -2) = -1;
            choice_LR_tmp(choice_LR_tmp == 2) = 1;
            defaultSide_tmp = choiceOptions_tmp.default_LR;
            RP_var = strcmp(choiceOptions_tmp.R_or_P,'R');
            money_nonDef_tmp = choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == -1);
            money_default_tmp = choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == 1);
            R_nonDef_tmp = choiceOptions_tmp.monetary_amount.left.*(RP_var == 1).*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var == 1).*(defaultSide_tmp == -1);
            R_default_tmp = choiceOptions_tmp.monetary_amount.left.*(RP_var == 1).*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var == 1).*(defaultSide_tmp == 1);
            P_nonDef_tmp = choiceOptions_tmp.monetary_amount.left.*(RP_var == 0).*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var == 0).*(defaultSide_tmp == -1);
            P_default_tmp = choiceOptions_tmp.monetary_amount.left.*(RP_var == 0).*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var == 0).*(defaultSide_tmp == 1);
            E_nonDef_tmp = choiceOptions_tmp.E.left.*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.E.right.*(defaultSide_tmp == -1);
            E_default_tmp = choiceOptions_tmp.E.left.*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.E.right.*(defaultSide_tmp == 1);
            % extract money levels
            money_level_nonDef_tmp = choiceOptions_tmp.R.left.*(defaultSide_tmp == 1).*(RP_var == 1) +...
                choiceOptions_tmp.R.right.*(defaultSide_tmp == -1).*(RP_var == 1) +...
                -choiceOptions_tmp.R.left.*(defaultSide_tmp == 1).*(RP_var == 0) +...
                -choiceOptions_tmp.R.right.*(defaultSide_tmp == -1).*(RP_var == 0);
            
            % extract relevant data
            choice_nonDef.(task_id)(runTrials_idx) = choice_LR_tmp == defaultSide_tmp;
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
    betas.(task_id).model1.kb0 = betaMdl1(1);
    betas.(task_id).model1.kMoney = betaMdl1(2);
    betas.(task_id).model1.kEffort = betaMdl1(3);
    fitMdl1.(task_id) = glmval(betaMdl1, xModel1, 'logit');
    
    % model 2: SV=kR*R-kE*E-kP*P (split R/P factor)
    % model 3: SV = kMoney*Money*(1+kE*E²)
    % model 4: SV=(kR*R+kP*P)*(1+kE*E²)
    % model 5: SV = kR*R-kE*E-kP*P+kF*E*(trial number-1) model with fatigue
    
    %% extract choice for each money and effort level
    % money level
    jM = 0;
    for iMoney = money_levels
        jM = jM + 1;
        choiceNonDef.perMoneylevel.(task_id)(jM) = nanmean(choice_nonDef.(task_id)( money_nonDef_levelPerTrial.(task_id) == iMoney));
        choiceFitNonDef.perMoneylevel.(task_id)(jM) = nanmean(fitMdl1.(task_id)( money_nonDef_levelPerTrial.(task_id) == iMoney));
    end
    % effort level
    jE = 0;
    for iEffort = E_levels
        jE = jE + 1;
        choiceNonDef.perEffortlevel.(task_id)(jE) = nanmean(choice_nonDef.(task_id)( E_nonDef_levelPerTrial.(task_id) == iEffort));
        choiceFitNonDef.perEffortlevel.(task_id)(jE) = nanmean(fitMdl1.(task_id)( E_nonDef_levelPerTrial.(task_id) == iEffort));
    end
    
    %% figures
    if figDisp == 1
        pSize = 30;
        lWidth = 3;
        %% model 1
        %% display choice = f(net value)
        fig;
        pointMdl1 = scatter(money_levels, choiceNonDef.perMoneylevel.(task_id));
        pointMdl1.MarkerEdgeColor = [0 0 0];
        pointMdl1.MarkerFaceColor = [143 143 143]./255;
        pointMdl1.SizeData = 100;
        pointMdl1.LineWidth = lWidth;
        hold on;
        lHdlMdl1 = plot(money_levels, choiceFitNonDef.perMoneylevel.(task_id));
        lHdlMdl1.LineStyle = '--';
        lHdlMdl1.LineWidth = lWidth;
        lHdlMdl1.Color = [0 0 0];
        ylim([-0.2 1.2]);
        xlabel('Money');
        ylabel('Choice non-default option');
        legend_size(pSize);
    end
    
end % physical/mental

end % function