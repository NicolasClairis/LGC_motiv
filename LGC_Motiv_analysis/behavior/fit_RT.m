function[betas_RT, RTstruct, pval] = fit_RT(computerRoot, study_nm, sub_nm,...
    figDisp, dispMoneyOrLevels, n_NV_bins, n_trialN_bins, RT_type)
% [betas_RT, RTstruct, pval] = fit_RT(computerRoot, study_nm, sub_nm,...
%       figDisp, dispMoneyOrLevels, n_NV_bins, n_trialN_bins, RT_type)
% logitfit_RT will perform a linear regression on the reaction times (RT)
% of the participants
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
% RT_type:
% 'raw': raw reaction times (in seconds)
% 'zscored': zscore reaction times per run
% 'log': take log(RT)
%
% OUTPUTS
% betas_RT: structure with reaction time betas
%
% RTstruct: structure with bins for RT
%
% pval: p.value for the betas

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

%% by default look at raw RT
if ~exist('RT_type','var') || isempty(RT_type)
    RT_type = 'raw';
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

%% load choices fit
[betas_choices, choices] = logitfit_choices(computerRoot, study_nm, sub_nm,...
    0, dispMoneyOrLevels, n_NV_bins, n_trialN_bins);
% check number of models in the choices
choice_models = fieldnames(betas_choices.Ep);
nModels = length(choice_models);
for iMdl = 1:nModels
    NV_perTrial.Ep.(['mdl_',num2str(iMdl)]) = choices.deltaNV.(['mdl_',num2str(iMdl)]).Ep;
    NV_perTrial.Em.(['mdl_',num2str(iMdl)]) = choices.deltaNV.(['mdl_',num2str(iMdl)]).Em;
end % model loop

%% initialize variables of interest
nTrialsPerRun = 54;

%% loop through physical and mental
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
            task_fullName = 'physical';
            nRunsPerTask = sum(runs_Ep);
        case 2
            task_id = 'Em';
            task_fullName = 'mental';
            nRunsPerTask = sum(runs_Em);
    end
    nTrials = nTrialsPerRun*nRunsPerTask;
    
    % initialize variables to store variables across all trials
    [RT_perTrial.(task_id),...
        trialN.(task_id),...
        R_nonDef_valuePerTrial.(task_id),...
        R_default_valuePerTrial.(task_id),...
        P_default_valuePerTrial.(task_id),...
        P_nonDef_valuePerTrial.(task_id),...
        money_default_valuePerTrial.(task_id),...
        money_nonDef_valuePerTrial.(task_id),...
        money_nonDef_levelPerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id),...
        E_default_levelPerTrial.(task_id),...
        confidence_perTrial.(task_id),...
        R_or_P_trial.(task_id)] = deal(NaN(nTrials,1));
    
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
            
            %% load data
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                '_task.mat']);
            %% load relevant data
            choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
            onsets_tmp = behaviorStruct_tmp.onsets;
            switch task_id
                case 'Ep'
                    onsets_tmp = behaviorStruct_tmp.physicalPerf.onsets;
                    choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
                case 'Em'
                    onsets_tmp = behaviorStruct_tmp.mentalE_perf.onsets;
                    choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
            end
            RT_tmp = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
            switch RT_type
                case 'raw'
                    RT_perTrial.(task_id)(runTrials_idx) = RT_tmp;
                case 'zscored'
                    RT_perTrial.(task_id)(runTrials_idx) = nanzscore(RT_tmp);
                case 'log'
                    RT_perTrial.(task_id)(runTrials_idx) = log(RT_tmp);
            end
            
            % extract confidence level
            conf_tmp = abs(choice_LR_tmp) == 2;
            % extract value variables
            defaultSide_tmp = choiceOptions_tmp.default_LR;
            RP_var_tmp = strcmp(choiceOptions_tmp.R_or_P,'R');
            money_nonDef_tmp = (choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == -1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
            money_default_tmp = (choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == 1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
            R_nonDef_tmp = choiceOptions_tmp.monetary_amount.left.*(RP_var_tmp == 1).*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var_tmp == 1).*(defaultSide_tmp == -1);
            R_default_tmp = choiceOptions_tmp.monetary_amount.left.*(RP_var_tmp == 1).*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var_tmp == 1).*(defaultSide_tmp == 1);
            P_nonDef_tmp = -(choiceOptions_tmp.monetary_amount.left.*(RP_var_tmp == 0).*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var_tmp == 0).*(defaultSide_tmp == -1));
            P_default_tmp = -(choiceOptions_tmp.monetary_amount.left.*(RP_var_tmp == 0).*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.monetary_amount.right.*(RP_var_tmp == 0).*(defaultSide_tmp == 1));
            E_nonDef_tmp = choiceOptions_tmp.E.left.*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.E.right.*(defaultSide_tmp == -1);
            E_default_tmp = choiceOptions_tmp.E.left.*(defaultSide_tmp == -1) +...
                choiceOptions_tmp.E.right.*(defaultSide_tmp == 1);
            % extract money levels
            money_level_nonDef_tmp = (choiceOptions_tmp.R.left.*(defaultSide_tmp == 1) +...
                choiceOptions_tmp.R.right.*(defaultSide_tmp == -1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
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
            confidence_perTrial.(task_id)(runTrials_idx) = conf_tmp;
            R_or_P_trial.(task_id)(runTrials_idx) = RP_var_tmp;
        end % run filter
        
    end % run loop
    
    %% perform the fits
    % basic fit: RT = f(R, P, R/P, E, Conf)
    x.(task_id).mdl_0 = [R_nonDef_valuePerTrial.(task_id) - R_default_valuePerTrial.(task_id),...
        P_nonDef_valuePerTrial.(task_id) - P_default_valuePerTrial.(task_id),...
        R_or_P_trial.(task_id),...
        E_nonDef_levelPerTrial.(task_id) - E_default_levelPerTrial.(task_id),...
        confidence_perTrial.(task_id)];
    [betas_RT_mdl0_tmp, ~, stats_tmp] = glmfit(x.(task_id).mdl_0, RT_perTrial.(task_id), 'normal');
    betas_RT.(task_id).mdl_0.b0 = betas_RT_mdl0_tmp(1);
    betas_RT.(task_id).mdl_0.bR = betas_RT_mdl0_tmp(2);
    betas_RT.(task_id).mdl_0.bP = betas_RT_mdl0_tmp(3);
    betas_RT.(task_id).mdl_0.bRP = betas_RT_mdl0_tmp(4);
    betas_RT.(task_id).mdl_0.bE = betas_RT_mdl0_tmp(5);
    betas_RT.(task_id).mdl_0.bConf = betas_RT_mdl0_tmp(6);
    pval.(task_id).mdl_0 = stats_tmp.p;
    RTfit_perTrial.(task_id).mdl_0 = glmval(betas_RT_mdl0_tmp,...
        x.(task_id).mdl_0,'identity');
     
    % RT = f(Net Value, Confidence) for each choice model
    for iMdl = 1:nModels
        mdl_nm = ['mdl_',num2str(iMdl)];
        x.(task_id).(mdl_nm) = [NV_perTrial.(task_id).(mdl_nm),...
            confidence_perTrial.(task_id)];
        [betas_RT_mdl_tmp, ~, stats_tmp] = glmfit(x.(task_id).(mdl_nm), RT_perTrial.(task_id), 'normal');
        pval.(task_id).(mdl_nm) = stats_tmp.p;
        betas_RT.(task_id).(mdl_nm).b0 = betas_RT_mdl_tmp(1);
        betas_RT.(task_id).(mdl_nm).bNV = betas_RT_mdl_tmp(2);
        betas_RT.(task_id).(mdl_nm).bConf = betas_RT_mdl_tmp(3);
        RTfit_perTrial.(task_id).(mdl_nm) = glmval(betas_RT_mdl_tmp,...
            x.(task_id).(mdl_nm),'identity');
    end % model loop
    
    %% extract the bins for real data
    % money level
    RT_bins.perMoneyLevel.(task_id) = deal(NaN(1,nMoneyLevels));
    jM = 0;
    for iMoney = money_levels
        jM = jM + 1;
        RT_bins.perMoneyLevel.(task_id)(jM) = nanmean(RT_perTrial.(task_id)( money_nonDef_levelPerTrial.(task_id) == iMoney));
    end
    
    % pool all reward and all punishment trial
    RT_bins.RP.(task_id) = NaN(1,2);
    RT_bins.RP.(task_id)(1) = nanmean(RT_perTrial.(task_id)( R_or_P_trial.(task_id) == 1));
    RT_bins.RP.(task_id)(2) = nanmean(RT_perTrial.(task_id)( R_or_P_trial.(task_id) == 0));
    
    % effort level
    RT_bins.perEffortLevel.(task_id) = deal(NaN(1,nELevels));
    jE = 0;
    for iEffort = E_levels
        jE = jE + 1;
        RT_bins.perEffortLevel.(task_id)(jE) = nanmean(RT_perTrial.(task_id)( E_nonDef_levelPerTrial.(task_id) == iEffort));
    end % effort level
    
    % trial number
    [RT_bins.perTrialN.(task_id),...
        trialN_bins.(task_id)] = do_bin2(RT_perTrial.(task_id),...
        trialN.(task_id), n_trialN_bins, 0);
    
    % confidence
    RT_bins.perConfidenceLevel.(task_id) = deal(NaN(1, 2));
    RT_bins.perConfidenceLevel.(task_id)(1) = nanmean(RT_perTrial.(task_id)( confidence_perTrial.(task_id) == 0));
    RT_bins.perConfidenceLevel.(task_id)(2) = nanmean(RT_perTrial.(task_id)( confidence_perTrial.(task_id) == 1));
    
    %% extract the bins for data based on the models
    for iMdl = 0:nModels
        mdl_nm = ['mdl_',num2str(iMdl)];
        
        % money level
        RTfit_bins.perMoneyLevel.(mdl_nm).(task_id) = deal(NaN(1,nMoneyLevels));
        jM = 0;
        for iMoney = money_levels
            jM = jM + 1;
            RTfit_bins.perMoneyLevel.(mdl_nm).(task_id)(jM) = nanmean(RTfit_perTrial.(task_id).(mdl_nm)( money_nonDef_levelPerTrial.(task_id) == iMoney));
        end
        
        % reward/punishment
        RTfit_bins.RP.(mdl_nm).(task_id) = NaN(1,2);
        RTfit_bins.RP.(mdl_nm).(task_id)(1) = nanmean(RTfit_perTrial.(task_id).(mdl_nm)( R_or_P_trial.(task_id) == 1));
        RTfit_bins.RP.(mdl_nm).(task_id)(2) = nanmean(RTfit_perTrial.(task_id).(mdl_nm)( R_or_P_trial.(task_id) == 0));
    
        % effort level
        RTfit_bins.perEffortLevel.(mdl_nm).(task_id) = deal(NaN(1,nELevels));
        jE = 0;
        for iEffort = E_levels
            jE = jE + 1;
            RTfit_bins.perEffortLevel.(mdl_nm).(task_id)(jE) = nanmean(RTfit_perTrial.(task_id).(mdl_nm)( E_nonDef_levelPerTrial.(task_id) == iEffort));
        end % effort level
        
        % net value based on the models
        if iMdl > 0
            [RT_bins.perNVlevel.(mdl_nm).(task_id),...
                RTfit_bins.perNVlevel.(mdl_nm).(task_id),...
                NV_bins.perNVlevel.(mdl_nm).(task_id)] = deal(NaN(1,n_NV_bins));
            [RT_bins.perNVlevel.(mdl_nm).(task_id), ~] = do_bin2(RT_perTrial.(task_id),...
                NV_perTrial.(task_id).(mdl_nm), n_NV_bins, 0);
            [RTfit_bins.perNVlevel.(mdl_nm).(task_id),...
                NV_bins.perNVlevel.(mdl_nm).(task_id)] = do_bin2(RTfit_perTrial.(task_id).(mdl_nm),...
                NV_perTrial.(task_id).(mdl_nm), n_NV_bins, 0);
        end
        
        % confidence
        RTfit_bins.perConfidenceLevel.(mdl_nm).(task_id) = deal(NaN(1, 2));
        RTfit_bins.perConfidenceLevel.(mdl_nm).(task_id)(1) = nanmean(RTfit_perTrial.(task_id).(mdl_nm)( confidence_perTrial.(task_id) == 0));
        RTfit_bins.perConfidenceLevel.(mdl_nm).(task_id)(2) = nanmean(RTfit_perTrial.(task_id).(mdl_nm)( confidence_perTrial.(task_id) == 1));
        
        % trial number
        RTfit_bins.perTrialN.(mdl_nm).(task_id) = NaN(1,n_trialN_bins);
        [RTfit_bins.perTrialN.(mdl_nm).(task_id)] = do_bin2(RTfit_perTrial.(task_id).(mdl_nm),...
            trialN_bins.(task_id), n_trialN_bins);
    end % model loop
    
    %% figures
    if figDisp == 1
        pSize = 30;
        lWidth = 3;
        lWidth_borders = 1;
        
        %% loop through models (for the fit)
        for iMdl = 0:nModels
            mdl_nm = ['mdl_',num2str(iMdl)];
            
            %% RT = f(monetary amounts)
            switch dispMoneyOrLevels
                case 'money'
                    money_or_levels = actualMoney_values.(task_id);
                case 'levels'
                    money_or_levels = money_levels;
            end
            fig;
            pointMdl = scatter(money_or_levels,...
                RT_bins.perMoneyLevel.(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            lHdlMdl = plot(money_or_levels,...
                RTfit_bins.perMoneyLevel.(['mdl_',num2str(iMdl)]).(task_id));
            lHdlMdl.LineStyle = '--';
            lHdlMdl.LineWidth = lWidth;
            lHdlMdl.Color = [0 0 0];
            switch dispMoneyOrLevels
                case 'money'
                    xlabel([task_fullName,' Money (€) - model ',num2str(iMdl)]);
                case 'levels'
                    xlabel([task_fullName,' Money level - model ',num2str(iMdl)]);
            end
            y_legend(RT_type);
            legend_size(pSize);
            
            %% RT = f(R/P)
            fig;
            pointMdl = scatter(1:2,...
                RT_bins.RP.(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            lHdlMdl = plot(1:2,...
                RTfit_bins.RP.(['mdl_',num2str(iMdl)]).(task_id));
            lHdlMdl.LineStyle = '--';
            lHdlMdl.LineWidth = lWidth;
            lHdlMdl.Color = [0 0 0];
            xticks(1:2);
            xticklabels({'R','P'});
            xlim([0.5 2.5]);
            xlabel(['Model ',num2str(iMdl)]);
            y_legend(RT_type);
            legend_size(pSize);
            
            %% RT = f(E levels)
            fig;
            pointMdl = scatter(E_levels,...
                RT_bins.perEffortLevel.(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            lHdlMdl = plot(E_levels,...
                RTfit_bins.perEffortLevel.(['mdl_',num2str(iMdl)]).(task_id));
            lHdlMdl.LineStyle = '--';
            lHdlMdl.LineWidth = lWidth;
            lHdlMdl.Color = [0 0 0];
            xlabel([task_fullName,' effort level - model ',num2str(iMdl)']);
            y_legend(RT_type);
            legend_size(pSize);
            
            %% RT = f(NV)
            if iMdl > 0
                fig;
                pointMdl = scatter(NV_bins.perNVlevel.(mdl_nm).(task_id),...
                    RT_bins.perNVlevel.(mdl_nm).(task_id));
                pointMdl.MarkerEdgeColor = [0 0 0];
                pointMdl.MarkerFaceColor = [143 143 143]./255;
                pointMdl.SizeData = 100;
                pointMdl.LineWidth = lWidth;
                hold on;
                lHdlMdl = plot(NV_bins.perNVlevel.(mdl_nm).(task_id),...
                    RTfit_bins.perNVlevel.(['mdl_',num2str(iMdl)]).(task_id));
                lHdlMdl.LineStyle = '--';
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = [0 0 0];
                xlabel([task_fullName,' net value - model ',num2str(iMdl)']);
                y_legend(RT_type);
                legend_size(pSize);
            end
            
            %% RT = f(confidence)
            fig;
            pointMdl = scatter(1:2,...
                RT_bins.perConfidenceLevel.(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            lHdlMdl = plot(1:2,...
                RTfit_bins.perConfidenceLevel.(['mdl_',num2str(iMdl)]).(task_id));
            lHdlMdl.LineStyle = '--';
            lHdlMdl.LineWidth = lWidth;
            lHdlMdl.Color = [0 0 0];
            xticks(1:2);
            xticklabels({'low','high'});
            xlabel([task_fullName,' Confidence - model ',num2str(iMdl)']);
            xlim([0.5 2.5]);
            y_legend(RT_type);
            legend_size(pSize);
            
            %% RT = f(trial number)
            fig;
            pointMdl = scatter(trialN_bins.(task_id),...
                RT_bins.perTrialN.(task_id));
            pointMdl.MarkerEdgeColor = [0 0 0];
            pointMdl.MarkerFaceColor = [143 143 143]./255;
            pointMdl.SizeData = 100;
            pointMdl.LineWidth = lWidth;
            hold on;
            lHdlMdl = plot(trialN_bins.(task_id),...
                RTfit_bins.perTrialN.(['mdl_',num2str(iMdl)]).(task_id));
            lHdlMdl.LineStyle = '--';
            lHdlMdl.LineWidth = lWidth;
            lHdlMdl.Color = [0 0 0];
            xlabel([task_fullName,' trial number - model ',num2str(iMdl)']);
            xlim([0 nTrialsPerRun]);
            y_legend(RT_type);
            legend_size(pSize);
        end % model loop
        
    end % figure display

    %% extract output
    RTstruct.RT_bins = RT_bins;
    RTstruct.RTfit_bins = RTfit_bins;
    RTstruct.trialN_bins = trialN_bins;
    RTstruct.NV_bins = NV_bins;
end % physical/mental loop

end % function

function[] = y_legend(RT_type)
switch RT_type
    case 'raw'
        ylabel('RT (s)');
    case 'zscored'
        ylabel('z(RT) (a.u.)');
    case 'log'
        ylabel('log(RT) (a.u.)');
end
end