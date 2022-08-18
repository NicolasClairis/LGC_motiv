function[betas, pval, betas_grp, pval_grp] = RT_GLM(figDisp)
% [betas, pval, betas_grp, pval_grp] = RT_GLM()
%RT_GLM will perform a GLM on the reaction times
% 
% INPUTS
% figDisp: display figure
%
% OUTPUTS
% betas: structure with betas corresponding to each regressor for each
% subject
%
% pval: structure with corresponding p.value for each beta for each subject
%
% betas_grp: structure with mean, SEM and SD across subjects for each
% regressor impact on RT
%
% pval_grp: structure with p.value indicating how significant each beta is
% at the group level
%
% See also which_RT_GLM

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% study names
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% by default, display group figures
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
    disp(['figDispGroup was not defined in the inputs so that by default ',...
        'group figures are displayed.']);
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
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% initialize variables of interest
GLM_str = inputdlg('Which GLM?');
GLM = str2double(GLM_str);
[GLMprm] = which_RT_GLM(GLM);
potentialRegressors = fieldnames(GLMprm.regs);
for iReg = 1:length(potentialRegressors)
    curr_reg_nm = potentialRegressors{iReg};
    isRegON = strcmp(GLMprm.regs.(curr_reg_nm), 'on');
    if isRegON == true
        switch curr_reg_nm
            case 'run_cstt'
                [betas.run1_cstt, betas.run2_cstt,...
                    betas.run3_cstt, betas.run4_cstt] = deal(NaN(1,NS));
            case 'task_cstt'
                [betas.Em_cstt, betas.Ep_cstt] = deal(NaN(1,NS));
            otherwise
                betas.(curr_reg_nm) = NaN(1,NS);
        end
    end
end % regressor loop
beta_names = fieldnames(betas);

nTrialsPerRun = 54;
nRuns = 4;
nTotalTrials = nTrialsPerRun*nRuns;
[RT_raw, RT, RT_raw_prevTrial, RT_fit,...
    run1_cstt, run2_cstt, run3_cstt, run4_cstt,...
    Ep_cstt, Em_cstt,...
    trialN, choice_highE, E_highE, Echosen,...
    R_or_P,...
    deltaMoney, deltaMoneyChosen,...
    conf, confRtg] = deal(NaN(nTotalTrials, NS));

nBins = 6;
[RT_f_conf, RTfit_f_conf, conf_f_conf] = deal(NaN(nBins, NS));

%% perform the correlation
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
    % extract model data if need to model confidence
    switch GLMprm.main.confMdlType
        case 'simple'
            [~, dataInferred] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                0, 'levels', 6, 6);
        case 'bayesian'
            error('bayesian model not ready yet');
    end
    
    % extract runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    okRuns = runsStruct.runsToKeep;
%     badRuns = runsStruct.runsToIgnore;
    taskNames = runsStruct.tasks;
    jRun = 0;
    for iRun = okRuns
        jRun = jRun + 1;
        task_nm_tmp = taskNames{jRun};
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(iRun-1);
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        
        % run constant
        switch iRun
            case 1
                run1_cstt(runTrials_idx, iS) = 1;
                run2_cstt(runTrials_idx, iS) = 0;
                run3_cstt(runTrials_idx, iS) = 0;
                run4_cstt(runTrials_idx, iS) = 0;
            case 2
                run1_cstt(runTrials_idx, iS) = 0;
                run2_cstt(runTrials_idx, iS) = 1;
                run3_cstt(runTrials_idx, iS) = 0;
                run4_cstt(runTrials_idx, iS) = 0;
            case 3
                run1_cstt(runTrials_idx, iS) = 0;
                run2_cstt(runTrials_idx, iS) = 0;
                run3_cstt(runTrials_idx, iS) = 1;
                run4_cstt(runTrials_idx, iS) = 0;
            case 4
                run1_cstt(runTrials_idx, iS) = 0;
                run2_cstt(runTrials_idx, iS) = 0;
                run3_cstt(runTrials_idx, iS) = 0;
                run4_cstt(runTrials_idx, iS) = 1;
        end
        
        % task constant
        switch task_nm_tmp
            case 'Em'
                Ep_cstt(runTrials_idx, iS) = 0;
                Em_cstt(runTrials_idx, iS) = 1;
            case 'Ep'
                Ep_cstt(runTrials_idx, iS) = 1;
                Em_cstt(runTrials_idx, iS) = 0;
        end
        %% load the data
        behaviorStruct_tmp = load([subBehaviorFolder,...
            'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
            '_task.mat']);
        choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
        switch task_nm_tmp
            case 'Em'
                choiceAndPerf_tmp = behaviorStruct_tmp.mentalE_perf;
            case 'Ep'
                choiceAndPerf_tmp = behaviorStruct_tmp.physicalPerf;
        end
        
        %% RT
        onsets_tmp = behaviorStruct_tmp.onsets;
        switch task_nm_tmp
            case 'Ep'
                onsets_tmp = behaviorStruct_tmp.physicalPerf.onsets;
                choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
            case 'Em'
                onsets_tmp = behaviorStruct_tmp.mentalE_perf.onsets;
                choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
        end
        RT_tmp = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
        RT_raw(runTrials_idx, iS) = RT_tmp;
        RT_raw_prevTrial(runTrials_idx, iS) = [0, RT_tmp(1:(end-1))];
        switch GLMprm.main.RT_format
            case {'zscorePerRun','zscorePerRunARuns'}
                RT(runTrials_idx, iS) = nanzscore(RT_tmp);
            case {'raw','zscoreARuns'}
                RT(runTrials_idx, iS) = RT_tmp;
            case {'log'}
                RT(runTrials_idx, iS) = log(RT_tmp);
        end
        
        %% extract confidence rating
        confRtg_tmp = NaN(1,length(choice_LR_tmp));
        confRtg_tmp(abs(choice_LR_tmp) == 2) = 1;
        confRtg_tmp(abs(choice_LR_tmp) == 1) = 0;
        confRtg(runTrials_idx, iS) = confRtg_tmp;
        
        %% extract confidence inferred by the model
        switch GLMprm.main.confMdlType
            case 'simple'
                run_nm = ['run',num2str(iRun)];
                conf(runTrials_idx, iS) = dataInferred.confidenceFitted.(['mdl_',GLMprm.main.confMdlN]).(run_nm);
            case 'bayesian'
                error('bayesian model not ready yet');
        end
        
        %% default side
        defaultSide_tmp = choiceOptions_tmp.default_LR;
        
        %% choice
        choice_LR_tmp = choiceAndPerf_tmp.choice;
        choice_highE_tmp = NaN(1,length(choice_LR_tmp));
        choice_highE_tmp(choice_LR_tmp == -defaultSide_tmp) = 1;
        choice_highE_tmp(choice_LR_tmp == defaultSide_tmp) = 0;
        choice_highE(runTrials_idx, iS) = choice_highE_tmp;
        
        %% trial number
        trialN(runTrials_idx, iS) = 1:nTrialsPerRun;
        
        %% extract R or P
        RP_var_tmp = strcmp(choiceOptions_tmp.R_or_P,'R');
        R_or_P(runTrials_idx, iS) = RP_var_tmp;
        
        %% money
        money_nonDef_tmp = (choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == 1) +...
            choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == -1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
        money_default_tmp = (choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == -1) +...
            choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == 1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
        deltaMoney(runTrials_idx, iS) = money_nonDef_tmp - money_default_tmp;
        deltaMoneyChosen(runTrials_idx, iS) = (money_nonDef_tmp - money_default_tmp).*(choice_highE_tmp == 1) +...
            -(money_nonDef_tmp - money_default_tmp).*(choice_highE_tmp == 0);
        
        %% effort
        E_highE(runTrials_idx, iS) = choiceOptions_tmp.E.left.*(defaultSide_tmp == 1) +...
            choiceOptions_tmp.E.right.*(defaultSide_tmp == -1);
        Echosen(runTrials_idx, iS) = E_highE(runTrials_idx, iS).*(choice_highE_tmp == 1)' +...
            0.*(choice_highE_tmp == 0)';
        
    end % run loop
    
    %% preparing GLM
    potentialRegs.run1_cstt = run1_cstt(:,iS);
    potentialRegs.run2_cstt = run2_cstt(:,iS);
    potentialRegs.run3_cstt = run3_cstt(:,iS);
    potentialRegs.run4_cstt = run4_cstt(:,iS);
    potentialRegs.Em_cstt = Em_cstt(:,iS);
    potentialRegs.Ep_cstt = Ep_cstt(:,iS);
    potentialRegs.RT_raw_prevTrial = RT_raw_prevTrial(:,iS);
    potentialRegs.trialN = trialN(:,iS);
    potentialRegs.choice_highE = choice_highE(:,iS);
    potentialRegs.deltaEffort = E_highE(:,iS);
    potentialRegs.deltaEchosen = Echosen(:,iS);
    potentialRegs.R_or_P = R_or_P(:,iS);
    potentialRegs.deltaMoney = deltaMoney(:,iS);
    potentialRegs.deltaMoneyChosen = deltaMoneyChosen(:,iS);
    potentialRegs.conf = conf(:,iS);
    % extract regressors
    [x_regs, reg_names] = RT_GLM_regs(GLMprm, potentialRegs);
    
    % transform reaction times if necessary
    switch GLMprm.main.RT_format
        case {'zscoreARuns','zscorePerRunARuns'}
            RT(:, iS) = nanzscore(RT(:, iS));
    end
    
    %% performing GLM
    okTrials = (~isnan(RT(:, iS))).*(~isnan(RT_raw_prevTrial(:, iS))) == 1;
    [betas_tmp,~,stats_tmp] = glmfit(x_regs(okTrials, :), RT(okTrials, iS), 'normal',...
        'constant','off');
    RT_fit(okTrials,iS) = glmval(betas_tmp, x_regs(okTrials, :), 'identity','constant','off');
    
    % store betas
    for iReg = 1:length(reg_names)
        curr_reg_nm = reg_names{iReg};
        betas.(curr_reg_nm)(iS) = betas_tmp(iReg);
        pval.(curr_reg_nm)(iS) = stats_tmp.p(iReg);
    end % regressor loop
    
    %% extract binned data
    [RT_f_conf(:, iS), conf_f_conf(:, iS)] = do_bin2(RT(:,iS), conf(:,iS), nBins, 0);
    [RTfit_f_conf(:, iS), conf_f_conf(:, iS)] = do_bin2(RT_fit(:,iS), conf(:,iS), nBins, 0);
end % subject loop

%% average resulting betas
for iReg = 1:length(beta_names)
        curr_reg_nm = beta_names{iReg};
        [betas_grp.mean.(curr_reg_nm),...
            betas_grp.sem.(curr_reg_nm),...
            betas_grp.sd.(curr_reg_nm)] = mean_sem_sd(betas.(curr_reg_nm), 2);
        % test how significant betas are at the group level
        [~,pval_grp.(curr_reg_nm)] = ttest(betas.(curr_reg_nm));
end

% average also bins
[binsRT.mean.conf_f_conf, binsRT.sem.conf_f_conf] = mean_sem_sd(conf_f_conf, 2);
[binsRT.mean.RT_f_conf, binsRT.sem.RT_f_conf] = mean_sem_sd(RT_f_conf, 2);
[binsRT.mean.RTfit_f_conf, binsRT.sem.RTfit_f_conf] = mean_sem_sd(RTfit_f_conf, 2);

%% display figure
if figDisp == 1
    pSize = 50;
    lWidth = 3;
    
    % show RT = f(conf)
    fig;
    hold on;
    %     RT_hdl = errorbar(binsRT.mean.conf_f_conf,...
    %         binsRT.mean.RT_f_conf,...
    %         binsRT.mean.RT_f_conf - binsRT.sem.RT_f_conf,...
    %         binsRT.mean.RT_f_conf + binsRT.sem.RT_f_conf,...
    %         binsRT.mean.conf_f_conf - binsRT.sem.conf_f_conf,...
    %         binsRT.mean.conf_f_conf + binsRT.sem.conf_f_conf);
    %     RTfit_hdl = errorbar(binsRT.mean.conf_f_conf,...
    %         binsRT.mean.RTfit_f_conf,...
    %         binsRT.mean.RTfit_f_conf - binsRT.sem.RTfit_f_conf,...
    %         binsRT.mean.RTfit_f_conf + binsRT.sem.RTfit_f_conf,...
    %         binsRT.mean.conf_f_conf - binsRT.sem.conf_f_conf,...
    %         binsRT.mean.conf_f_conf + binsRT.sem.conf_f_conf);
    RT_hdl = errorbar(binsRT.mean.conf_f_conf,...
        binsRT.mean.RT_f_conf, binsRT.sem.RT_f_conf);
    RTfit_hdl = errorbar(binsRT.mean.conf_f_conf,...
        binsRT.mean.RTfit_f_conf, binsRT.sem.RTfit_f_conf);
    % improve visualization
    RT_hdl.LineWidth = lWidth;
    RTfit_hdl.LineWidth = lWidth;
    legend([RT_hdl, RTfit_hdl],{'RT','RT fitted'},'Location','Northeast');
    legend('boxoff');
    xlabel('Confidence');
    switch GLMprm.main.RT_format
        case 'raw'
            ylabel('RT (s)');
        case 'log'
            ylabel('log(RT)');
        otherwise
            ylabel('RT (u.a.)');
    end
    legend_size(pSize);
end % figure display

end % function