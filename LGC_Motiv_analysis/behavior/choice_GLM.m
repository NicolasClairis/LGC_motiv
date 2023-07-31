function[betas, pval, betas_grp, pval_grp, choice_fit_perRun] = choice_GLM(figDisp, computerRoot, study_nm, subject_id, GLM)
% [betas, pval, betas_grp, pval_grp, choice_fit_perRun] = choice_GLM(figDisp, computerRoot, study_nm, subject_id, GLM)
%choice_GLM will perform a GLM on the choices. The different with
%logitfit_choices is that choice_GLM pools all runs and tasks together,
%while logitfit_choices looks at each task independently.
% 
% INPUTS
% figDisp: display figure (yes by default if left empty)
%
% computerRoot: path to computer root (will be asked by default if left empty)
%
% study_nm: study name 'study1' or 'study2' ? (study1 by default if left empty)
%
% subject_id: list of subjects (will be asked by default if left empty)
%
% GLM: which choice GLM to use?
%
% OUTPUTS
% betas: structure with betas corresponding to each regressor for each
% subject
%
% pval: structure with corresponding p.value for each beta for each subject
%
% betas_grp: structure with mean, SEM and SD across subjects for each
% regressor impact on choices
%
% pval_grp: structure with p.value indicating how significant each beta is
% at the group level
%
% choice_fit_perRun: structure with 1 vector framed as nTrialsPerRun*(number of
% subjects) for each run allowing to extract and re-use the values of the
% fit in other scripts (see for ex. choice_RT_heatmap_f_dR_dE.m)
%
% See also which_choice_GLM

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
if ~exist('subject_id','var') || isempty(subject_id)
    condition = subject_condition;
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end

%% initialize variables of interest
if ~exist('GLM','var') || isempty(GLM) || GLM <=0
    GLM_str = inputdlg('Which choice GLM?');
    GLM = str2double(GLM_str);
end
[GLMprm] = which_choice_GLM(GLM);
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
[choice_hE, choice_hE_fit, NV_fit,...
    run1_cstt, run2_cstt, run3_cstt, run4_cstt,...
    Ep_cstt, Em_cstt,...
    choice_bias,...
    trialN,...
    RP, deltaR, deltaP,...
    deltaMoney, deltaMoneyChosen,...
    E_highE, Echosen, deltaEp, deltaEm,...
    physical_Fatigue, mental_Facilitation] = deal(NaN(nTotalTrials, NS));
[choice_fit_perRun.run1,...
    choice_fit_perRun.run2,...
    choice_fit_perRun.run3,...
    choice_fit_perRun.run4] = deal(NaN(nTrialsPerRun, NS));

nBins = 6;
[choice_f_NV, choice_fit_f_NV, NV_f_NV] = deal(NaN(nBins, NS));
%% perform the correlation
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder,...
        'CID',sub_nm, filesep, 'behavior',filesep];
    
    % extract runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    okRuns = runsStruct.runsToKeep;
%     badRuns = runsStruct.runsToIgnore;
    taskNames = runsStruct.tasks;
    jRun = 0;
    for iRun = okRuns
        run_nm = num2str(iRun);
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

        % choice bias
        choice_bias(runTrials_idx, iS) = 1;

        %% load the data
        behaviorStruct_tmp = load([subBehaviorFolder,...
            'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
            '_task.mat']);
        choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
        
        %% choice
        choice_highE_tmp = extract_choice_hE(subBehaviorFolder,sub_nm,run_nm,task_fullName);
        choice_hE(runTrials_idx, iS) = choice_highE_tmp;

        %% trial number
        trialN(runTrials_idx, iS) = 1:nTrialsPerRun;
        
        %% extract R or P
        [RP_var_tmp] = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        % binarize RP_trial (because equal to +1 for R and -1 for P)
        RP_var_tmp(RP_var_tmp == -1) = 0;
        RP(runTrials_idx, iS) = RP_var_tmp;
        
        %% extract deltaR and deltaP
        [deltaR(runTrials_idx, iS)] = extract_deltaR_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [deltaP(runTrials_idx, iS)] = extract_deltaP_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
%         [deltaR(runTrials_idx, iS)] = extract_hR_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
%         [deltaP(runTrials_idx, iS)] = extract_hP_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        
        %% money
        defaultSide_tmp = choiceOptions_tmp.default_LR;
        money_nonDef_tmp = (choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == 1) +...
            choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == -1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
        money_default_tmp = (choiceOptions_tmp.monetary_amount.left.*(defaultSide_tmp == -1) +...
            choiceOptions_tmp.monetary_amount.right.*(defaultSide_tmp == 1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
        deltaMoney(runTrials_idx, iS) = money_nonDef_tmp - money_default_tmp;
        deltaMoneyChosen(runTrials_idx, iS) = (money_nonDef_tmp - money_default_tmp).*(choice_highE_tmp == 1) +...
            -(money_nonDef_tmp - money_default_tmp).*(choice_highE_tmp == 0);
        
        %% effort
        E_highE(runTrials_idx, iS) = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        Echosen(runTrials_idx, iS) = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        switch task_fullName
            case 'physical'
                deltaEp(runTrials_idx, iS) = E_highE(runTrials_idx, iS);
                deltaEm(runTrials_idx, iS) = 0;
                % physical fatigue
                [physical_Fatigue(runTrials_idx, iS)] = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                % mental facilitation
                mental_Facilitation(runTrials_idx, iS) = 0;
            case 'mental'
                deltaEp(runTrials_idx, iS) = 0;
                deltaEm(runTrials_idx, iS) = E_highE(runTrials_idx, iS);
                % physical fatigue
                physical_Fatigue(runTrials_idx, iS) = 0;
                % mental facilitation
                [~,mental_Facilitation(runTrials_idx, iS),~,~] = extract_mental_previous_efficacy(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        end
    end % run loop
    
    %% preparing GLM
    potentialRegs.run1_cstt = run1_cstt(:,iS);
    potentialRegs.run2_cstt = run2_cstt(:,iS);
    potentialRegs.run3_cstt = run3_cstt(:,iS);
    potentialRegs.run4_cstt = run4_cstt(:,iS);
    potentialRegs.Em_cstt = Em_cstt(:,iS);
    potentialRegs.Ep_cstt = Ep_cstt(:,iS);
    potentialRegs.choice_bias = choice_bias(:,iS);
    potentialRegs.trialN = trialN(:,iS);
    potentialRegs.deltaEffort = E_highE(:,iS);
    potentialRegs.deltaEchosen = Echosen(:,iS);
    potentialRegs.deltaEp = deltaEp(:,iS);
    potentialRegs.deltaEm = deltaEm(:,iS);
    potentialRegs.RP = RP(:,iS);
    potentialRegs.deltaR = deltaR(:,iS);
    potentialRegs.deltaP = deltaP(:,iS);
    potentialRegs.deltaMoney = deltaMoney(:,iS);
    potentialRegs.deltaMoneyChosen = deltaMoneyChosen(:,iS);
    potentialRegs.physical_Fatigue = physical_Fatigue(:,iS).*deltaEp(:,iS);
    potentialRegs.mental_Facilitation = mental_Facilitation(:,iS).*deltaEm(:,iS);
    % extract regressors
    [x_regs, reg_names] = choice_GLM_regs(GLMprm, potentialRegs);
    
    %% performing GLM
    okTrials = ~isnan(choice_hE(:, iS));
    [betas_tmp,~,stats_tmp] = glmfit(x_regs(okTrials, :), choice_hE(okTrials, iS),...
        'binomial','constant','off','link','logit');
    NV_fit(okTrials,iS) = glmval(betas_tmp, x_regs(okTrials, :),...
        'identity','constant','off'); % fit before going through the sigmoid
    choice_hE_fit(okTrials,iS) = glmval(betas_tmp, x_regs(okTrials, :),...
        'logit','constant','off');
    
    % distribute data for each run
    for iRun = 1:nRuns
        run_trial_idx = (1:nTrialsPerRun) + nTrialsPerRun*(iRun - 1);
        choice_fit_perRun.(['run',num2str(iRun)])(:,iS) = choice_hE_fit(run_trial_idx,iS);
    end % run loop
    
    % store betas
    for iReg = 1:length(reg_names)
        curr_reg_nm = reg_names{iReg};
        betas.(curr_reg_nm)(iS) = betas_tmp(iReg);
        pval.(curr_reg_nm)(iS) = stats_tmp.p(iReg);
    end % regressor loop
    
    %% extract binned data with net value
    [choice_f_NV(:, iS), NV_f_NV(:, iS)] = do_bin2(choice_hE(:,iS), NV_fit(:,iS), nBins, 0);
    [choice_fit_f_NV(:, iS), NV_f_NV(:, iS)] = do_bin2(choice_hE_fit(:,iS), NV_fit(:,iS), nBins, 0);
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
[bins_choice.mean.NV_f_NV, bins_choice.sem.NV_f_NV] = mean_sem_sd(NV_f_NV, 2);
[bins_choice.mean.choice_f_NV, bins_choice.sem.choice_f_NV] = mean_sem_sd(choice_f_NV, 2);
[bins_choice.mean.choice_fit_f_NV, bins_choice.sem.choice_fit_f_NV] = mean_sem_sd(choice_fit_f_NV, 2);

%% display figure
if figDisp == 1
    pSize = 50;
    lWidth = 3;

    % show choice = f(NV)
    fig;
    hold on;
    choice_hdl = errorbar(bins_choice.mean.NV_f_NV,...
        bins_choice.mean.choice_f_NV, bins_choice.sem.choice_f_NV);
    choice_fit_hdl = errorbar(bins_choice.mean.NV_f_NV,...
        bins_choice.mean.choice_fit_f_NV, bins_choice.sem.choice_fit_f_NV);
    % improve visualization
    choice_hdl.LineWidth = lWidth;
    choice_fit_hdl.LineWidth = lWidth;
    legend([choice_hdl, choice_fit_hdl],{'choice','choice fitted'},'Location','Northeast');
    legend('boxoff');
    xlabel('Net value');
    ylabel('Choice (%)');
    legend_size(pSize);
end % figure display

end % function