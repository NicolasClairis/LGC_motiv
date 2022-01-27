function[avg_nonDefaultChoice, sd_nonDefaultChoice, avg_conf, avg_RT, avg_Eperformance] = checkChoiceNonDefaultProportion(subid, figDisp, n_bins)
%[avg_defaultChoice, sd_defaultChoice, avg_conf, avg_RT, avg_Eperformance] = checkChoiceNonDefaultProportion(subid, figDisp, n_bins)
% check choices of the non-default option in function of the other parameters of the task
% (fatigue, difficulty, physical/mental, reward/punishment, etc.) for each
% individual.
%
% INPUTS
% subid: string with subject identification number
%
% figDisp:
% (0) do not display individual figures data
% (1) display all figures for the current individual
%
% n_bins: number of bins for the analysis of choice = f(fatigue)
%
% OUTPUTS
% avg_nonDefaultChoice: structure with average choice of the non-default
% option in function of the other task parameters
%
% sd_nonDefaultChoice: same as avg_nonDefaultChoice but for standard
% deviation
%
% avg_conf: structure with average confidence on the choice in function of
% other task parameters
%
% avg_RT: structure with average reaction time (RT) of choices in function
% of other task parameters
%
% avg_Eperformance: structure with average effort performance in function
% of other task parameters


%% initialize all variables of interest
if ~exist('subid','var') || isempty(subid)
    subid = '074';
    disp(['data for CID',subid]);
end
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
    disp('display of individual figures');
end
nRuns = 4;
n_R_levels = 4;
n_E_levels = 4;
nTrials = 54;
if ~exist('n_bins','var') || isempty(n_bins)
    n_bins = 6;
    disp(['n_bins = ',num2str(n_bins),' by default']);
end
nBinTrials = nTrials/n_bins;
[percentageChoiceNonDefault,...
    percentage_R_choiceNonDefault,...
    percentage_P_choiceNonDefault,...
    n_failedChoices] = deal(NaN(1,nRuns));
[choiceNonDefault, choice,...
    confidence, RT, Eperformance] = deal( NaN(nRuns,nTrials));
[choiceNonDefault_R,...
    choiceNonDefault_P] = deal( NaN(nRuns,nTrials/2));

% prepare analysis per effort level
nRepeatsPerEffortLevel = nTrials/(n_E_levels - 1);
for iE = 1:(n_E_levels - 1)
    E_nm = ['E_level_',num2str(iE)];
    [choiceNonDefault_perElevel.(E_nm),...
        conf_perElevel.(E_nm),...
        RT_perElevel.(E_nm),...
        conf_perElevel.(E_nm),...
        Eperformance_perElevel.(E_nm)] = deal(NaN(nRuns, nRepeatsPerEffortLevel));
    [choiceNonDefault_perElevel.(['avg_E_',num2str(iE)]),...
        conf_perElevel.(['avg_E_',num2str(iE)]),...
        RT_perElevel.(['avg_E_',num2str(iE)]),...
        Eperformance_perElevel.(['avg_E_',num2str(iE)])] = deal(NaN(1,nRuns));
end
% prepare analysis per money level for unsigned (R+P) and signed (P/R)
% money levels
% 1) unsigned money levels
nRepeatsPerAbsDeltaMoneyLevel = nTrials/(n_R_levels - 1);
for iAbsDeltaMoney = 1:(n_R_levels - 1)
    deltaMoney_nm = ['deltaMoney_level_',num2str(iAbsDeltaMoney)];
    [choiceNonDefault_perDeltaMoneylevel.(deltaMoney_nm),...
        conf_perDeltaMoneylevel.(deltaMoney_nm),...
        RT_perDeltaMoneylevel.(deltaMoney_nm),...
        Eperformance_perDeltaMoneylevel.(deltaMoney_nm)] = deal(NaN(nRuns, nRepeatsPerAbsDeltaMoneyLevel));
    [choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)]),...
        conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)]),...
        RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)]),...
        Eperformance_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])] = deal(NaN(1,nRuns));
end
% 2) signed money levels
nRepeatsPerSignedMoneyLevel = (nTrials/(n_E_levels - 1))*(1/2);
for iMoney = [-(n_R_levels - 1):(-1), 1:(n_R_levels - 1)]
    if iMoney < 0
        signedMoney_nm = ['P_level_',num2str(abs(iMoney))];
    elseif iMoney > 0
        signedMoney_nm = ['R_level_',num2str(iMoney)];
    end
    jSignedMoney.(signedMoney_nm) = 0;
    [choiceNonDefault_perSignedMoneylevel.(signedMoney_nm),...
        conf_perSignedMoneylevel.(signedMoney_nm),...
        RT_perSignedMoneylevel.(signedMoney_nm)] = deal(NaN(nRuns, nRepeatsPerSignedMoneyLevel));
    [choiceNonDefault_perSignedMoneylevel.(['avg_',signedMoney_nm]),...
        conf_perSignedMoneylevel.(signedMoney_nm),...
        RT_perSignedMoneylevel.(signedMoney_nm),...
        Eperformance_perSignedMoneylevel.(signedMoney_nm)] = deal(NaN(1,nRuns));
end

%% extract data trial by trial
[choiceNonDefault_f_time,...
    conf_f_time,...
    RT_f_time,...
    Eperformance_f_time] = deal(NaN(nRuns, n_bins));
for iRun = 1:nRuns
    filenm = ls(['CID',subid,'_session',num2str(iRun),'_*_task_messyAllStuff.mat']);
    runData = load(filenm);
    defaultSide = runData.choiceOptions.default_LR;
    [absMoney_NonDefault_level, E_nonDefault_level] = deal(NaN(1,nTrials));
    R_or_P = runData.choiceOptions.R_or_P;
    %% extract choice
    choice(iRun,:) = runData.perfSummary.choice;
    choiceSide = choice(iRun,:);
    choiceSide(choice(iRun,:) >= 1) = 1;
    choiceSide(choice(iRun,:) <= -1) = -1;
    %% extract confidence
    conf_tmp = runData.perfSummary.confidence.lowOrHigh;
    confidence(iRun,:) = conf_tmp;
    confidence(iRun, conf_tmp == 0) = NaN; % no answer provided
    confidence(iRun, conf_tmp == 1) = 0; % low confidence
    confidence(iRun, conf_tmp == 2) = 1; % high confidence
    %% extract RT
    RT_tmp = runData.perfSummary.onsets.choice - runData.perfSummary.onsets.dispChoiceOptions;
    RT(iRun,:) = RT_tmp;
    %% extract effort performance per trial and rescale between 0 and 100%
    Eperf_tmp = runData.perfSummary.percentagePerf.*100;
    Eperformance(iRun,:) = Eperf_tmp;
    
    %% reinitialize all counters
    [jR, jP] = deal(0);
    for iE = 1:(n_E_levels - 1)
        E_nm = ['E_level_',num2str(iE)];
        jE.(E_nm) = 0;
    end
    for iAbsDeltaMoney = 1:(n_R_levels - 1)
        deltaMoney_nm = ['deltaMoney_level_',num2str(iAbsDeltaMoney)];
        jDeltaMoney.(deltaMoney_nm) = 0;
        jSignedMoney.(['R_level_',num2str(iAbsDeltaMoney)]) = 0;
        jSignedMoney.(['P_level_',num2str(iAbsDeltaMoney)]) = 0;
    end
    
    for iTrial = 1:nTrials
        %% extract effort and money level of the non-default option for each
        % trial
        switch defaultSide(iTrial)
            case -1
                E_nonDefault_level(iTrial) = runData.choiceOptions.E.right(iTrial);
                absMoney_NonDefault_level(iTrial) = runData.choiceOptions.R.right(iTrial);
            case 1
                E_nonDefault_level(iTrial) = runData.choiceOptions.E.left(iTrial);
                absMoney_NonDefault_level(iTrial) = runData.choiceOptions.R.left(iTrial);
        end
        E_nm = ['E_level_',num2str(E_nonDefault_level(iTrial))];
        jE.(E_nm) = jE.(E_nm) + 1;
        switch R_or_P{iTrial}
            case 'R'
                jR = jR + 1;
                signedMoney_nm = ['R_level_',num2str(absMoney_NonDefault_level(iTrial))];
                deltaMoney_nm = ['deltaMoney_level_',num2str(absMoney_NonDefault_level(iTrial))];
            case 'P'
                jP = jP + 1;
                signedMoney_nm = ['P_level_',num2str(absMoney_NonDefault_level(iTrial))];
                deltaMoney_nm = ['deltaMoney_level_',num2str(n_R_levels-absMoney_NonDefault_level(iTrial))];
        end
        jSignedMoney.(signedMoney_nm) = jSignedMoney.(signedMoney_nm) + 1;
        jDeltaMoney.(deltaMoney_nm) = jDeltaMoney.(deltaMoney_nm) + 1;
        
        %% did the participant choose the default or the non-default option
        % for the current trial?
        if choiceSide(iTrial) == -defaultSide(iTrial) % choice non-default option
            % extract average selection of the default option
            choiceNonDefault(iRun, iTrial) = 1;
            
            % extract default choices per money level
            choiceNonDefault_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = 1;
            choiceNonDefault_perDeltaMoneylevel.(deltaMoney_nm)(iRun, jDeltaMoney.(deltaMoney_nm)) = 1;
            
            % extract average reward (or respectively punishment) default choices
            switch R_or_P{iTrial}
                case 'R'
                    choiceNonDefault_R(iRun,jR) = 1;
                case 'P'
                    choiceNonDefault_P(iRun,jP) = 1;
            end
            
            % extract choice and perf per effort level
            choiceNonDefault_perElevel.(E_nm)(iRun, jE.(E_nm)) = 1;
        elseif choiceSide(iTrial) == defaultSide(iTrial) % choice default option
            % extract average selection of the default option
            choiceNonDefault(iRun, iTrial) = 0;
            
            % extract default choices per money level
            choiceNonDefault_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = 0;
            choiceNonDefault_perDeltaMoneylevel.(deltaMoney_nm)(iRun, jDeltaMoney.(deltaMoney_nm)) = 0;
            
            % extract average reward (or respectively punishment) default choices
            switch R_or_P{iTrial}
                case 'R'
                    choiceNonDefault_R(iRun,jR) = 0;
                case 'P'
                    choiceNonDefault_P(iRun,jP) = 0;
            end
            
            % extract choice and perf per effort level
            choiceNonDefault_perElevel.(E_nm)(iRun, jE.(E_nm)) = 0;
        else % no choice was made = ignore the trial
            
        end
        
        %% level of confidence depending on money/effort
        conf_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = confidence(iRun,iTrial);
        conf_perDeltaMoneylevel.(deltaMoney_nm)(iRun, jDeltaMoney.(deltaMoney_nm)) =  confidence(iRun, iTrial);
        conf_perElevel.(E_nm)(iRun, jE.(E_nm)) = confidence(iRun,iTrial);
        
        %% level of reaction times (RT) depending on money/effort
        RT_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = RT(iRun,iTrial);
        RT_perDeltaMoneylevel.(deltaMoney_nm)(iRun, jDeltaMoney.(deltaMoney_nm)) =  RT(iRun, iTrial);
        RT_perElevel.(E_nm)(iRun, jE.(E_nm)) = RT(iRun,iTrial);
        
        %% level of performance during effort period depending on money/effort
        Eperformance_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = Eperformance(iRun,iTrial);
        Eperformance_perDeltaMoneylevel.(deltaMoney_nm)(iRun, jDeltaMoney.(deltaMoney_nm)) =  Eperformance(iRun, iTrial);
        Eperformance_perElevel.(E_nm)(iRun, jE.(E_nm)) = Eperformance(iRun,iTrial);
    end % trial loop
    
    %% ratio of chosen the default option per run
    percentageChoiceNonDefault(iRun) = sum(choiceNonDefault(iRun,:),'omitnan')/nTrials;
    percentage_R_choiceNonDefault(iRun) = sum(choiceNonDefault_R(iRun,:),'omitnan')/(nTrials/2);
    percentage_P_choiceNonDefault(iRun) = sum(choiceNonDefault_P(iRun,:),'omitnan')/(nTrials/2);
    
    %% default choice/confidence = f(time)
    for iBin = 1:n_bins
        trial_idx = (1:nBinTrials) + nBinTrials*(iBin - 1);
        choiceNonDefault_f_time(iRun, iBin) = mean( choiceNonDefault(iRun,trial_idx), 2,'omitnan' );
        conf_f_time(iRun, iBin) = mean( confidence(iRun,trial_idx), 2,'omitnan' );
        RT_f_time(iRun, iBin) = mean( RT(iRun,trial_idx), 2,'omitnan' );
        Eperformance_f_time(iRun, iBin) = mean( Eperformance(iRun,trial_idx), 2,'omitnan' );
    end
    
    %% sum of failed choices to check
    n_failedChoices(iRun) = sum(choice(iRun,:) == 0,'omitnan');
    
    %% percentage choice of default per effort level
    for iE = 1:(n_E_levels - 1)
        choiceNonDefault_perElevel.(['avg_E_',num2str(iE)])(iRun) = 100*sum(choiceNonDefault_perElevel.(['E_level_',num2str(iE)])(iRun, :),'omitnan')/nRepeatsPerEffortLevel;
        conf_perElevel.(['avg_E_',num2str(iE)])(iRun) = 100*mean(conf_perElevel.(['E_level_',num2str(iE)])(iRun, :),2,'omitnan');
        RT_perElevel.(['avg_E_',num2str(iE)])(iRun) = mean(RT_perElevel.(['E_level_',num2str(iE)])(iRun, :),2,'omitnan');
        Eperformance_perElevel.(['avg_E_',num2str(iE)])(iRun) = mean(Eperformance_perElevel.(['E_level_',num2str(iE)])(iRun, :),2,'omitnan');
    end
    
    %% percentage choice of default per money level
    for iAbsDeltaMoney = 1:(n_R_levels - 1)
        choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(iRun) = 100*sum(choiceNonDefault_perDeltaMoneylevel.(['deltaMoney_level_',num2str(iAbsDeltaMoney)])(iRun, :),'omitnan')/nRepeatsPerAbsDeltaMoneyLevel;
        choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(iRun) = 100*sum(choiceNonDefault_perSignedMoneylevel.(['R_level_',num2str(iAbsDeltaMoney)])(iRun, :),'omitnan')/nRepeatsPerSignedMoneyLevel;
        choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(iRun) = 100*sum(choiceNonDefault_perSignedMoneylevel.(['P_level_',num2str(iAbsDeltaMoney)])(iRun, :),'omitnan')/nRepeatsPerSignedMoneyLevel;
        % same for confidence
        conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(iRun) = mean(100*conf_perDeltaMoneylevel.(['deltaMoney_level_',num2str(iAbsDeltaMoney)])(iRun, :),2,'omitnan');
        conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(iRun) = mean(100*conf_perSignedMoneylevel.(['R_level_',num2str(iAbsDeltaMoney)])(iRun, :),2,'omitnan');
        conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(iRun) = mean(100*conf_perSignedMoneylevel.(['P_level_',num2str(iAbsDeltaMoney)])(iRun, :),2,'omitnan');
        % same for RT
        RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(iRun) = mean(RT_perDeltaMoneylevel.(['deltaMoney_level_',num2str(iAbsDeltaMoney)])(iRun, :),2,'omitnan');
        RT_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(iRun) = mean(RT_perSignedMoneylevel.(['R_level_',num2str(iAbsDeltaMoney)])(iRun, :),2,'omitnan');
        RT_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(iRun) = mean(RT_perSignedMoneylevel.(['P_level_',num2str(iAbsDeltaMoney)])(iRun, :),2,'omitnan');
        % same for effort performance
        Eperformance_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(iRun) = mean(Eperformance_perDeltaMoneylevel.(['deltaMoney_level_',num2str(iAbsDeltaMoney)])(iRun, :),2,'omitnan');
        Eperformance_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(iRun) = mean(Eperformance_perSignedMoneylevel.(['R_level_',num2str(iAbsDeltaMoney)])(iRun, :),2,'omitnan');
        Eperformance_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(iRun) = mean(Eperformance_perSignedMoneylevel.(['P_level_',num2str(iAbsDeltaMoney)])(iRun, :),2,'omitnan');
    end
end % run loop

%% extract task type order
if exist(['CID',subid,'_session1_mental_task_behavioral_tmp.mat'],'file') &&...
        exist(['CID',subid,'_session2_physical_task_behavioral_tmp.mat'],'file') &&...
        exist(['CID',subid,'_session3_mental_task_behavioral_tmp.mat'],'file') &&...
        exist(['CID',subid,'_session4_physical_task_behavioral_tmp.mat'],'file')
    Em_runs = [1,3];
    Ep_runs = [2,4];
elseif exist(['CID',subid,'_session1_physical_task_behavioral_tmp.mat'],'file') &&...
        exist(['CID',subid,'_session2_mental_task_behavioral_tmp.mat'],'file') &&...
        exist(['CID',subid,'_session3_physical_task_behavioral_tmp.mat'],'file') &&...
        exist(['CID',subid,'_session4_mental_task_behavioral_tmp.mat'],'file')
    Em_runs = [2,4];
    Ep_runs = [1,3];
else
    error('couldn''t find participant task types');
end

%% extract average
% ratio of choosing default per task
avg_nonDefaultChoice.Em = mean(percentageChoiceNonDefault(Em_runs).*100);
avg_nonDefaultChoice.Ep = mean(percentageChoiceNonDefault(Ep_runs).*100);
% ratio of choosing default per R/P condition per task
avg_nonDefaultChoice.Em_R = mean(percentage_R_choiceNonDefault(Em_runs).*100,2,'omitnan');
avg_nonDefaultChoice.Ep_R = mean(percentage_R_choiceNonDefault(Ep_runs).*100,2,'omitnan');
avg_nonDefaultChoice.Em_P = mean(percentage_P_choiceNonDefault(Em_runs).*100,2,'omitnan');
avg_nonDefaultChoice.Ep_P = mean(percentage_P_choiceNonDefault(Ep_runs).*100,2,'omitnan');
% average across tasks
avg_nonDefaultChoice.R = mean(percentage_R_choiceNonDefault.*100,2,'omitnan');
avg_nonDefaultChoice.P = mean(percentage_P_choiceNonDefault.*100,2,'omitnan');

% per task per effort level
for iE = 1:(n_E_levels - 1)
    avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]) = mean(choiceNonDefault_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]) = mean(choiceNonDefault_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
    avg_conf.perElevel.(['Em_',num2str(iE)]) = mean(conf_perElevel.(['avg_E_',num2str(iE)])(Em_runs),'omitnan');
    avg_conf.perElevel.(['Ep_',num2str(iE)]) = mean(conf_perElevel.(['avg_E_',num2str(iE)])(Ep_runs),'omitnan');
    avg_RT.perElevel.(['Em_',num2str(iE)]) = mean(RT_perElevel.(['avg_E_',num2str(iE)])(Em_runs),'omitnan');
    avg_RT.perElevel.(['Ep_',num2str(iE)]) = mean(RT_perElevel.(['avg_E_',num2str(iE)])(Ep_runs),'omitnan');
    avg_Eperformance.perElevel.(['Em_',num2str(iE)]) = mean(Eperformance_perElevel.(['avg_E_',num2str(iE)])(Em_runs),'omitnan');
    avg_Eperformance.perElevel.(['Ep_',num2str(iE)]) = mean(Eperformance_perElevel.(['avg_E_',num2str(iE)])(Ep_runs),'omitnan');
    % STD
    sd_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]) = std(choiceNonDefault_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]) = std(choiceNonDefault_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
    sd_conf.perElevel.(['Em_',num2str(iE)]) = std(conf_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_conf.perElevel.(['Ep_',num2str(iE)]) = std(conf_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
    sd_RT.perElevel.(['Em_',num2str(iE)]) = std(RT_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_RT.perElevel.(['Ep_',num2str(iE)]) = std(RT_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
    sd_Eperformance.perElevel.(['Em_',num2str(iE)]) = std(Eperformance_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_Eperformance.perElevel.(['Ep_',num2str(iE)]) = std(Eperformance_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
end

% per task per money level
for iAbsDeltaMoney = 1:(n_R_levels - 1)
    avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsDeltaMoney)]) = mean(choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsDeltaMoney)]) = mean(choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_nonDefaultChoice.perSignedMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = mean(choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = mean(choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_nonDefaultChoice.perSignedMoneylevel.(['Em_P_',num2str(iAbsDeltaMoney)]) = mean(choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_P_',num2str(iAbsDeltaMoney)]) = mean(choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_conf.perDeltaMoneylevel.(['Em_',num2str(iAbsDeltaMoney)]) = mean(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_conf.perDeltaMoneylevel.(['Ep_',num2str(iAbsDeltaMoney)]) = mean(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_conf.perSignedMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = mean(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_conf.perSignedMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = mean(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_conf.perSignedMoneylevel.(['Em_P_',num2str(iAbsDeltaMoney)]) = mean(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_conf.perSignedMoneylevel.(['Ep_P_',num2str(iAbsDeltaMoney)]) = mean(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_RT.perDeltaMoneylevel.(['Em_',num2str(iAbsDeltaMoney)]) = mean(RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_RT.perDeltaMoneylevel.(['Ep_',num2str(iAbsDeltaMoney)]) = mean(RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_RT.perSignedMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = mean(RT_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_RT.perSignedMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = mean(RT_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_RT.perSignedMoneylevel.(['Em_P_',num2str(iAbsDeltaMoney)]) = mean(RT_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_RT.perSignedMoneylevel.(['Ep_P_',num2str(iAbsDeltaMoney)]) = mean(RT_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_Eperformance.perDeltaMoneylevel.(['Em_',num2str(iAbsDeltaMoney)]) = mean(Eperformance_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_Eperformance.perDeltaMoneylevel.(['Ep_',num2str(iAbsDeltaMoney)]) = mean(Eperformance_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_Eperformance.perSignedMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = mean(Eperformance_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_Eperformance.perSignedMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = mean(Eperformance_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Ep_runs));
    avg_Eperformance.perSignedMoneylevel.(['Em_P_',num2str(iAbsDeltaMoney)]) = mean(Eperformance_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Em_runs));
    avg_Eperformance.perSignedMoneylevel.(['Ep_P_',num2str(iAbsDeltaMoney)]) = mean(Eperformance_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Ep_runs));
    % STD
    sd_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsDeltaMoney)]) = std(choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsDeltaMoney)]) = std(choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_nonDefaultChoice.perSignedMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = std(choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_nonDefaultChoice.perSignedMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = std(choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_nonDefaultChoice.perSignedMoneylevel.(['Em_P_',num2str(iAbsDeltaMoney)]) = std(choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_nonDefaultChoice.perSignedMoneylevel.(['Ep_P_',num2str(iAbsDeltaMoney)]) = std(choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_conf.perDeltaMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = std(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_conf.perDeltaMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = std(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_conf.perSignedMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = std(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_conf.perSignedMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = std(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_conf.perSignedMoneylevel.(['Em_P_',num2str(iAbsDeltaMoney)]) = std(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_conf.perSignedMoneylevel.(['Ep_P_',num2str(iAbsDeltaMoney)]) = std(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_RT.perDeltaMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = std(RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_RT.perDeltaMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = std(RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_RT.perSignedMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = std(RT_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_RT.perSignedMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = std(RT_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_RT.perSignedMoneylevel.(['Em_P_',num2str(iAbsDeltaMoney)]) = std(RT_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_RT.perSignedMoneylevel.(['Ep_P_',num2str(iAbsDeltaMoney)]) = std(RT_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_Eperformance.perDeltaMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = std(Eperformance_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_Eperformance.perDeltaMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = std(Eperformance_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_Eperformance.perSignedMoneylevel.(['Em_R_',num2str(iAbsDeltaMoney)]) = std(Eperformance_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_Eperformance.perSignedMoneylevel.(['Ep_R_',num2str(iAbsDeltaMoney)]) = std(Eperformance_perSignedMoneylevel.(['avg_R_',num2str(iAbsDeltaMoney)])(Ep_runs));
    sd_Eperformance.perSignedMoneylevel.(['Em_P_',num2str(iAbsDeltaMoney)]) = std(Eperformance_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Em_runs));
    sd_Eperformance.perSignedMoneylevel.(['Ep_P_',num2str(iAbsDeltaMoney)]) = std(Eperformance_perSignedMoneylevel.(['avg_P_',num2str(iAbsDeltaMoney)])(Ep_runs));
end

avg_nonDefaultChoice.Em_f_time = mean(choiceNonDefault_f_time(Em_runs,:),1,'omitnan');
avg_nonDefaultChoice.Ep_f_time = mean(choiceNonDefault_f_time(Ep_runs,:),1,'omitnan');
avg_conf.Em_f_time = mean(conf_f_time(Em_runs,:),1,'omitnan');
avg_conf.Ep_f_time = mean(conf_f_time(Ep_runs,:),1,'omitnan');
avg_RT.Em_f_time = mean(RT_f_time(Em_runs,:),1,'omitnan');
avg_RT.Ep_f_time = mean(RT_f_time(Ep_runs,:),1,'omitnan');
avg_Eperformance.Em_f_time = mean(Eperformance_f_time(Em_runs,:),1,'omitnan');
avg_Eperformance.Ep_f_time = mean(Eperformance_f_time(Ep_runs,:),1,'omitnan');

%% figures
if figDisp == 1
    % figure parameters
    lWidth = 3;
    pSize = 30;
    bWidth = 0.4;
    bDist = 0.2;
    
    % check choices = f(fatigue, R/P, level of R, level of E)
    %% check choices = f(fatigue)
    fig;
    % mark the 50% trait
    plot(0:n_bins+1, 50*ones(1,n_bins+2),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    runType = cell(1,nRuns);
    for iRun = 1:nRuns
        % different colour for mental and physical effort and for each run
        switch iRun
            case Em_runs(1)
                lColour = [0 255 0]./255;
                runType{iRun} = ['run ',num2str(iRun),' mental'];
            case Em_runs(2)
                lColour = [0 143 0]./255;
                runType{iRun} = ['run ',num2str(iRun),' mental'];
            case Ep_runs(1)
                lColour = [0 0 255]./255;
                runType{iRun} = ['run ',num2str(iRun),' physical'];
            case Ep_runs(2)
                lColour = [0 0 143]./255;
                runType{iRun} = ['run ',num2str(iRun),' physical'];
        end
        bar_hdl.(['run',num2str(iRun)]) = plot(1:n_bins, choiceNonDefault_f_time(iRun, :)*100,...
            'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    end
    ylim([0 100]);
    xlim([0 n_bins+1]);
    xlabel('trial bins');
    ylabel('Choice non-default option (%)');
    legend([bar_hdl.run1, bar_hdl.run2, bar_hdl.run3, bar_hdl.run4],...
        runType);
    legend('boxoff');
    legend('Location','SouthWest');
    legend_size(pSize);
    
    %% check choices = f(R/P)
    fig;
    % mark the 50% trait
    plot(0:3, 50*ones(1,4),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    bar(1, mean(percentage_R_choiceNonDefault.*100),'FaceColor',[0 153/255 1]);
    errorbar(1, mean(percentage_R_choiceNonDefault.*100), std(percentage_R_choiceNonDefault.*100),'k');
    bar(2, mean(percentage_P_choiceNonDefault.*100),'FaceColor',[1 153/255 0]);
    errorbar(2, mean(percentage_P_choiceNonDefault.*100), std(percentage_P_choiceNonDefault.*100),'k');
    xticks(1:2);
    xticklabels({'R','P'});
    ylim([0 100]);
    xlim([0 3]);
    ylabel('Choice non-default option (%)');
    legend_size(pSize);
    
    %% check choices = f(|high money level - default money level|)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iDeltaMoney = 1:(n_R_levels - 1)
        Em_bar_hdl = bar(iDeltaMoney-bDist, avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iDeltaMoney)]),...
            'FaceColor','g','BarWidth',bWidth);
        Ep_bar_hdl = bar(iDeltaMoney+bDist, avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iDeltaMoney)]),...
            'FaceColor','b','BarWidth',bWidth);
        errorbar(iDeltaMoney-bDist,...
            avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iDeltaMoney)]),...
            sd_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iDeltaMoney)]),...
            'k')
        errorbar(iDeltaMoney+bDist,...
            avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iDeltaMoney)]),...
            sd_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iDeltaMoney)]),...
            'k')
    end
    ylim([0 100]);
    ylabel('Choice non-default option (%)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_R_levels]);
    xlabel('|Δ money| level');
    legend([Em_bar_hdl, Ep_bar_hdl],'mental','physical');
    legend('boxoff');
    legend('Location','NorthWest');
    legend_size(pSize);
    
    %% check choices = f(money levels (splitting R and P trials))
    fig;
    % mark the 50% trait
    plot(-n_bins:n_bins, 50*ones(1,(n_bins*2)+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
        jMoney = abs(iMoney);
        if iMoney < 0
            taskCond = 'P';
        elseif iMoney > 0
            taskCond = 'R';
        end
        Em_bar_hdl = bar(iMoney-bDist,...
            avg_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            'FaceColor','g','BarWidth',bWidth);
        Ep_bar_hdl = bar(iMoney+bDist,...
            avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            'FaceColor','b','BarWidth',bWidth);
        errorbar(iMoney-bDist,...
            avg_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            sd_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            'k')
        errorbar(iMoney+bDist,...
            avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            sd_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            'k')
    end
    ylim([0 100]);
    ylabel('Choice non-default option (%)');
    xticks([-3:-(1), 1:3]);
    xticklabels({'-3','-2','-1','1','2','3'});
    xlim([-n_R_levels n_R_levels]);
    xlabel('Money level');
    legend_size(pSize);
    ylim([0 100]);
    ylabel('Choice non-default option (%)');
    legend([Em_bar_hdl, Ep_bar_hdl],'mental','physical');
    legend('boxoff');
    legend('Location','NorthWest');
    legend_size(pSize);
    
    %% check choices = f(E level)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iE = 1:(n_E_levels - 1)
        Em_bar_hdl = bar(iE-bDist, avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]),...
            'FaceColor','g','BarWidth',bWidth);
        Ep_bar_hdl = bar(iE+bDist, avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]),...
            'FaceColor','b','BarWidth',bWidth);
        errorbar(iE-bDist, avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]),...
            sd_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]),'k')
        errorbar(iE+bDist, avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]),...
            sd_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]),'k')
    end
    ylim([0 100]);
    ylabel('Choice non-default option (%)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_E_levels]);
    xlabel('Effort level');
    legend([Em_bar_hdl, Ep_bar_hdl],'mental','physical');
    legend('boxoff');
    legend('Location','NorthWest');
    legend_size(pSize);
    
    %% check confidence = f(R/P levels)
    fig;
    % mark the 50% trait
    plot(-n_bins:n_bins, 50*ones(1,(n_bins*2)+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
        jMoney = abs(iMoney);
        if iMoney < 0
            taskCond = 'P';
        elseif iMoney > 0
            taskCond = 'R';
        end
        Em_bar_hdl = bar(iMoney-bDist,...
            avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            'FaceColor','g','BarWidth',bWidth);
        Ep_bar_hdl = bar(iMoney+bDist,...
            avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            'FaceColor','b','BarWidth',bWidth);
        errorbar(iMoney-bDist,...
            avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            sd_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            'k')
        errorbar(iMoney+bDist,...
            avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            sd_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            'k')
    end
    xticks([-3:-(1), 1:3]);
    xticklabels({'-3','-2','-1','1','2','3'});
    xlim([-n_R_levels n_R_levels]);
    ylim([0 100]);
    ylabel('Level of confidence');
    xlabel('Money level');
    legend([Em_bar_hdl, Ep_bar_hdl],'mental','physical');
    legend('boxoff');
    legend('Location','NorthWest');
    legend_size(pSize);
    
    %% check confidence = f(E levels)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iE = 1:(n_E_levels - 1)
        Em_bar_hdl = bar(iE-bDist, avg_conf.perElevel.(['Em_',num2str(iE)]),...
            'FaceColor','g','BarWidth',bWidth);
        Ep_bar_hdl = bar(iE+bDist, avg_conf.perElevel.(['Ep_',num2str(iE)]),...
            'FaceColor','b','BarWidth',bWidth);
        errorbar(iE-bDist,...
            avg_conf.perElevel.(['Em_',num2str(iE)]),...
            sd_conf.perElevel.(['Em_',num2str(iE)]),'k')
        errorbar(iE+bDist,...
            avg_conf.perElevel.(['Ep_',num2str(iE)]),...
            sd_conf.perElevel.(['Ep_',num2str(iE)]),'k')
    end
    ylim([0 100]);
    ylabel('Confidence');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_E_levels]);
    xlabel('Effort level');
    legend_size(pSize);
    legend([Em_bar_hdl, Ep_bar_hdl],'mental','physical');
    legend('boxoff');
    legend('Location','NorthWest');
    
    %% check confidence = f(time)
    fig;
    % mark the 50% trait
    plot(1:n_bins, 50*ones(1,n_bins),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    runType = cell(1,nRuns);
    for iRun = 1:nRuns
        % different colour for mental and physical effort and for each run
        switch iRun
            case Em_runs(1)
                lColour = [0 255 0]./255;
                runType{iRun} = ['run ',num2str(iRun),' mental'];
            case Em_runs(2)
                lColour = [0 143 0]./255;
                runType{iRun} = ['run ',num2str(iRun),' mental'];
            case Ep_runs(1)
                lColour = [0 0 255]./255;
                runType{iRun} = ['run ',num2str(iRun),' physical'];
            case Ep_runs(2)
                lColour = [0 0 143]./255;
                runType{iRun} = ['run ',num2str(iRun),' physical'];
        end
        bar_hdl.(['run',num2str(iRun)]) = plot(1:n_bins, conf_f_time(iRun, :)*100,...
            'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    end
    ylim([0 100]);
    xlim([0 n_bins+1]);
    xlabel('trial bins');
    ylabel('Confidence (%)');
    legend([bar_hdl.run1, bar_hdl.run2, bar_hdl.run3, bar_hdl.run4],...
        runType);
    legend('boxoff');
    legend('Location','SouthEast');
    legend_size(pSize);
    
    %% check RT = f(R/P levels)
    fig;
    hold on;
    for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
        jMoney = abs(iMoney);
        if iMoney < 0
            taskCond = 'P';
        elseif iMoney > 0
            taskCond = 'R';
        end
        Em_bar_hdl = bar(iMoney-bDist,...
            avg_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            'FaceColor','g','BarWidth',bWidth);
        Ep_bar_hdl = bar(iMoney+bDist,...
            avg_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            'FaceColor','b','BarWidth',bWidth);
        errorbar(iMoney-bDist,...
            avg_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            sd_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            'k')
        errorbar(iMoney+bDist,...
            avg_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            sd_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            'k')
    end
    xticks([-3:-(1), 1:3]);
    xticklabels({'-3','-2','-1','1','2','3'});
    xlim([-n_R_levels n_R_levels]);
    ylim([0 5]);
    ylabel('RT (s)');
    xlabel('Money level');
    legend([Em_bar_hdl, Ep_bar_hdl],'mental','physical');
    legend('boxoff');
    legend('Location','NorthWest');
    legend_size(pSize);
    
    %% check RT = f(E levels)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iE = 1:(n_E_levels - 1)
        Em_bar_hdl = bar(iE-bDist, avg_RT.perElevel.(['Em_',num2str(iE)]),...
            'FaceColor','g','BarWidth',bWidth);
        Ep_bar_hdl = bar(iE+bDist, avg_RT.perElevel.(['Ep_',num2str(iE)]),...
            'FaceColor','b','BarWidth',bWidth);
        errorbar(iE-bDist, avg_RT.perElevel.(['Em_',num2str(iE)]),...
            sd_RT.perElevel.(['Em_',num2str(iE)]),'k')
        errorbar(iE+bDist, avg_RT.perElevel.(['Ep_',num2str(iE)]),...
            sd_RT.perElevel.(['Ep_',num2str(iE)]),'k')
    end
    ylim([0 5]);
    ylabel('RT (s)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_E_levels]);
    xlabel('Effort level');
    legend([Em_bar_hdl, Ep_bar_hdl],'mental','physical');
    legend('boxoff');
    legend('Location','NorthWest');
    legend_size(pSize);
    
    %% check RT = f(time)
    fig;
    % mark the 50% trait
    plot(1:n_bins, 50*ones(1,n_bins),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    runType = cell(1,nRuns);
    for iRun = 1:nRuns
        % different colour for mental and physical effort and for each run
        switch iRun
            case Em_runs(1)
                lColour = [0 255 0]./255;
                runType{iRun} = ['run ',num2str(iRun),' mental'];
            case Em_runs(2)
                lColour = [0 143 0]./255;
                runType{iRun} = ['run ',num2str(iRun),' mental'];
            case Ep_runs(1)
                lColour = [0 0 255]./255;
                runType{iRun} = ['run ',num2str(iRun),' physical'];
            case Ep_runs(2)
                lColour = [0 0 143]./255;
                runType{iRun} = ['run ',num2str(iRun),' physical'];
        end
        bar_hdl.(['run',num2str(iRun)]) = plot(1:n_bins, RT_f_time(iRun, :),...
            'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    end
    ylim([0 5]);
    xlim([0 n_bins+1]);
    xlabel('trial bins');
    ylabel('RT (s)');
    legend([bar_hdl.run1, bar_hdl.run2, bar_hdl.run3, bar_hdl.run4],...
        runType);
    legend('boxoff');
    legend('Location','NorthEast');
    legend_size(pSize);
    
    %% check effort performance = f(R/P levels)
    fig;
    hold on;
    for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
        jMoney = abs(iMoney);
        if iMoney < 0
            taskCond = 'P';
        elseif iMoney > 0
            taskCond = 'R';
        end
        Em_bar_hdl = bar(iMoney-bDist,...
            avg_Eperformance.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            'FaceColor','g','BarWidth',bWidth);
        Ep_bar_hdl = bar(iMoney+bDist,...
            avg_Eperformance.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            'FaceColor','b','BarWidth',bWidth);
        errorbar(iMoney-bDist,...
            avg_Eperformance.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            sd_Eperformance.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
            'k')
        errorbar(iMoney+bDist,...
            avg_Eperformance.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            sd_Eperformance.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
            'k')
    end
    xticks([-3:-(1), 1:3]);
    xticklabels({'-3','-2','-1','1','2','3'});
    xlim([-n_R_levels n_R_levels]);
    ylim([0 100]);
    ylabel('Effort performance (%)');
    xlabel('Money level');
    legend([Em_bar_hdl, Ep_bar_hdl],'mental','physical');
    legend('boxoff');
    legend('Location','NorthWest');
    legend_size(pSize);
    
    %% check effort performance = f(E levels)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iE = 1:(n_E_levels - 1)
        Em_bar_hdl = bar(iE-bDist,...
            avg_Eperformance.perElevel.(['Em_',num2str(iE)]),...
            'FaceColor','g','BarWidth',bWidth);
        Ep_bar_hdl = bar(iE+bDist,...
            avg_Eperformance.perElevel.(['Ep_',num2str(iE)]),...
            'FaceColor','b','BarWidth',bWidth);
        errorbar(iE-bDist,...
            avg_Eperformance.perElevel.(['Em_',num2str(iE)]),...
            sd_Eperformance.perElevel.(['Em_',num2str(iE)]),...
            'k')
        errorbar(iE+bDist,...
            avg_Eperformance.perElevel.(['Ep_',num2str(iE)]),...
            sd_Eperformance.perElevel.(['Ep_',num2str(iE)]),...
            'k')
    end
    ylim([0 100]);
    ylabel('Effort performance (%)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_E_levels]);
    xlabel('Effort level');
    legend([Em_bar_hdl, Ep_bar_hdl],'mental','physical');
    legend('boxoff');
    legend('Location','NorthWest');
    legend_size(pSize);
    
    %% check effort performance = f(time)
    fig;
    % mark the 50% trait
    plot(1:n_bins, 50*ones(1,n_bins),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    runType = cell(1,nRuns);
    for iRun = 1:nRuns
        % different colour for mental and physical effort and for each run
        switch iRun
            case Em_runs(1)
                lColour = [0 255 0]./255;
                runType{iRun} = ['run ',num2str(iRun),' mental'];
            case Em_runs(2)
                lColour = [0 143 0]./255;
                runType{iRun} = ['run ',num2str(iRun),' mental'];
            case Ep_runs(1)
                lColour = [0 0 255]./255;
                runType{iRun} = ['run ',num2str(iRun),' physical'];
            case Ep_runs(2)
                lColour = [0 0 143]./255;
                runType{iRun} = ['run ',num2str(iRun),' physical'];
        end
        bar_hdl.(['run',num2str(iRun)]) = plot(1:n_bins, Eperformance_f_time(iRun, :),...
            'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    end
    ylim([0 100]);
    xlim([0 n_bins+1]);
    xlabel('trial bins');
    ylabel('Effort performance (%)');
    legend([bar_hdl.run1, bar_hdl.run2, bar_hdl.run3, bar_hdl.run4],...
        runType);
    legend('boxoff');
    legend('Location','SouthEast');
    legend_size(pSize);
    
end % figure display

end % function display