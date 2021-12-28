function[avg_nonDefaultChoice, sd_nonDefaultChoice, avg_conf, avg_RT] = checkChoiceNonDefaultProportion(subid, figDisp, n_bins)
%[avg_defaultChoice, sd_defaultChoice, avg_conf, avg_RT] = checkChoiceNonDefaultProportion(subid, figDisp, n_bins)
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
[choiceNonDefault, choice, confidence, RT] = deal( NaN(nRuns,nTrials));
[choiceNonDefault_R,...
    choiceNonDefault_P] = deal( NaN(nRuns,nTrials/2));

% prepare analysis per effort level
nRepeatsPerEffortLevel = nTrials/(n_E_levels - 1);
for iE = 1:(n_E_levels - 1)
    E_nm = ['E_level_',num2str(iE)];
    [choiceNonDefault_perElevel.(E_nm),...
        conf_perElevel.(E_nm),...
        RT_perElevel.(E_nm)] = deal(NaN(nRuns, nRepeatsPerEffortLevel));
    [choiceNonDefault_perElevel.(['avg_E_',num2str(iE)]),...
        conf_perElevel.(['avg_E_',num2str(iE)]),...
        RT_perElevel.(['avg_E_',num2str(iE)])] = deal(NaN(1,nRuns));
end
% prepare analysis per money level for unsigned (R+P) and signed (P/R)
% money levels
% 1) unsigned money levels
nRepeatsPerAbsMoneyLevel = nTrials/(n_R_levels - 1);
for iAbsMoney = 1:(n_R_levels - 1)
    absMoney_nm = ['absMoney_level_',num2str(iAbsMoney)];
    deltaMoney_nm = ['deltaMoney_level_',num2str(iAbsMoney)];
    [choiceNonDefault_perAbsMoneylevel.(absMoney_nm),...
        conf_perAbsMoneylevel.(absMoney_nm),...
        RT_perAbsMoneylevel.(absMoney_nm),...
        choiceNonDefault_perDeltaMoneylevel.(deltaMoney_nm),...
        conf_perDeltaMoneylevel.(deltaMoney_nm),...
        RT_perDeltaMoneylevel.(deltaMoney_nm)] = deal(NaN(nRuns, nRepeatsPerAbsMoneyLevel));
    [choiceNonDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)]),...
        conf_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)]),...
        RT_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)]),...
        choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)]),...
        conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)]),...
        RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])] = deal(NaN(1,nRuns));
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
        RT_perSignedMoneylevel.(signedMoney_nm)] = deal(NaN(1,nRuns));
end

%% extract data trial by trial
[choiceNonDefault_f_time,...
    conf_f_time,...
    RT_f_time] = deal(NaN(nRuns, n_bins));
for iRun = 1:nRuns
    filenm = ls(['CID',subid,'_session',num2str(iRun),'_*_task_messyAllStuff.mat']);
    runData = load(filenm);
    defaultSide = runData.choiceOptions.default_LR;
    [absMoney_level, E_level] = deal(NaN(1,nTrials));
    R_or_P = runData.choiceOptions.R_or_P;
    % extract choice
    choice(iRun,:) = runData.perfSummary.choice;
    choiceSide = choice(iRun,:);
    choiceSide(choice(iRun,:) >= 1) = 1;
    choiceSide(choice(iRun,:) <= -1) = -1;
    % extract confidence
    conf_tmp = runData.perfSummary.confidence.lowOrHigh;
    confidence(iRun,:) = conf_tmp;
    confidence(iRun, conf_tmp == 0) = NaN; % no answer provided
    confidence(iRun, conf_tmp == 1) = 0; % low confidence
    confidence(iRun, conf_tmp == 2) = 1; % high confidence
    % extract RT
    RT_tmp = runData.perfSummary.onsets.choice - runData.perfSummary.onsets.dispChoiceOptions;
    RT(iRun,:) = RT_tmp;
    % reinitialize all counters
    [jR, jP] = deal(0);
    for iE = 1:(n_E_levels - 1)
        E_nm = ['E_level_',num2str(iE)];
        jE.(E_nm) = 0;
    end
    for iAbsMoney = 1:(n_R_levels - 1)
        absMoney_nm = ['absMoney_level_',num2str(iAbsMoney)];
        jAbsMoney.(absMoney_nm) = 0;
        deltaMoney_nm = ['deltaMoney_level_',num2str(iAbsMoney)];
        jDeltaMoney.(deltaMoney_nm) = 0;
        jSignedMoney.(['R_level_',num2str(iAbsMoney)]) = 0;
        jSignedMoney.(['P_level_',num2str(iAbsMoney)]) = 0;
    end
    
    for iTrial = 1:nTrials
        %% extract effort and money level of the non-default option for each
        % trial
        switch defaultSide(iTrial)
            case -1
                E_level(iTrial) = runData.choiceOptions.E.right(iTrial);
                absMoney_level(iTrial) = runData.choiceOptions.R.right(iTrial);
            case 1
                E_level(iTrial) = runData.choiceOptions.E.left(iTrial);
                absMoney_level(iTrial) = runData.choiceOptions.R.left(iTrial);
        end
        E_nm = ['E_level_',num2str(E_level(iTrial))];
        jE.(E_nm) = jE.(E_nm) + 1;
        absMoney_nm = ['absMoney_level_',num2str(absMoney_level(iTrial))];
        jAbsMoney.(absMoney_nm) = jAbsMoney.(absMoney_nm) + 1;
        switch R_or_P{iTrial}
            case 'R'
                jR = jR + 1;
                signedMoney_nm = ['R_level_',num2str(absMoney_level(iTrial))];
                deltaMoney_nm = ['deltaMoney_level_',num2str(absMoney_level(iTrial))];
            case 'P'
                jP = jP + 1;
                signedMoney_nm = ['P_level_',num2str(absMoney_level(iTrial))];
                deltaMoney_nm = ['deltaMoney_level_',num2str(n_R_levels-absMoney_level(iTrial))];
        end
        jSignedMoney.(signedMoney_nm) = jSignedMoney.(signedMoney_nm) + 1;
        jDeltaMoney.(deltaMoney_nm) = jDeltaMoney.(deltaMoney_nm) + 1;

        %% did the participant choose the default or the non-default option
        % for the current trial?
        if choiceSide(iTrial) == -defaultSide(iTrial) % choice non-default option
            % extract average selection of the default option
            choiceNonDefault(iRun, iTrial) = 1;

            % extract default choices per money level
            choiceNonDefault_perAbsMoneylevel.(absMoney_nm)(iRun, jAbsMoney.(absMoney_nm)) = 1;
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
            choiceNonDefault_perAbsMoneylevel.(absMoney_nm)(iRun, jAbsMoney.(absMoney_nm)) = 0;
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
        conf_perAbsMoneylevel.(absMoney_nm)(iRun, jAbsMoney.(absMoney_nm)) =  confidence(iRun, iTrial);
        conf_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = confidence(iRun,iTrial);
        conf_perDeltaMoneylevel.(deltaMoney_nm)(iRun, jDeltaMoney.(deltaMoney_nm)) =  confidence(iRun, iTrial);
        conf_perElevel.(E_nm)(iRun, jE.(E_nm)) = confidence(iRun,iTrial);
        
        %% level of reaction times (RT) depending on money/effort
        RT_perAbsMoneylevel.(absMoney_nm)(iRun, jAbsMoney.(absMoney_nm)) =  RT(iRun, iTrial);
        RT_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = RT(iRun,iTrial);
        RT_perDeltaMoneylevel.(deltaMoney_nm)(iRun, jDeltaMoney.(deltaMoney_nm)) =  RT(iRun, iTrial);
        RT_perElevel.(E_nm)(iRun, jE.(E_nm)) = RT(iRun,iTrial);
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
    end
    
    %% sum of failed choices to check
    n_failedChoices(iRun) = sum(choice(iRun,:) == 0,'omitnan');
    
    %% percentage choice of default per effort level
    for iE = 1:(n_E_levels - 1)
        choiceNonDefault_perElevel.(['avg_E_',num2str(iE)])(iRun) = 100*sum(choiceNonDefault_perElevel.(['E_level_',num2str(iE)])(iRun, :),'omitnan')/nRepeatsPerEffortLevel;
        conf_perElevel.(['avg_E_',num2str(iE)])(iRun) = 100*mean(conf_perElevel.(['E_level_',num2str(iE)])(iRun, :),2,'omitnan');
        RT_perElevel.(['avg_E_',num2str(iE)])(iRun) = mean(RT_perElevel.(['E_level_',num2str(iE)])(iRun, :),2,'omitnan');
    end

    %% percentage choice of default per money level
    for iAbsMoney = 1:(n_R_levels - 1)
        choiceNonDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(iRun) = 100*sum(choiceNonDefault_perAbsMoneylevel.(['absMoney_level_',num2str(iAbsMoney)])(iRun, :),'omitnan')/nRepeatsPerAbsMoneyLevel;
        choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(iRun) = 100*sum(choiceNonDefault_perDeltaMoneylevel.(['deltaMoney_level_',num2str(iAbsMoney)])(iRun, :),'omitnan')/nRepeatsPerAbsMoneyLevel;
        choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(iRun) = 100*sum(choiceNonDefault_perSignedMoneylevel.(['R_level_',num2str(iAbsMoney)])(iRun, :),'omitnan')/nRepeatsPerSignedMoneyLevel;
        choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(iRun) = 100*sum(choiceNonDefault_perSignedMoneylevel.(['P_level_',num2str(iAbsMoney)])(iRun, :),'omitnan')/nRepeatsPerSignedMoneyLevel;
        % same for confidence
        conf_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(iRun) = mean(100*conf_perAbsMoneylevel.(['absMoney_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
        conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(iRun) = mean(100*conf_perDeltaMoneylevel.(['deltaMoney_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
        conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(iRun) = mean(100*conf_perSignedMoneylevel.(['R_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
        conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(iRun) = mean(100*conf_perSignedMoneylevel.(['P_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
        % same for RT
        RT_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(iRun) = mean(RT_perAbsMoneylevel.(['absMoney_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
        RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(iRun) = mean(RT_perDeltaMoneylevel.(['deltaMoney_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
        RT_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(iRun) = mean(RT_perSignedMoneylevel.(['R_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
        RT_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(iRun) = mean(RT_perSignedMoneylevel.(['P_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
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
    % STD
    sd_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]) = std(choiceNonDefault_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]) = std(choiceNonDefault_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
    sd_conf.perElevel.(['Em_',num2str(iE)]) = std(conf_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_conf.perElevel.(['Ep_',num2str(iE)]) = std(conf_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
    sd_RT.perElevel.(['Em_',num2str(iE)]) = std(RT_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_RT.perElevel.(['Ep_',num2str(iE)]) = std(RT_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
end

% per task per money level
for iAbsMoney = 1:(n_R_levels - 1)
    avg_nonDefaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]) = mean(choiceNonDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Em_runs));
    avg_nonDefaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]) = mean(choiceNonDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Ep_runs));
    avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]) = mean(choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Em_runs));
    avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]) = mean(choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Ep_runs));
    avg_nonDefaultChoice.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = mean(choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = mean(choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    avg_nonDefaultChoice.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = mean(choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = mean(choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
    avg_conf.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]) = mean(conf_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Em_runs));
    avg_conf.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]) = mean(conf_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Ep_runs));
    avg_conf.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]) = mean(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Em_runs));
    avg_conf.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]) = mean(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Ep_runs));
    avg_conf.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = mean(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    avg_conf.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = mean(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    avg_conf.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = mean(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    avg_conf.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = mean(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
    avg_RT.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]) = mean(RT_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Em_runs));
    avg_RT.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]) = mean(RT_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Ep_runs));
    avg_RT.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]) = mean(RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Em_runs));
    avg_RT.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]) = mean(RT_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Ep_runs));
    avg_RT.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = mean(RT_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    avg_RT.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = mean(RT_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    avg_RT.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = mean(RT_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    avg_RT.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = mean(RT_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
    % STD
    sd_nonDefaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]) = std(choiceNonDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Em_runs));
    sd_nonDefaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]) = std(choiceNonDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Ep_runs));
    sd_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]) = std(choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Em_runs));
    sd_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]) = std(choiceNonDefault_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Ep_runs));
    sd_nonDefaultChoice.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = std(choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    sd_nonDefaultChoice.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = std(choiceNonDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    sd_nonDefaultChoice.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = std(choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    sd_nonDefaultChoice.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = std(choiceNonDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
    sd_conf.perAbsMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = std(conf_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Em_runs));
    sd_conf.perAbsMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = std(conf_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Ep_runs));
    sd_conf.perDeltaMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = std(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Em_runs));
    sd_conf.perDeltaMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = std(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Ep_runs));
    sd_conf.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    sd_conf.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    sd_conf.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    sd_conf.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
    sd_RT.perAbsMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = std(conf_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Em_runs));
    sd_RT.perAbsMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = std(conf_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Ep_runs));
    sd_RT.perDeltaMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = std(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Em_runs));
    sd_RT.perDeltaMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = std(conf_perDeltaMoneylevel.(['avg_deltaMoney_',num2str(iAbsMoney)])(Ep_runs));
    sd_RT.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    sd_RT.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    sd_RT.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    sd_RT.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
end

avg_nonDefaultChoice.Em_f_time = mean(choiceNonDefault_f_time(Em_runs,:),1,'omitnan');
avg_nonDefaultChoice.Ep_f_time = mean(choiceNonDefault_f_time(Ep_runs,:),1,'omitnan');
avg_conf.Em_f_time = mean(conf_f_time(Em_runs,:),1,'omitnan');
avg_conf.Ep_f_time = mean(conf_f_time(Ep_runs,:),1,'omitnan');
avg_RT.Em_f_time = mean(RT_f_time(Em_runs,:),1,'omitnan');
avg_RT.Ep_f_time = mean(RT_f_time(Ep_runs,:),1,'omitnan');

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
    plot(1:n_bins, 50*ones(1,n_bins),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iRun = 1:nRuns
        % different colour for mental and physical effort
        if ismember(iRun, Em_runs)
            lColour = 'g';
        else
            lColour = 'r';
        end
        plot(1:n_bins, choiceNonDefault_f_time(iRun, :)*100,...
            'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    end
    ylim([0 100]);
    xlim([0 n_bins+1]);
    xlabel('trial bins');
    ylabel('Choice non-default option (%)');
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
    
    
    %% check choices = f(|money level|)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iAbsMoney = 1:(n_R_levels - 1)
        bar(iAbsMoney-bDist, avg_nonDefaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iAbsMoney+bDist, avg_nonDefaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iAbsMoney-bDist, avg_nonDefaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]), sd_nonDefaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]),'k')
        errorbar(iAbsMoney+bDist, avg_nonDefaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]), sd_nonDefaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]),'k')
    end
    ylim([0 100]);
    ylabel('Choice non-default option (%)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_R_levels]);
    xlabel('Money level');
    legend_size(pSize);
    ylim([0 100]);
    ylabel('Choice non-default option (%)');
    legend_size(pSize);
    
    %% check choices = f(|high money level - default money level|)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iDeltaMoney = 1:(n_R_levels - 1)
        bar(iDeltaMoney-bDist, avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iDeltaMoney)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iDeltaMoney+bDist, avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iDeltaMoney)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iDeltaMoney-bDist, avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iDeltaMoney)]), sd_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iDeltaMoney)]),'k')
        errorbar(iDeltaMoney+bDist, avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iDeltaMoney)]), sd_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iDeltaMoney)]),'k')
    end
    ylim([0 100]);
    ylabel('Choice non-default option (%)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_R_levels]);
    xlabel('|Δ money| level');
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
        bar(iMoney-bDist, avg_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iMoney+bDist, avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iMoney-bDist, avg_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), sd_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),'k')
        errorbar(iMoney+bDist, avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), sd_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),'k')
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
    legend_size(pSize);
    
    %% check choices = f(E level)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iE = 1:(n_E_levels - 1)
        bar(iE-bDist, avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iE+bDist, avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iE-bDist, avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]), sd_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]),'k')
        errorbar(iE+bDist, avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]), sd_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]),'k')
    end
    ylim([0 100]);
    ylabel('Choice non-default option (%)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_E_levels]);
    xlabel('Effort level');
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
        bar(iMoney-bDist, avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iMoney+bDist, avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iMoney-bDist, avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), sd_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),'k')
        errorbar(iMoney+bDist, avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), sd_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),'k')
    end
    xticks([-3:-(1), 1:3]);
    xticklabels({'-3','-2','-1','1','2','3'});
    xlim([-n_R_levels n_R_levels]);
    ylim([0 100]);
    ylabel('Level of confidence');
    xlabel('Money level');
    legend_size(pSize);
    
    %% check confidence = f(E levels)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iE = 1:(n_E_levels - 1)
        bar(iE-bDist, avg_conf.perElevel.(['Em_',num2str(iE)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iE+bDist, avg_conf.perElevel.(['Ep_',num2str(iE)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iE-bDist, avg_conf.perElevel.(['Em_',num2str(iE)]), sd_conf.perElevel.(['Em_',num2str(iE)]),'k')
        errorbar(iE+bDist, avg_conf.perElevel.(['Ep_',num2str(iE)]), sd_conf.perElevel.(['Ep_',num2str(iE)]),'k')
    end
    ylim([0 100]);
    ylabel('Confidence');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_E_levels]);
    xlabel('Effort level');
    legend_size(pSize);
    
    
    ylabel('Level of confidence');
    xlabel('Effort level');
    legend_size(pSize);
    
    %% check confidence = f(time)
    fig;
    % mark the 50% trait
    plot(1:n_bins, 50*ones(1,n_bins),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iRun = 1:nRuns
        % different colour for mental and physical effort
        if ismember(iRun, Em_runs)
            lColour = 'g';
        else
            lColour = 'r';
        end
        plot(1:n_bins, conf_f_time(iRun, :)*100,...
            'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    end
    ylim([0 100]);
    xlim([0 n_bins+1]);
    xlabel('trial bins');
    ylabel('Confidence (%)');
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
        bar(iMoney-bDist, avg_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iMoney+bDist, avg_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iMoney-bDist, avg_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), sd_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),'k')
        errorbar(iMoney+bDist, avg_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), sd_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),'k')
    end
    xticks([-3:-(1), 1:3]);
    xticklabels({'-3','-2','-1','1','2','3'});
    xlim([-n_R_levels n_R_levels]);
    ylim([0 5]);
    ylabel('RT (s)');
    xlabel('Money level');
    legend_size(pSize);
    
    %% check RT = f(E levels)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iE = 1:(n_E_levels - 1)
        bar(iE-bDist, avg_RT.perElevel.(['Em_',num2str(iE)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iE+bDist, avg_RT.perElevel.(['Ep_',num2str(iE)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iE-bDist, avg_RT.perElevel.(['Em_',num2str(iE)]), sd_RT.perElevel.(['Em_',num2str(iE)]),'k')
        errorbar(iE+bDist, avg_RT.perElevel.(['Ep_',num2str(iE)]), sd_RT.perElevel.(['Ep_',num2str(iE)]),'k')
    end
    ylim([0 5]);
    ylabel('RT (s)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_E_levels]);
    xlabel('Effort level');
    legend_size(pSize);
    
    %% check RT = f(time)
    fig;
    % mark the 50% trait
    plot(1:n_bins, 50*ones(1,n_bins),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iRun = 1:nRuns
        % different colour for mental and physical effort
        if ismember(iRun, Em_runs)
            lColour = 'g';
        else
            lColour = 'r';
        end
        plot(1:n_bins, RT_f_time(iRun, :),...
            'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    end
    ylim([0 5]);
    xlim([0 n_bins+1]);
    xlabel('trial bins');
    ylabel('RT (s)');
    legend_size(pSize);
    
end % figure display

end % function display