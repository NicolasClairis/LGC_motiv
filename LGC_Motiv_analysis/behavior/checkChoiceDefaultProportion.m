% function[avg_defaultChoice, sd_defaultChoice] = checkChoiceDefaultProportion(subid, figDisp, n_bins)
%[avg_defaultChoice, sd_defaultChoice] = checkChoiceDefaultProportion(subid, figDisp, n_bins)
% check choices of the default option in function of the other parameters of the task
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
% avg_defaultChoice: structure with average choice of the default option in
% function of the other task parameters
%
% sd_defaultChoice: same as avg_defaultChoice but for standard deviation


%% initialize all variables of interest
if ~exist('subid','var') || isempty(subid)
    subid = '204';
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
[percentageChoiceDefault,...
    percentage_R_choiceDefault,...
    percentage_P_choiceDefault,...
    n_failedChoices] = deal(NaN(1,nRuns));
[choiceDefault, choice, confidence] = deal( NaN(nRuns,nTrials));
[choiceDefault_R,...
    choiceDefault_P] = deal( NaN(nRuns,nTrials/2));

% prepare analysis per effort level
nRepeatsPerEffortLevel = nTrials/(n_E_levels - 1);
for iE = 1:(n_E_levels - 1)
    E_nm = ['E_level_',num2str(iE)];
    [choiceDefault_perElevel.(E_nm),...
        conf_perElevel.(E_nm)] = deal(NaN(nRuns, nRepeatsPerEffortLevel));
    [choiceDefault_perElevel.(['avg_E_',num2str(iE)]),...
        conf_perElevel.(['avg_E_',num2str(iE)])] = deal(NaN(1,nRuns));
end
% prepare analysis per money level for unsigned (R+P) and signed (P/R)
% money levels
% 1) unsigned money levels
nRepeatsPerAbsMoneyLevel = nTrials/(n_R_levels - 1);
for iAbsMoney = 1:(n_R_levels - 1)
    absMoney_nm = ['absMoney_level_',num2str(iAbsMoney)];
    [choiceDefault_perAbsMoneylevel.(absMoney_nm),...
        conf_perAbsMoneylevel.(absMoney_nm)] = deal(NaN(nRuns, nRepeatsPerAbsMoneyLevel));
    [choiceDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)]),...
        conf_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])] = deal(NaN(1,nRuns));
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
    [choiceDefault_perSignedMoneylevel.(signedMoney_nm),...
        conf_perSignedMoneylevel.(signedMoney_nm)] = deal(NaN(nRuns, nRepeatsPerSignedMoneyLevel));
    [choiceDefault_perSignedMoneylevel.(['avg_',signedMoney_nm]),...
        conf_perSignedMoneylevel.(signedMoney_nm)] = deal(NaN(1,nRuns));
end

%% extract data trial by trial
[choiceDefault_f_time,...
    conf_f_time] = deal(NaN(nRuns, n_bins));
for iRun = 1:nRuns
    filenm = ls(['CID',subid,'_session',num2str(iRun),'_*_task_messyAllStuff.mat']);
    runData = load(filenm);
    defaultSide = runData.choiceOptions.default_LR;
    [absMoney_level, E_level] = deal(NaN(1,nTrials));
    R_or_P = runData.choiceOptions.R_or_P;
    choice(iRun,:) = runData.perfSummary.choice;
    choiceSide = choice(iRun,:);
    choiceSide(choice(iRun,:) >= 1) = 1;
    choiceSide(choice(iRun,:) <= -1) = -1;
    conf_tmp = runData.perfSummary.confidence.lowOrHigh;
    confidence(iRun,:) = conf_tmp;
    confidence(iRun, conf_tmp == 0) = NaN; % no answer provided
    confidence(iRun, conf_tmp == 1) = 0; % low confidence
    confidence(iRun, conf_tmp == 2) = 1; % high confidence
    % reinitialize all counters
    [jR, jP] = deal(0);
    for iE = 1:(n_E_levels - 1)
        E_nm = ['E_level_',num2str(iE)];
        jE.(E_nm) = 0;
    end
    for iAbsMoney = 1:(n_R_levels - 1)
        absMoney_nm = ['absMoney_level_',num2str(iAbsMoney)];
        jAbsMoney.(absMoney_nm) = 0;
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
                jSignedMoney.(signedMoney_nm) = jSignedMoney.(signedMoney_nm) + 1;
            case 'P'
                jP = jP + 1;
                signedMoney_nm = ['P_level_',num2str(absMoney_level(iTrial))];
                jSignedMoney.(signedMoney_nm) = jSignedMoney.(signedMoney_nm) + 1;
        end

        %% did the participant choose the default or the non-default option
        % for the current trial?
        if choiceSide(iTrial) == defaultSide(iTrial)
            % extract average selection of the default option
            choiceDefault(iRun, iTrial) = 1;

            % extract default choices per money level
            choiceDefault_perAbsMoneylevel.(absMoney_nm)(iRun, jAbsMoney.(absMoney_nm)) = 1;
            choiceDefault_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = 1;

            % extract average reward (or respectively punishment) default choices
            switch R_or_P{iTrial}
                case 'R'
                    choiceDefault_R(iRun,jR) = 1;
                case 'P'
                    choiceDefault_P(iRun,jP) = 1;
            end
            
            % extract choice and perf per effort level
            choiceDefault_perElevel.(E_nm)(iRun, jE.(E_nm)) = 1;
        else
            % extract average selection of the default option
            choiceDefault(iRun, iTrial) = 0;

            % extract default choices per money level
            choiceDefault_perAbsMoneylevel.(absMoney_nm)(iRun, jAbsMoney.(absMoney_nm)) = 0;
            choiceDefault_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = 0;

            % extract average reward (or respectively punishment) default choices
            switch R_or_P{iTrial}
                case 'R'
                    choiceDefault_R(iRun,jR) = 0;
                case 'P'
                    choiceDefault_P(iRun,jP) = 0;
            end
            
            % extract choice and perf per effort level
            choiceDefault_perElevel.(E_nm)(iRun, jE.(E_nm)) = 0;
        end
        
        %% level of confidence depending on money/effort
        conf_perSignedMoneylevel.(signedMoney_nm)(iRun, jSignedMoney.(signedMoney_nm)) = confidence(iRun,iTrial);
        conf_perElevel.(E_nm)(iRun, jE.(E_nm)) = confidence(iRun,iTrial);
    end % trial loop
    
    %% ratio of chosen the default option per run
    percentageChoiceDefault(iRun) = sum(choiceDefault(iRun,:))/nTrials;
    percentage_R_choiceDefault(iRun) = sum(choiceDefault_R(iRun,:))/(nTrials/2);
    percentage_P_choiceDefault(iRun) = sum(choiceDefault_P(iRun,:))/(nTrials/2);
    
    %% default choice/confidence = f(time)
    for iBin = 1:n_bins
        trial_idx = (1:nBinTrials) + nBinTrials*(iBin - 1);
        choiceDefault_f_time(iRun, iBin) = mean( choiceDefault(iRun,trial_idx) );
        conf_f_time(iRun, iBin) = mean( confidence(iRun,trial_idx), 2,'omitnan' );
    end
    
    %% sum of failed choices to check
    n_failedChoices(iRun) = sum(choice(iRun,:) == 0);
    
    %% percentage choice of default per effort level
    for iE = 1:(n_E_levels - 1)
        choiceDefault_perElevel.(['avg_E_',num2str(iE)])(iRun) = 100*sum(choiceDefault_perElevel.(['E_level_',num2str(iE)])(iRun, :))/nRepeatsPerEffortLevel;
        conf_perElevel.(['avg_E_',num2str(iE)])(iRun) = 100*mean(conf_perElevel.(['E_level_',num2str(iE)])(iRun, :),2,'omitnan');
    end

    %% percentage choice of default per money level
    for iAbsMoney = 1:(n_R_levels - 1)
        choiceDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(iRun) = 100*sum(choiceDefault_perAbsMoneylevel.(['absMoney_level_',num2str(iAbsMoney)])(iRun, :))/nRepeatsPerAbsMoneyLevel;
        choiceDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(iRun) = 100*sum(choiceDefault_perSignedMoneylevel.(['R_level_',num2str(iAbsMoney)])(iRun, :))/nRepeatsPerSignedMoneyLevel;
        choiceDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(iRun) = 100*sum(choiceDefault_perSignedMoneylevel.(['P_level_',num2str(iAbsMoney)])(iRun, :))/nRepeatsPerSignedMoneyLevel;
        % same for confidence
        conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(iRun) = mean(100*conf_perSignedMoneylevel.(['R_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
        conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(iRun) = mean(100*conf_perSignedMoneylevel.(['P_level_',num2str(iAbsMoney)])(iRun, :),2,'omitnan');
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

%% ratio of choosing default per task
avg_defaultChoice.Em = mean(percentageChoiceDefault(Em_runs));
avg_defaultChoice.Ep = mean(percentageChoiceDefault(Ep_runs));

% per task per effort level
for iE = 1:(n_E_levels - 1)
    avg_defaultChoice.perElevel.(['Em_',num2str(iE)]) = mean(choiceDefault_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    avg_defaultChoice.perElevel.(['Ep_',num2str(iE)]) = mean(choiceDefault_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
    avg_conf.perElevel.(['Em_',num2str(iE)]) = mean(conf_perElevel.(['avg_E_',num2str(iE)])(Em_runs),'omitnan');
    avg_conf.perElevel.(['Ep_',num2str(iE)]) = mean(conf_perElevel.(['avg_E_',num2str(iE)])(Ep_runs),'omitnan');
    % STD
    sd_defaultChoice.perElevel.(['Em_',num2str(iE)]) = std(choiceDefault_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_defaultChoice.perElevel.(['Ep_',num2str(iE)]) = std(choiceDefault_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
    sd_conf.perElevel.(['Em_',num2str(iE)]) = std(conf_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_conf.perElevel.(['Ep_',num2str(iE)]) = std(conf_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
end

% per task per money level
for iAbsMoney = 1:(n_R_levels - 1)
    avg_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]) = mean(choiceDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Em_runs));
    avg_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]) = mean(choiceDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Ep_runs));
    avg_defaultChoice.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = mean(choiceDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    avg_defaultChoice.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = mean(choiceDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    avg_defaultChoice.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = mean(choiceDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    avg_defaultChoice.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = mean(choiceDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
    avg_conf.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = mean(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    avg_conf.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = mean(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    avg_conf.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = mean(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    avg_conf.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = mean(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
    % STD
    sd_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]) = std(choiceDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Em_runs));
    sd_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]) = std(choiceDefault_perAbsMoneylevel.(['avg_absMoney_',num2str(iAbsMoney)])(Ep_runs));
    sd_defaultChoice.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = std(choiceDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    sd_defaultChoice.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = std(choiceDefault_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    sd_defaultChoice.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = std(choiceDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    sd_defaultChoice.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = std(choiceDefault_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
    sd_conf.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Em_runs));
    sd_conf.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_R_',num2str(iAbsMoney)])(Ep_runs));
    sd_conf.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Em_runs));
    sd_conf.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = std(conf_perSignedMoneylevel.(['avg_P_',num2str(iAbsMoney)])(Ep_runs));
end

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
        plot(1:n_bins, choiceDefault_f_time(iRun, :)*100,...
            'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    end
    ylim([0 100]);
    xlim([0 n_bins+1]);
    xlabel('trial bins');
    ylabel('Choice = default option (%)');
    legend_size(pSize);
    
    %% check choices = f(R/P)
    fig;
    % mark the 50% trait
    plot(0:3, 50*ones(1,4),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    bar(1, mean(percentage_R_choiceDefault.*100),'FaceColor',[0 153/255 1]);
    errorbar(1, mean(percentage_R_choiceDefault.*100), std(percentage_R_choiceDefault.*100),'k');
    bar(2, mean(percentage_P_choiceDefault.*100),'FaceColor',[1 153/255 0]);
    errorbar(2, mean(percentage_P_choiceDefault.*100), std(percentage_P_choiceDefault.*100),'k');
    xticks(1:2);
    xticklabels({'R','P'});
    ylim([0 100]);
    xlim([0 3]);
    ylabel('Choice = default option (%)');
    legend_size(pSize);
    
    
    %% check choices = f(money level)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iAbsMoney = 1:(n_R_levels - 1)
        bar(iAbsMoney-bDist, avg_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iAbsMoney+bDist, avg_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iAbsMoney-bDist, avg_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]), sd_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]),'k')
        errorbar(iAbsMoney+bDist, avg_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]), sd_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]),'k')
    end
    ylim([0 100]);
    ylabel('Choice = default option (%)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_R_levels]);
    xlabel('Money level');
    legend_size(pSize);
    ylim([0 100]);
    ylabel('Choice = default option (%)');
    legend_size(pSize);
    
    
    %% check choices = f(E level)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iE = 1:(n_E_levels - 1)
        bar(iE-bDist, avg_defaultChoice.perElevel.(['Em_',num2str(iE)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iE+bDist, avg_defaultChoice.perElevel.(['Ep_',num2str(iE)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iE-bDist, avg_defaultChoice.perElevel.(['Em_',num2str(iE)]), sd_defaultChoice.perElevel.(['Em_',num2str(iE)]),'k')
        errorbar(iE+bDist, avg_defaultChoice.perElevel.(['Ep_',num2str(iE)]), sd_defaultChoice.perElevel.(['Ep_',num2str(iE)]),'k')
    end
    ylim([0 100]);
    ylabel('Choice = default option (%)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_E_levels]);
    xlabel('Effort level');
    legend_size(pSize);
    
    %% check choices = f(R-P levels)
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
        bar(iMoney-bDist, avg_defaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iMoney+bDist, avg_defaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iMoney-bDist, avg_defaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), sd_defaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),'k')
        errorbar(iMoney+bDist, avg_defaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), sd_defaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),'k')
    end
    ylim([0 100]);
    ylabel('Choice = default option (%)');
    xticks([-3:-(1), 1:3]);
    xticklabels({'-3','-2','-1','1','2','3'});
    xlim([-n_R_levels n_R_levels]);
    xlabel('Money level');
    legend_size(pSize);
    ylim([0 100]);
    ylabel('Choice = default option (%)');
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
    
end % figure display

% end % function display