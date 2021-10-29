%% check choices of the default option in function of the other parameters of the task
% (fatigue, difficulty, physical/mental, reward/punishment, etc.)

subid = '203';
figDisp = 1;
nRuns = 4;
n_R_levels = 4;
n_E_levels = 4;
nTrials = 54;
n_bins = 6;
nBinTrials = nTrials/n_bins;
[percentageChoiceDefault,...
    percentage_R_choiceDefault,...
    percentage_P_choiceDefault,...
    n_failedChoices] = deal(NaN(1,nRuns));
[choiceDefault, choice, performance] = deal( NaN(nRuns,nTrials));
[choiceDefault_R,...
    choiceDefault_P] = deal( NaN(nRuns,nTrials/2));

nRepeatsPerEffortLevel = nTrials/(n_E_levels - 1);
for iE = 1:(n_E_levels - 1)
    E_nm = ['E_level_',num2str(iE)];
    jE.(E_nm) = 0;
    choiceDefault_perElevel.(E_nm) = NaN(nRuns, nRepeatsPerEffortLevel);
    choiceDefault_perElevel.(['avg_E_',num2str(iE)]) = NaN(1,nRuns);
end

choiceDefault_f_fatigue = NaN(nRuns, n_bins);
for iRun = 1:nRuns
    run_nm = ['run_',num2str(iRun)];
    filenm = ls(['CID',subid,'_session',num2str(iRun),'_*_task_messyAllStuff.mat']);
    runData = load(filenm);
    defaultSide = runData.choiceOptions.default_LR;
    E_level = NaN(1,nTrials);
    R_or_P = runData.choiceOptions.R_or_P;
    choice(iRun,:) = runData.perfSummary.choice;
    choiceSide = choice(iRun,:);
    choiceSide(choice(iRun,:) >= 1) = 1;
    choiceSide(choice(iRun,:) <= -1) = -1;
    % reinitialize all counters
    [jR, jP] = deal(0);
    for iE = 1:(n_E_levels - 1)
        E_nm = ['E_level_',num2str(iE)];
        jE.(E_nm) = 0;
    end
    
    for iTrial = 1:nTrials
        switch defaultSide(iTrial)
            case -1
                E_level(iTrial) = runData.choiceOptions.E.right(iTrial);
            case 1
                E_level(iTrial) = runData.choiceOptions.E.left(iTrial);
        end
        E_nm = ['E_level_',num2str(E_level(iTrial))];
        if choiceSide(iTrial) == defaultSide(iTrial)
            choiceDefault(iRun, iTrial) = 1;
            switch R_or_P{iTrial}
                case 'R'
                    jR = jR + 1;
                    choiceDefault_R(iRun,jR) = 1;
                case 'P'
                    jP = jP + 1;
                    choiceDefault_P(iRun,jP) = 1;
            end
            
            % extract choice and perf per effort level
            jE.(E_nm) = jE.(E_nm) + 1;
            choiceDefault_perElevel.(E_nm)(iRun, jE.(E_nm)) = 1;
        else
            choiceDefault(iRun, iTrial) = 0;
            switch R_or_P{iTrial}
                case 'R'
                    jR = jR + 1;
                    choiceDefault_R(iRun,jR) = 0;
                case 'P'
                    jP = jP + 1;
                    choiceDefault_P(iRun,jP) = 0;
            end
            
            % extract choice and perf per effort level
            jE.(E_nm) = jE.(E_nm) + 1;
            choiceDefault_perElevel.(E_nm)(iRun, jE.(E_nm)) = 0;
        end
    end % trial loop
    
    %% ratio of chosen the default option per run
    percentageChoiceDefault(iRun) = sum(choiceDefault(iRun,:))/nTrials;
    percentage_R_choiceDefault(iRun) = sum(choiceDefault_R(iRun,:))/(nTrials/2);
    percentage_P_choiceDefault(iRun) = sum(choiceDefault_P(iRun,:))/(nTrials/2);
    
    %% default choice = f(fatigue)
    for iBin = 1:n_bins
        trial_idx = (1:nBinTrials) + nBinTrials*(iBin - 1);
        choiceDefault_f_fatigue(iRun, iBin) = mean( choiceDefault(iRun,trial_idx) );
    end
    
    %% sum of failed choices to check
    n_failedChoices(iRun) = sum(choice(iRun,:) == 0);
    
    %% percentage choice of default per effort level
    for iE = 1:(n_E_levels - 1)
        choiceDefault_perElevel.(['avg_E_',num2str(iE)])(iRun) = 100*sum(choiceDefault_perElevel.(['E_level_',num2str(iE)])(iRun, :))/nRepeatsPerEffortLevel;
    end
end

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
Em_avg_defaultChoice = mean(percentageChoiceDefault(Em_runs));
Ep_avg_defaultChoice = mean(percentageChoiceDefault(Ep_runs));

% per task per effort level
for iE = 1:(n_E_levels - 1)
    avg_choiceDefault_perElevel.(['Em_',num2str(iE)]) = mean(choiceDefault_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    avg_choiceDefault_perElevel.(['Ep_',num2str(iE)]) = mean(choiceDefault_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
    % STD
    sd_choiceDefault_perElevel.(['Em_',num2str(iE)]) = std(choiceDefault_perElevel.(['avg_E_',num2str(iE)])(Em_runs));
    sd_choiceDefault_perElevel.(['Ep_',num2str(iE)]) = std(choiceDefault_perElevel.(['avg_E_',num2str(iE)])(Ep_runs));
end

%% figures
if figDisp == 1
    % figure parameters
    lWidth = 3;
    pSize = 30;
    bWidth = 0.4;
    bDist = 0.2;
    %% check choices = f(fatigue, R/P, level of R, level of E)
    
    % check choices = f(fatigue)
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
        plot(1:n_bins, choiceDefault_f_fatigue(iRun, :)*100,...
            'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    end
    ylim([0 100]);
    xlim([0 n_bins+1]);
    xlabel('trial bins');
    ylabel('Choice = default option (%)');
    legend_size(pSize);
    
    % check choices = f(R/P)
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
    
    
    % check choices = f(money level)
    fig;
    % mark the 50% trait
    plot(1:n_bins, 50*ones(1,n_bins),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    
    ylim([0 100]);
    ylabel('Choice = default option (%)');
    legend_size(pSize);
    
    
    % check choices = f(E level)
    fig;
    % mark the 50% trait
    plot(0:n_bins, 50*ones(1,n_bins+1),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    for iE = 1:(n_E_levels - 1)
        bar(iE-bDist, avg_choiceDefault_perElevel.(['Em_',num2str(iE)]), 'FaceColor','g','BarWidth',bWidth);
        bar(iE+bDist, avg_choiceDefault_perElevel.(['Ep_',num2str(iE)]), 'FaceColor','b','BarWidth',bWidth);
        errorbar(iE-bDist, avg_choiceDefault_perElevel.(['Em_',num2str(iE)]), sd_choiceDefault_perElevel.(['Em_',num2str(iE)]),'k')
        errorbar(iE+bDist, avg_choiceDefault_perElevel.(['Ep_',num2str(iE)]), sd_choiceDefault_perElevel.(['Ep_',num2str(iE)]),'k')
    end
    ylim([0 100]);
    ylabel('Choice = default option (%)');
    xticks(1:3);
    xticklabels({'1','2','3'});
    xlim([0 n_E_levels]);
    xlabel('Effort level');
    legend_size(pSize);
    
    % check choices = f(R-P levels)
    fig;
    % mark the 50% trait
    plot(1:n_bins, 50*ones(1,n_bins),...
        'LineWidth',2,'Color','k','LineStyle',':');
    hold on;
    
    ylim([0 100]);
    ylabel('Choice = default option (%)');
    legend_size(pSize);
    
end % figure display