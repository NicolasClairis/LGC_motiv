%% check performance = f(fatigue, R/P, level of R, level of E)

subid = '202';
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
for iE = 0:(n_E_levels - 1)
    E_nm = ['E_level_',num2str(iE)];
    jE.(E_nm) = 0;
    if iE > 0
        choiceDefault_perElevel.(E_nm) = NaN(nRuns, nRepeatsPerEffortLevel);
    end
    for iRun = 1:nRuns
        sess_nm = ['run_',num2str(iRun)];
        performance_perElevel.(['E_level_',num2str(iE)]).(sess_nm) = [];
    end
end

choiceDefault_f_fatigue = NaN(nRuns, n_bins);
for iRun = 1:nRuns
    run_nm = ['run_',num2str(iRun)];
    filenm = ls(['CID',subid,'_session',num2str(iRun),'_*_task_messyAllStuff.mat']);
    runData = load(filenm);
    defaultSide = runData.choiceOptions.default_LR;
    R_or_P = runData.choiceOptions.R_or_P;
    choice(iRun,:) = runData.perfSummary.choice;
    choiceSide = choice(iRun,:);
    choiceSide(choice(iRun,:) >= 1) = 1;
    choiceSide(choice(iRun,:) <= -1) = -1;
    [jR, jP] = deal(0);
    for iTrial = 1:nTrials
%         performance(iRun, iTrial) = runData.perfSummary.perfSummary.performance{1,iTrial};
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
            E_nm = 'E_level_0';
            jE.(E_nm) = jE.(E_nm) + 1;
            choiceDefault_perElevel.(E_nm)(iRun, jE.(E_nm)) = 1;
%             performance_perElevel.(E_nm).(run_nm) = [performance_perElevel.(E_nm).(run_nm), performance(iRun, iTrial)];
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
%             E_nm = ['E_level_',num2str(E_level(iTrial))];
%             jE.(E_nm) = jE.(E_nm) + 1;
%             choiceDefault_perElevel.(E_nm)(iRun, jE.(E_nm)) = 1;
%             performance_perElevel.(E_nm).(run_nm) = [performance_perElevel.(E_nm).(run_nm), performance(iRun, iTrial)];
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
end

%% extract task type order
lastTaskType = getfield(load(['training_data_CID',subid,'.mat'],'p_or_m'),'p_or_m');
switch lastTaskType
    case 'p'
        Em_runs = [1,3];
        Ep_runs = [2,4];
    case 'm'
        Em_runs = [2,4];
        Ep_runs = [1,3];
end

%% figures
if figDisp == 1
    % figure parameters
    lWidth = 3;
    pSize = 30;
    
    % performance = f(fatigue)
    fig;
    
    ylim([0 100]);
    ylabel('Performance (%)');
    legend_size(pSize);
    
    % performance = f(R/P)
    fig;
    
    ylim([0 100]);
    ylabel('Performance (%)');
    legend_size(pSize);
    
    % performance = f(R level)
    fig;
    
    ylim([0 100]);
    ylabel('Performance (%)');
    legend_size(pSize);
    
    % performance = f(E level)
    fig;
    
    ylim([0 100]);
    ylabel('Performance (%)');
    legend_size(pSize);
    
    % performance = f(R/P level)
    fig;
    
    ylim([0 100]);
    ylabel('Performance (%)');
    legend_size(pSize);
    
end % figure display