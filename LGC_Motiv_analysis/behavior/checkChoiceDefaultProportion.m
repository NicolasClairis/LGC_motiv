subid = '202';
nRuns = 4;
nTrials = 54;
n_bins = 6;
nBinTrials = nTrials/n_bins;
[percentageChoiceDefault,...
    percentage_R_choiceDefault,...
    percentage_P_choiceDefault] = deal(NaN(1,nRuns));
[choiceDefault, choice] = deal( NaN(nRuns,nTrials));
[choiceDefault_R,...
    choiceDefault_P] = deal( NaN(nRuns,nTrials/2));
choiceDefault_f_fatigue = NaN(nRuns, n_bins);
for iRun = 1:nRuns
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
        end
    end
    %% ratio of chosen the default option per run
    percentageChoiceDefault(iRun) = sum(choiceDefault(iRun,:))/nTrials;
    percentage_R_choiceDefault(iRun) = sum(choiceDefault_R(iRun,:))/(nTrials/2);
    percentage_P_choiceDefault(iRun) = sum(choiceDefault_P(iRun,:))/(nTrials/2);
    
    %% default choice = f(fatigue)
    for iBin = 1:n_bins
        trial_idx = (1:nBinTrials) + nBinTrials*(iBin - 1);
        choiceDefault_f_fatigue(iRun, iBin) = mean( choiceDefault(iRun,trial_idx) );
    end
end

%% ratio of choosing default per task
lastTaskType = getfield(load(['training_data_CID',subid,'.mat'],'p_or_m'),'p_or_m');
switch lastTaskType
    case 'p'
        Em_runs = [1,3];
        Ep_runs = [2,4];
    case 'm'
        Em_runs = [2,4];
        Ep_runs = [1,3];
end

Em_avg_defaultChoice = mean(percentageChoiceDefault(Em_runs));
Ep_avg_defaultChoice = mean(percentageChoiceDefault(Ep_runs));

%% check with fatigue
fig;
lWidth = 3;
pSize = 30;
for iRun = 1:nRuns
    % different colour for mental and physical effort
    if ismember(iRun, Em_runs)
        lColour = 'g';
    else
        lColour = 'r';
    end
    plot(1:n_bins, choiceDefault_f_fatigue(iRun, :)*100,...
        'LineWidth',lWidth, 'Color',lColour,'LineStyle','--');
    hold on;
end
ylim([0 100]);
xlim([0 n_bins+1]);
xlabel('Bins');
ylabel('Choice of the default option (%)');
legend_size(pSize);

%% check proportion of choice failures to see if time ok

