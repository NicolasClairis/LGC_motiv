%% check proportion of choices where people selected the highest effort option


taskTypes = {'mental','physical'};
nTasks = length(taskTypes);
subjectNames = {'MC','AC','HC','AR'};
NS = length(subjectNames);
nSessions = 2;
n_trials = 44;
n_RP_trials = n_trials/2;

[choiceHighRPerTrial.mental, choiceLowPPerTrial.mental,...
    choiceHighRPerTrial.mental, choiceLowPPerTrial.mental] = deal(NaN(n_RP_trials, nSessions, NS));
[choiceHighEPerTrial.mental, choiceHighEPerTrial.physical] = deal(NaN(n_trials, nSessions, NS));
[choiceHighR.mental, choiceLowP.mental, choiceHighE.mental,...
    choiceHighR.physical, choiceLowP.physical, choiceHighE.physical] = deal(NaN(nSessions, NS));
[choiceHighR_total.mental, choiceLowP_total.mental, choiceHighE_total.mental,...
    choiceHighR_total.physical, choiceLowP_total.physical, choiceHighE_total.physical] = deal(NaN(1, NS));
for iS = 1:NS
    subject_nm =  subjectNames{iS};
    for iTask = 1:nTasks
        taskType = taskTypes{iTask};
        for iSession = 1:nSessions
            session_summary_tmp = getfield(load(['pilot_data',subject_nm,'_sub_',num2str(iS),'_',taskType,'_session',num2str(iSession),'_behavioral_tmp.mat'],'summary'),'summary');
            
            % for each trial check if the option chosen was the highest or
            % not
            jR_Trial = 0;
            jP_Trial = 0;
            for iTrial = 1:n_trials
                
                R_maxOfTheTrial = nanmax([session_summary_tmp.choiceOptions.R.left(iTrial), session_summary_tmp.choiceOptions.R.right(iTrial)]);
                R_minOfTheTrial = nanmin([session_summary_tmp.choiceOptions.R.left(iTrial), session_summary_tmp.choiceOptions.R.right(iTrial)]);
                E_maxOfTheTrial = nanmax([session_summary_tmp.choiceOptions.E.left(iTrial), session_summary_tmp.choiceOptions.E.right(iTrial)]);
                switch session_summary_tmp.choiceOptions.R_or_P{iTrial}
                    case 'R'
                        jR_Trial = jR_Trial + 1;
                        if session_summary_tmp.R_chosen(iTrial) == R_maxOfTheTrial
                            choiceHighRPerTrial.(taskType)(jR_Trial,iSession,iS) = 1;
                        else
                            choiceHighRPerTrial.(taskType)(jR_Trial,iSession,iS) = 0;
                        end
                    case 'P'
                        jP_Trial = jP_Trial + 1;
                        if session_summary_tmp.R_chosen(iTrial) == R_minOfTheTrial
                            choiceLowPPerTrial.(taskType)(jP_Trial,iSession,iS) = 1;
                        else
                            choiceLowPPerTrial.(taskType)(jP_Trial,iSession,iS) = 0;
                        end
                end
                
                if session_summary_tmp.E_chosen(iTrial) == E_maxOfTheTrial
                    choiceHighEPerTrial.(taskType)(iTrial,iSession,iS) = 1;
                else
                    choiceHighEPerTrial.(taskType)(iTrial,iSession,iS) = 0;
                end
            end % trial loop
            
            % measure proportion of choice of high effort
            choiceHighR.(taskType)(iSession, iS) = (sum(choiceHighRPerTrial.(taskType)(:,iSession,iS))/n_RP_trials)*100;
            choiceLowP.(taskType)(iSession, iS) = (sum(choiceLowPPerTrial.(taskType)(:,iSession,iS))/n_RP_trials)*100;
            choiceHighE.(taskType)(iSession, iS) = (sum(choiceHighEPerTrial.(taskType)(:,iSession,iS))/n_trials)*100;
        end % session loop
        
        choiceHighR_total.(taskType)(iS) = (sum(choiceHighRPerTrial.(taskType)(:,1,iS) + choiceHighRPerTrial.(taskType)(:,2,iS))/n_trials)*100;
        choiceLowP_total.(taskType)(iS) = (sum(choiceLowPPerTrial.(taskType)(:,1,iS) + choiceLowPPerTrial.(taskType)(:,2,iS))/n_trials)*100;
        choiceHighE_total.(taskType)(iS) = (sum(choiceHighEPerTrial.(taskType)(:,1,iS) + choiceHighEPerTrial.(taskType)(:,2,iS))/n_trials*2)*100;
    end % task loop
end % subject loop

% mental effort
fig;
for iS = 1:NS
    bar(iS, choiceHighR_total.mental(iS));
    hold on;
end
line([0 NS+1],[100 100],'LineWidth',3,'Color','k');
ylim([0 110]);
xticks();
xticklabels({'','','',''});
xlabel('pilot subjects');
ylabel('Choice high R (%)');
legend_size(40);

% physical effort
fig;
for iS = 1:NS
    bar(iS, choiceHighR_total.physical(iS));
    hold on;
end
line([0 NS+1],[100 100],'LineWidth',3,'Color','k');
ylim([0 110]);
xticks();
xticklabels({'','','',''});
xlabel('pilot subjects');
ylabel('Choice high R (%)');
legend_size(40);

% mental effort
fig;
for iS = 1:NS
    bar(iS, choiceLowP_total.mental(iS));
    hold on;
end
line([0 NS+1],[100 100],'LineWidth',3,'Color','k');
ylim([0 110]);
xticks();
xticklabels({'','','',''});
xlabel('pilot subjects');
ylabel('Choice low P (%)');
legend_size(40);

% physical effort
fig;
for iS = 1:NS
    bar(iS, choiceLowP_total.physical(iS));
    hold on;
end
line([0 NS+1],[100 100],'LineWidth',3,'Color','k');
ylim([0 110]);
xticks();
xticklabels({'','','',''});
xlabel('pilot subjects');
ylabel('Choice ¨low P (%)');
legend_size(40);