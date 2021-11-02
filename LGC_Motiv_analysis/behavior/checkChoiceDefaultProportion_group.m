%% average choice results out

%% working directories
root = fullfile('C:','Users','clairis','Desktop','GitHub','LGC_motiv','LGC_Motiv_results');

%% list of subjects
subs = {'201','202','203','204','205','206','207','208','209'}; % behavioral pilots
NS = length(subs);
% fMRI participants

%% main parameters
figDisp = 0;
n_bins = 6;
[avg_defaultChoice_data_perSub, avg_conf_data_perSub] = deal(cell(1,NS));
n_R_levels = 4;
n_E_levels = 4;

% ratio of choosing default per task
[defaultChoice_perSub.Em, defaultChoice_perSub.Ep] = deal(NaN(1,NS));
[defaultChoice_perSub.Em_f_time, defaultChoice_perSub.Ep_f_time,...
    conf_perSub.Em_f_time, conf_perSub.Ep_f_time] = deal(NaN(n_bins, NS));
% per task per effort level
for iE = 1:(n_E_levels - 1)
    [defaultChoice_perSub.perElevel.(['Em_',num2str(iE)]),...
        defaultChoice_perSub.perElevel.(['Ep_',num2str(iE)]),...
        conf_perSub.perElevel.(['Em_',num2str(iE)]),...
        conf_perSub.perElevel.(['Ep_',num2str(iE)])] = deal(NaN(1,NS));
end
% per task per money level
for iAbsMoney = 1:(n_R_levels - 1)
    [defaultChoice_perSub.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]),...
        defaultChoice_perSub.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
        defaultChoice_perSub.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]),...
        defaultChoice_perSub.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]),...
        defaultChoice_perSub.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]),...
        defaultChoice_perSub.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]),...
        conf_perSub.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]),...
        conf_perSub.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]),...
        conf_perSub.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]),...
        conf_perSub.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)])] = deal(NaN(1,NS));
end

%% extract individual data
for iS = 1:NS
    sub_nm = subs{iS};
    subFolder = [root,filesep,'CID',sub_nm,filesep,'behavior'];
    cd(subFolder);
    [avg_defaultChoice_data_perSub{iS}, ~, avg_conf_data_perSub{iS}] = checkChoiceDefaultProportion(sub_nm, figDisp, n_bins);
    % extract the data
    defaultChoice_perSub.Em(iS) = avg_defaultChoice_data_perSub{iS}.Em;
    defaultChoice_perSub.Ep(iS) = avg_defaultChoice_data_perSub{iS}.Ep;
    defaultChoice_perSub.Em_f_time(:,iS) = avg_defaultChoice_data_perSub{iS}.Em_f_time;
    defaultChoice_perSub.Ep_f_time(:,iS) = avg_defaultChoice_data_perSub{iS}.Ep_f_time;
    conf_perSub.Em_f_time(:,iS) = avg_conf_data_perSub{iS}.Em_f_time;
    conf_perSub.Ep_f_time(:,iS) = avg_conf_data_perSub{iS}.Ep_f_time;
    % per task per effort level
    for iE = 1:(n_E_levels - 1)
        defaultChoice_perSub.perElevel.(['Em_',num2str(iE)])(iS) = avg_defaultChoice_data_perSub{iS}.perElevel.(['Em_',num2str(iE)]);
        defaultChoice_perSub.perElevel.(['Ep_',num2str(iE)])(iS) = avg_defaultChoice_data_perSub{iS}.perElevel.(['Ep_',num2str(iE)]);
        conf_perSub.perElevel.(['Em_',num2str(iE)])(iS) = avg_conf_data_perSub{iS}.perElevel.(['Em_',num2str(iE)]);
        conf_perSub.perElevel.(['Ep_',num2str(iE)])(iS) = avg_conf_data_perSub{iS}.perElevel.(['Ep_',num2str(iE)]);
    end
    % per task per money level
    for iAbsMoney = 1:(n_R_levels - 1)
        defaultChoice_perSub.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]);
        defaultChoice_perSub.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]);
        defaultChoice_perSub.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]);
        defaultChoice_perSub.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]);
        defaultChoice_perSub.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]);
        defaultChoice_perSub.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]);
        conf_perSub.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]);
        conf_perSub.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]);
        conf_perSub.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]);
        conf_perSub.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]);
    end
end % subject loop

%% extract the average across participants
avg_defaultChoice.Em = mean(defaultChoice_perSub.Em,2);
avg_defaultChoice.Ep = mean(defaultChoice_perSub.Ep,2);
% SEM
sem_defaultChoice.Em = sem(defaultChoice_perSub.Em,2);
sem_defaultChoice.Ep = sem(defaultChoice_perSub.Ep,2);

% fatigue
avg_defaultChoice.Em_f_time = mean(defaultChoice_perSub.Em_f_time,2);
avg_defaultChoice.Ep_f_time = mean(defaultChoice_perSub.Ep_f_time,2);
avg_conf.Em_f_time = mean(conf_perSub.Em_f_time,2);
avg_conf.Ep_f_time = mean(conf_perSub.Ep_f_time,2);
sem_defaultChoice.Em_f_time = sem(defaultChoice_perSub.Em_f_time,2);
sem_defaultChoice.Ep_f_time = sem(defaultChoice_perSub.Ep_f_time,2);
sem_conf.Em_f_time = sem(conf_perSub.Em_f_time,2);
sem_conf.Ep_f_time = sem(conf_perSub.Ep_f_time,2);

% per task per effort level
for iE = 1:(n_E_levels - 1)
    avg_defaultChoice.perElevel.(['Em_',num2str(iE)]) = mean(defaultChoice_perSub.perElevel.(['Em_',num2str(iE)]),2);
    avg_defaultChoice.perElevel.(['Ep_',num2str(iE)]) = mean(defaultChoice_perSub.perElevel.(['Ep_',num2str(iE)]),2);
    avg_conf.perElevel.(['Em_',num2str(iE)]) = mean(conf_perSub.perElevel.(['Em_',num2str(iE)]),2,'omitnan');
    avg_conf.perElevel.(['Ep_',num2str(iE)]) = mean(conf_perSub.perElevel.(['Ep_',num2str(iE)]),2,'omitnan');
    % SEM
    sem_defaultChoice.perElevel.(['Em_',num2str(iE)]) = sem(defaultChoice_perSub.perElevel.(['Em_',num2str(iE)]),2);
    sem_defaultChoice.perElevel.(['Ep_',num2str(iE)]) = sem(defaultChoice_perSub.perElevel.(['Ep_',num2str(iE)]),2);
    sem_conf.perElevel.(['Em_',num2str(iE)]) = sem(conf_perSub.perElevel.(['Em_',num2str(iE)]),2);
    sem_conf.perElevel.(['Ep_',num2str(iE)]) = sem(conf_perSub.perElevel.(['Ep_',num2str(iE)]),2);
end

% per task per money level
for iAbsMoney = 1:(n_R_levels - 1)
    avg_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]),2);
    avg_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]),2);
    avg_defaultChoice.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]),2);
    avg_defaultChoice.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]),2);
    avg_defaultChoice.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]),2);
    avg_defaultChoice.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]),2);
    avg_conf.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = mean(conf_perSub.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]),2);
    avg_conf.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = mean(conf_perSub.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]),2);
    avg_conf.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = mean(conf_perSub.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]),2);
    avg_conf.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = mean(conf_perSub.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]),2);
    % sem
    sem_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]),2);
    sem_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]),2);
    sem_defaultChoice.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]),2);
    sem_defaultChoice.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]),2);
    sem_defaultChoice.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]),2);
    sem_defaultChoice.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]),2);
    sem_conf.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]) = sem(conf_perSub.perSignedMoneylevel.(['Em_R_',num2str(iAbsMoney)]),2);
    sem_conf.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]) = sem(conf_perSub.perSignedMoneylevel.(['Ep_R_',num2str(iAbsMoney)]),2);
    sem_conf.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]) = sem(conf_perSub.perSignedMoneylevel.(['Em_P_',num2str(iAbsMoney)]),2);
    sem_conf.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]) = sem(conf_perSub.perSignedMoneylevel.(['Ep_P_',num2str(iAbsMoney)]),2);
end

%% figures
warning('add SEM across participants for visual display');
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
% mental effort
jbfill(1:n_bins,...
    (avg_defaultChoice.Em_f_time + sem_defaultChoice.Em_f_time)'.*100,...
    (avg_defaultChoice.Em_f_time - sem_defaultChoice.Em_f_time)'.*100,...
    avg_defaultChoice.Em_f_time'.*100,...
    'g');
% physical effort
jbfill(1:n_bins,...
    (avg_defaultChoice.Ep_f_time + sem_defaultChoice.Ep_f_time)'.*100,...
    (avg_defaultChoice.Ep_f_time - sem_defaultChoice.Ep_f_time)'.*100,...
    avg_defaultChoice.Ep_f_time'.*100,...
    'b');
ylim([0 100]);
xlim([0 n_bins+1]);
xlabel('trial bins');
ylabel('Choice = default option (%)');
legend_size(pSize);

% %% check choices = f(R/P)
% fig;
% % mark the 50% trait
% plot(0:3, 50*ones(1,4),...
%     'LineWidth',2,'Color','k','LineStyle',':');
% hold on;
% bar(1, mean(percentage_R_choiceDefault.*100),'FaceColor',[0 153/255 1]);
% errorbar(1, mean(percentage_R_choiceDefault.*100), sem(percentage_R_choiceDefault.*100),'k');
% bar(2, mean(percentage_P_choiceDefault.*100),'FaceColor',[1 153/255 0]);
% errorbar(2, mean(percentage_P_choiceDefault.*100), sem(percentage_P_choiceDefault.*100),'k');
% xticks(1:2);
% xticklabels({'R','P'});
% ylim([0 100]);
% xlim([0 3]);
% ylabel('Choice = default option (%)');
% legend_size(pSize);


%% check choices = f(money level)
fig;
% mark the 50% trait
plot(0:n_bins, 50*ones(1,n_bins+1),...
    'LineWidth',2,'Color','k','LineStyle',':');
hold on;
for iAbsMoney = 1:(n_R_levels - 1)
    bar(iAbsMoney-bDist, avg_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]), 'FaceColor','g','BarWidth',bWidth);
    bar(iAbsMoney+bDist, avg_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]), 'FaceColor','b','BarWidth',bWidth);
    errorbar(iAbsMoney-bDist, avg_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]), sem_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]),'k')
    errorbar(iAbsMoney+bDist, avg_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]), sem_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]),'k')
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
    errorbar(iE-bDist, avg_defaultChoice.perElevel.(['Em_',num2str(iE)]), sem_defaultChoice.perElevel.(['Em_',num2str(iE)]),'k')
    errorbar(iE+bDist, avg_defaultChoice.perElevel.(['Ep_',num2str(iE)]), sem_defaultChoice.perElevel.(['Ep_',num2str(iE)]),'k')
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
    errorbar(iMoney-bDist, avg_defaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), sem_defaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),'k')
    errorbar(iMoney+bDist, avg_defaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), sem_defaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),'k')
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
    errorbar(iMoney-bDist, avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), sem_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),'k')
    errorbar(iMoney+bDist, avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), sem_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),'k')
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
    errorbar(iE-bDist, avg_conf.perElevel.(['Em_',num2str(iE)]), sem_conf.perElevel.(['Em_',num2str(iE)]),'k')
    errorbar(iE+bDist, avg_conf.perElevel.(['Ep_',num2str(iE)]), sem_conf.perElevel.(['Ep_',num2str(iE)]),'k')
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
% mental effort
jbfill(1:n_bins,...
    (avg_conf.Em_f_time + sem_conf.Em_f_time)'.*100,...
    (avg_conf.Em_f_time - sem_conf.Em_f_time)'.*100,...
    avg_conf.Em_f_time'.*100,...
    'g');
% physical effort
jbfill(1:n_bins,...
    (avg_conf.Ep_f_time + sem_conf.Ep_f_time)'.*100,...
    (avg_conf.Ep_f_time - sem_conf.Ep_f_time)'.*100,...
    avg_conf.Ep_f_time'.*100,...
    'b');
ylim([0 100]);
xlim([0 n_bins+1]);
xlabel('trial bins');
ylabel('Confidence (%)');
legend_size(pSize);
