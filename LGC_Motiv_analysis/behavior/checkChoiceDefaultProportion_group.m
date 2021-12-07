function[avg_defaultChoice, sem_defaultChoice,...
    avg_conf, sem_conf] = checkChoiceDefaultProportion_group(study_nm, n_bins)
% [avg_defaultChoice, sem_defaultChoice,...
%     avg_conf, sem_conf] = checkChoiceDefaultProportion_group(study_nm, n_bins)
%% checkChoiceDefaultProportion_group checks average proportion of choosing 
% the default option and the average level of confidence across participants.
% 
% INPUTS
% study_nm: which group of subjects do you want to check?
%
% n_bins: number of bins for your graphs?
%
% OUTPUTS
% avg_defaultChoice: sem_defaultChoice: structure with information about
% average and standard error of the mean (across participants) proportion 
% of choosing the default option in function of the other variables of
% interest of the experiment.
%
% avg_conf: sem_conf: structure with information about
% average and standard error of the mean (across participants) confidence
% over the choice to be made in function of the other variables of
% interest of the experiment.
%
% See also checkChoiceDefaultProportion.m (same but for individual subject)

%% list of subjects
if ~exist('study_nm','var') || isempty(study_nm)
    study_names = {'fMRI_pilots','study1','study2'};
    study_nm_idx = listdlg('ListString',study_names);
    study_nm = study_names{study_nm_idx};
end
[subs, NS] = LGCM_subject_selection(study_nm);

%% working directories
switch study_nm
    case 'behavioral_pilots'
        root = fullfile('C:','Users','clairis','Desktop','GitHub','LGC_motiv','LGC_Motiv_results');
    otherwise
        computer_root = LGCM_root_paths();
        switch study_nm
            case 'fMRI_pilots'
                root = fullfile(computer_root,'fMRI_pilots');
            case 'study1'
                root = fullfile(computer_root,'study1');
            case 'study2'
                root = fullfile(computer_root,'study2');
        end
end

%% main parameters
figDisp = 0;
if ~exist('n_bins','var') || isempty(n_bins)
    n_bins = 6; % default number of bins if not entered in the inputs
end

[avg_defaultChoice_data_perSub, avg_conf_data_perSub] = deal(cell(1,NS));
n_R_levels = 4;
n_E_levels = 4;

taskTypes = {'Ep','Em'};
nTaskTypes = length(taskTypes);

for iTask = 1:nTaskTypes
    task_nm = taskTypes{iTask};
    % ratio of choosing default per task
    [defaultChoice_perSub.(task_nm)] = deal(NaN(1,NS));
    [defaultChoice_perSub.([task_nm,'_R']),...
        defaultChoice_perSub.([task_nm,'_P'])] = deal(NaN(1,NS));
    [defaultChoice_perSub.([task_nm,'_f_time']),...
        conf_perSub.([task_nm,'_f_time'])] = deal(NaN(n_bins, NS));
    % per task per effort level
    for iE = 1:(n_E_levels - 1)
        [defaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)]),...
            conf_perSub.perElevel.([task_nm,'_',num2str(iE)])] = deal(NaN(1,NS));
    end
    % per task per money level
    for iAbsMoney = 1:(n_R_levels - 1)
        [defaultChoice_perSub.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),...
            defaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),...
            defaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),...
            defaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),...
            conf_perSub.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),...
            conf_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),...
            conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),...
            conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])] = deal(NaN(1,NS));
    end
    
    %% extract individual data
    for iS = 1:NS
        sub_nm = subs{iS};
        subFolder = [root,filesep,'CID',sub_nm,filesep,'behavior'];
        cd(subFolder);
        [avg_defaultChoice_data_perSub{iS}, ~, avg_conf_data_perSub{iS}] = checkChoiceDefaultProportion(sub_nm, figDisp, n_bins);
        % extract the data
        defaultChoice_perSub.(task_nm)(iS) = avg_defaultChoice_data_perSub{iS}.(task_nm);
        defaultChoice_perSub.([task_nm,'_R'])(iS) = avg_defaultChoice_data_perSub{iS}.([task_nm,'_R']);
        defaultChoice_perSub.([task_nm,'_P'])(iS) = avg_defaultChoice_data_perSub{iS}.([task_nm,'_P']);
        defaultChoice_perSub.([task_nm,'_f_time'])(:,iS) = avg_defaultChoice_data_perSub{iS}.([task_nm,'_f_time']);
        conf_perSub.([task_nm,'_f_time'])(:,iS) = avg_conf_data_perSub{iS}.([task_nm,'_f_time']);
        % per task per effort level
        for iE = 1:(n_E_levels - 1)
            defaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)])(iS) = avg_defaultChoice_data_perSub{iS}.perElevel.([task_nm,'_',num2str(iE)]);
            conf_perSub.perElevel.([task_nm,'_',num2str(iE)])(iS) = avg_conf_data_perSub{iS}.perElevel.([task_nm,'_',num2str(iE)]);
        end
        % per task per money level
        for iAbsMoney = 1:(n_R_levels - 1)
            % default choice proportion
            defaultChoice_perSub.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]);
            defaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]);
            defaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]);
            defaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(iS) = avg_defaultChoice_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]);
            % same for confidence
            conf_perSub.perAbsMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]);
            conf_perSub.perDeltaMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]);
            conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]);
            conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]);
        end
    end % subject loop
    
    %% extract the average across participants
    avg_defaultChoice.(task_nm) = mean(defaultChoice_perSub.(task_nm),2);
    avg_defaultChoice.([task_nm,'_R']) = mean(defaultChoice_perSub.([task_nm,'_R']),2);
    avg_defaultChoice.([task_nm,'_P']) = mean(defaultChoice_perSub.([task_nm,'_P']),2);
    % SEM
    sem_defaultChoice.(task_nm) = sem(defaultChoice_perSub.(task_nm),2);
    sem_defaultChoice.([task_nm,'_R']) = sem(defaultChoice_perSub.([task_nm,'_R']),2);
    sem_defaultChoice.([task_nm,'_P']) = sem(defaultChoice_perSub.([task_nm,'_P']),2);
    
    % fatigue
    avg_defaultChoice.([task_nm,'_f_time']) = mean(defaultChoice_perSub.([task_nm,'_f_time']),2);
    avg_conf.([task_nm,'_f_time']) = mean(conf_perSub.([task_nm,'_f_time']),2);
    sem_defaultChoice.([task_nm,'_f_time']) = sem(defaultChoice_perSub.([task_nm,'_f_time']),2);
    sem_conf.([task_nm,'_f_time']) = sem(conf_perSub.([task_nm,'_f_time']),2);
    
    % per task per effort level
    for iE = 1:(n_E_levels - 1)
        avg_defaultChoice.perElevel.([task_nm,'_',num2str(iE)]) = mean(defaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)]),2);
        avg_conf.perElevel.([task_nm,'_',num2str(iE)]) = mean(conf_perSub.perElevel.([task_nm,'_',num2str(iE)]),2,'omitnan');
        % SEM
        sem_defaultChoice.perElevel.([task_nm,'_',num2str(iE)]) = sem(defaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)]),2);
        sem_conf.perElevel.([task_nm,'_',num2str(iE)]) = sem(conf_perSub.perElevel.([task_nm,'_',num2str(iE)]),2);
    end
    
    % per task per money level
    for iAbsMoney = 1:(n_R_levels - 1)
        % mean default choice proportion
        avg_defaultChoice.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        avg_defaultChoice.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        avg_defaultChoice.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2);
        avg_defaultChoice.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = mean(defaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2);
        % mean confidence
        avg_conf.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = mean(conf_perSub.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        avg_conf.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = mean(conf_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        avg_conf.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = mean(conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2);
        avg_conf.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = mean(conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2);
        % sem choice proportion
        sem_defaultChoice.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        sem_defaultChoice.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        sem_defaultChoice.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2);
        sem_defaultChoice.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = sem(defaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2);
        % sem confidence
        sem_conf.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = sem(conf_perSub.perAbsMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        sem_conf.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = sem(conf_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        sem_conf.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = sem(conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2);
        sem_conf.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = sem(conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2);
    end
    
end % loop through physical and mental effort

%% figures
% figure parameters
lWidth_50percentTrait = 2;
pSize = 30;
bWidth = 0.4;
bDist = 0.2;
Em_col = [0 1 0];
Ep_col = [0 153/255 1];

% check choices = f(fatigue, R/P, level of R, level of E)
%% check choices = f(fatigue)
fig;
% mark the 50% trait
plot(1:n_bins, 50*ones(1,n_bins),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
% mental effort
bar_hdl.Em = jbfill(1:n_bins,...
    (avg_defaultChoice.Em_f_time + sem_defaultChoice.Em_f_time)'.*100,...
    (avg_defaultChoice.Em_f_time - sem_defaultChoice.Em_f_time)'.*100,...
    avg_defaultChoice.Em_f_time'.*100,...
    Em_col);
% physical effort
bar_hdl.Ep = jbfill(1:n_bins,...
    (avg_defaultChoice.Ep_f_time + sem_defaultChoice.Ep_f_time)'.*100,...
    (avg_defaultChoice.Ep_f_time - sem_defaultChoice.Ep_f_time)'.*100,...
    avg_defaultChoice.Ep_f_time'.*100,...
    Ep_col);
ylim([0 100]);
xlim([0 n_bins+1]);
xlabel('trial bins');
ylabel('Choice = default option (%)');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check choices = f(R/P) per task
fig;
% mark the 50% trait
plot(0:5, 50*ones(1,6),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
% R mental effort
bar_hdl.Em = bar(1, avg_defaultChoice.Em_R,'FaceColor',Em_col);
errorbar(1, avg_defaultChoice.Em_R, sem_defaultChoice.Em_R,'k');
% P mental effort
bar(2, avg_defaultChoice.Em_P,'FaceColor',Em_col);
errorbar(2, avg_defaultChoice.Em_P, sem_defaultChoice.Em_P,'k');
% R physical effort
bar_hdl.Ep = bar(3, avg_defaultChoice.Ep_R,'FaceColor',Ep_col);
errorbar(3, avg_defaultChoice.Ep_R, sem_defaultChoice.Ep_R,'k');
% P physical effort
bar(4, avg_defaultChoice.Ep_P,'FaceColor',Ep_col);
errorbar(4, avg_defaultChoice.Ep_P, sem_defaultChoice.Ep_P,'k');
xticks(1:4);
xticklabels({'R','P','R','P'});
ylim([0 100]);
xlim([0 5]);
ylabel('Choice = default option (%)');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check choices = f(|money level|)
fig;
% mark the 50% trait
plot(0:n_bins, 50*ones(1,n_bins+1),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
for iAbsMoney = 1:(n_R_levels - 1)
    bar_hdl.Em = bar(iAbsMoney-bDist, avg_defaultChoice.perAbsMoneylevel.(['Em_',num2str(iAbsMoney)]), 'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iAbsMoney+bDist, avg_defaultChoice.perAbsMoneylevel.(['Ep_',num2str(iAbsMoney)]), 'FaceColor',Ep_col,'BarWidth',bWidth);
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
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check choices = f(|delta money| level)
fig;
% mark the 50% trait
plot(0:n_bins, 50*ones(1,n_bins+1),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
for iAbsMoney = 1:(n_R_levels - 1)
    bar_hdl.Em = bar(iAbsMoney-bDist, avg_defaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]), 'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iAbsMoney+bDist, avg_defaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]), 'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iAbsMoney-bDist, avg_defaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]), sem_defaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),'k')
    errorbar(iAbsMoney+bDist, avg_defaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]), sem_defaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),'k')
end
ylim([0 100]);
ylabel('Choice = default option (%)');
xticks(1:3);
xticklabels({'1','2','3'});
xlim([0 n_R_levels]);
xlabel('|Δ money| level');
legend_size(pSize);
ylim([0 100]);
ylabel('Choice = default option (%)');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check choices = f(money levels (splitting R and P trials))
fig;
% mark the 50% trait
plot(-n_bins:n_bins, 50*ones(1,(n_bins*2)+1),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
    jMoney = abs(iMoney);
    if iMoney < 0
        taskCond = 'P';
    elseif iMoney > 0
        taskCond = 'R';
    end
    bar_hdl.Em = bar(iMoney-bDist, avg_defaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), 'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iMoney+bDist, avg_defaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), 'FaceColor',Ep_col,'BarWidth',bWidth);
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
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check choices = f(E level)
fig;
% mark the 50% trait
plot(0:n_bins, 50*ones(1,n_bins+1),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
for iE = 1:(n_E_levels - 1)
    bar_hdl.Em = bar(iE-bDist, avg_defaultChoice.perElevel.(['Em_',num2str(iE)]), 'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iE+bDist, avg_defaultChoice.perElevel.(['Ep_',num2str(iE)]), 'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iE-bDist, avg_defaultChoice.perElevel.(['Em_',num2str(iE)]), sem_defaultChoice.perElevel.(['Em_',num2str(iE)]),'k')
    errorbar(iE+bDist, avg_defaultChoice.perElevel.(['Ep_',num2str(iE)]), sem_defaultChoice.perElevel.(['Ep_',num2str(iE)]),'k')
end
ylim([0 100]);
ylabel('Choice = default option (%)');
xticks(1:3);
xticklabels({'1','2','3'});
xlim([0 n_E_levels]);
xlabel('Effort level');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check confidence = f(money levels (splitting R and P trials))
fig;
% mark the 50% trait
plot(-n_bins:n_bins, 50*ones(1,(n_bins*2)+1),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
    jMoney = abs(iMoney);
    if iMoney < 0
        taskCond = 'P';
    elseif iMoney > 0
        taskCond = 'R';
    end
    bar_hdl.Em = bar(iMoney-bDist, avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), 'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iMoney+bDist, avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), 'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iMoney-bDist, avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), sem_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),'k')
    errorbar(iMoney+bDist, avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), sem_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),'k')
end
xticks([-3:-(1), 1:3]);
xticklabels({'-3','-2','-1','1','2','3'});
xlim([-n_R_levels n_R_levels]);
ylim([0 100]);
ylabel('Level of confidence');
xlabel('Money level');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check confidence = f(E levels)
fig;
% mark the 50% trait
plot(0:n_bins, 50*ones(1,n_bins+1),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
for iE = 1:(n_E_levels - 1)
    bar_hdl.Em = bar(iE-bDist, avg_conf.perElevel.(['Em_',num2str(iE)]), 'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iE+bDist, avg_conf.perElevel.(['Ep_',num2str(iE)]), 'FaceColor',Ep_col,'BarWidth',bWidth);
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
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check confidence = f(time)
fig;
% mark the 50% trait
plot(1:n_bins, 50*ones(1,n_bins),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
% mental effort
bar_hdl.Em = jbfill(1:n_bins,...
    (avg_conf.Em_f_time + sem_conf.Em_f_time)'.*100,...
    (avg_conf.Em_f_time - sem_conf.Em_f_time)'.*100,...
    avg_conf.Em_f_time'.*100,...
    Em_col);
% physical effort
bar_hdl.Ep = jbfill(1:n_bins,...
    (avg_conf.Ep_f_time + sem_conf.Ep_f_time)'.*100,...
    (avg_conf.Ep_f_time - sem_conf.Ep_f_time)'.*100,...
    avg_conf.Ep_f_time'.*100,...
    Ep_col);
ylim([0 100]);
xlim([0 n_bins+1]);
xlabel('trial bins');
ylabel('Confidence (%)');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

end % function