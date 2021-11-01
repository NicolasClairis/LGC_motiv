%% average choice results out

%% list of subjects
subs = {'201','202','203','204','205','206','207','208'}; % behavioral pilots
NS = length(subs);
% fMRI participants

%% main parameters
figDisp = 0;
n_bins = 6;
[avg_defaultChoice, sd_defaultChoice] = cell(1,NS);

%% extract individual data
for iS = 1:NS
    sub_nm = subs{iS};
    subFolder = ['LGC_Motiv_results',filesep,'CID',sub_nm,filesep,'behavior'];
    cd(subFolder);
    [avg_defaultChoice{iS}, sd_defaultChoice{iS}] = checkChoiceDefaultProportion(sub_nm, figDisp, n_bins);
    
end % subject loop

%% display the average

%% figures
if figDisp == 1
    error('add SEM across participants for visual display');
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