function[avg_nonDefaultChoice, sem_nonDefaultChoice,...
    avg_conf, sem_conf,...
    avg_RT, sem_RT,...
    avg_Eperformance, sem_Eperformance] = checkChoiceNonDefaultProportion_group(study_nm, n_bins)
% [avg_nonDefaultChoice, sem_nonDefaultChoice,...
%     avg_conf, sem_conf,...
%     avg_RT, sem_RT,...
%     avg_Eperformance, sem_Eperformance] = checkChoiceNonDefaultProportion_group(study_nm, n_bins)
%% checkChoiceNonDefaultProportion_group checks average proportion of choosing 
% the non-default option and the average level of confidence across participants.
% 
% INPUTS
% study_nm: which group of subjects do you want to check?
%
% n_bins: number of bins for your graphs?
%
% OUTPUTS
% avg_nonDefaultChoice: sem_nonDefaultChoice: structure with information about
% average and standard error of the mean (across participants) proportion 
% of choosing the non-default option in function of the other variables of
% interest of the experiment.
%
% avg_conf, sem_conf: structure with information about
% average and standard error of the mean (across participants) confidence
% over the choice to be made in function of the other variables of
% interest of the experiment.
%
% avg_RT, sem_RT: structure with information about
% average and standard error of the mean (across participants) reaction
% times over the choice to be made in function of the other variables of
% interest of the experiment.
%
% avg_Eperformance, sem_Eperformance: structure with information about
% average and standard error of the mean (across participants) for effort
% performance in function of the other variables of interest of the experiment.
%
% See also checkChoiceNonDefaultProportion.m (same but for individual subject)

%% list of subjects
[study_nm,~,~,subject_id,NS] = sub_id;

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

[avg_nonDefaultChoice_data_perSub,...
    avg_conf_data_perSub,...
    avg_RT_data_perSub,...
    avg_Eperformance_data_perSub] = deal(cell(1,NS));
n_R_levels = 4;
n_E_levels = 4;

nTrialsPerRun = 54;

taskTypes = {'Ep','Em'};
nTaskTypes = length(taskTypes);

for iTask = 1:nTaskTypes
    task_nm = taskTypes{iTask};
    % ratio of choosing default per task
    [nonDefaultChoice_perSub.(task_nm)] = deal(NaN(1,NS));
    [nonDefaultChoice_perSub.([task_nm,'_R']),...
        nonDefaultChoice_perSub.([task_nm,'_P'])] = deal(NaN(1,NS));
    [nonDefaultChoice_perSub.([task_nm,'_f_time']),...
        conf_perSub.([task_nm,'_f_time']),...
        RT_perSub.([task_nm,'_f_time']),...
        Eperformance_perSub.([task_nm,'_f_time'])] = deal(NaN(n_bins, NS));
    % per task per effort level
    for iE = 1:(n_E_levels - 1)
        [nonDefaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)]),...
            conf_perSub.perElevel.([task_nm,'_',num2str(iE)]),...
            RT_perSub.perElevel.([task_nm,'_',num2str(iE)]),...
            Eperformance_perSub.perElevel.([task_nm,'_',num2str(iE)])] = deal(NaN(1,NS));
    end
    % per task per money level
    for iAbsMoney = 1:(n_R_levels - 1)
        [nonDefaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),...
            nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),...
            nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),...
            conf_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),...
            conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),...
            conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),...
            RT_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),...
            RT_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),...
            RT_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),...
            Eperformance_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),...
            Eperformance_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),...
            Eperformance_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])] = deal(NaN(1,NS));
    end
    
    %% extract individual data
    for iS = 1:NS
        sub_nm = subject_id{iS};
        subFolder = [root,filesep,'CID',sub_nm,filesep,'behavior'];
        cd(subFolder);
        [avg_nonDefaultChoice_data_perSub{iS},...
            ~,...
            avg_conf_data_perSub{iS},...
            avg_RT_data_perSub{iS},...
            avg_Eperformance_data_perSub{iS}] = checkChoiceNonDefaultProportion(sub_nm, figDisp, n_bins);
        % extract the data
        nonDefaultChoice_perSub.(task_nm)(iS) = avg_nonDefaultChoice_data_perSub{iS}.(task_nm);
        nonDefaultChoice_perSub.([task_nm,'_R'])(iS) = avg_nonDefaultChoice_data_perSub{iS}.([task_nm,'_R']);
        nonDefaultChoice_perSub.([task_nm,'_P'])(iS) = avg_nonDefaultChoice_data_perSub{iS}.([task_nm,'_P']);
        nonDefaultChoice_perSub.([task_nm,'_f_time'])(:,iS) = avg_nonDefaultChoice_data_perSub{iS}.([task_nm,'_f_time']);
        conf_perSub.([task_nm,'_f_time'])(:,iS) = avg_conf_data_perSub{iS}.([task_nm,'_f_time']);
        RT_perSub.([task_nm,'_f_time'])(:,iS) = avg_RT_data_perSub{iS}.([task_nm,'_f_time']);
        Eperformance_perSub.([task_nm,'_f_time'])(:,iS) = avg_Eperformance_data_perSub{iS}.([task_nm,'_f_time']);
        % per task per effort level
        for iE = 1:(n_E_levels - 1)
            nonDefaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)])(iS) = avg_nonDefaultChoice_data_perSub{iS}.perElevel.([task_nm,'_',num2str(iE)]);
            conf_perSub.perElevel.([task_nm,'_',num2str(iE)])(iS) = avg_conf_data_perSub{iS}.perElevel.([task_nm,'_',num2str(iE)]);
            RT_perSub.perElevel.([task_nm,'_',num2str(iE)])(iS) = avg_RT_data_perSub{iS}.perElevel.([task_nm,'_',num2str(iE)]);
            Eperformance_perSub.perElevel.([task_nm,'_',num2str(iE)])(iS) = avg_Eperformance_data_perSub{iS}.perElevel.([task_nm,'_',num2str(iE)]);
        end
        % per task per money level
        for iAbsMoney = 1:(n_R_levels - 1)
            % default choice proportion
            nonDefaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(iS) = avg_nonDefaultChoice_data_perSub{iS}.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]);
            nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(iS) = avg_nonDefaultChoice_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]);
            nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(iS) = avg_nonDefaultChoice_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]);
            % same for confidence
            conf_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]);
            conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]);
            conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(iS) = avg_conf_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]);
            % same for RT
            RT_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(iS) = avg_RT_data_perSub{iS}.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]);
            RT_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(iS) = avg_RT_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]);
            RT_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(iS) = avg_RT_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]);
            % same for effort performance
            Eperformance_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(iS) = avg_Eperformance_data_perSub{iS}.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]);
            Eperformance_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(iS) = avg_Eperformance_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]);
            Eperformance_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(iS) = avg_Eperformance_data_perSub{iS}.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]);
        end
    end % subject loop
    
    %% extract the average across participants
    avg_nonDefaultChoice.(task_nm) = mean(nonDefaultChoice_perSub.(task_nm),2,'omitnan');
    avg_nonDefaultChoice.([task_nm,'_R']) = mean(nonDefaultChoice_perSub.([task_nm,'_R']),2,'omitnan');
    avg_nonDefaultChoice.([task_nm,'_P']) = mean(nonDefaultChoice_perSub.([task_nm,'_P']),2,'omitnan');
    % SEM
    sem_nonDefaultChoice.(task_nm) = sem(nonDefaultChoice_perSub.(task_nm),2);
    sem_nonDefaultChoice.([task_nm,'_R']) = sem(nonDefaultChoice_perSub.([task_nm,'_R']),2);
    sem_nonDefaultChoice.([task_nm,'_P']) = sem(nonDefaultChoice_perSub.([task_nm,'_P']),2);
    
    % fatigue
    % mean
    avg_nonDefaultChoice.([task_nm,'_f_time']) = mean(nonDefaultChoice_perSub.([task_nm,'_f_time']),2,'omitnan');
    avg_conf.([task_nm,'_f_time']) = mean(conf_perSub.([task_nm,'_f_time']),2,'omitnan');
    avg_RT.([task_nm,'_f_time']) = mean(RT_perSub.([task_nm,'_f_time']),2,'omitnan');
    avg_Eperformance.([task_nm,'_f_time']) = mean(Eperformance_perSub.([task_nm,'_f_time']),2,'omitnan');
    % sem
    sem_nonDefaultChoice.([task_nm,'_f_time']) = sem(nonDefaultChoice_perSub.([task_nm,'_f_time']),2);
    sem_conf.([task_nm,'_f_time']) = sem(conf_perSub.([task_nm,'_f_time']),2);
    sem_RT.([task_nm,'_f_time']) = sem(RT_perSub.([task_nm,'_f_time']),2);
    sem_Eperformance.([task_nm,'_f_time']) = sem(Eperformance_perSub.([task_nm,'_f_time']),2);
    
    % per task per effort level
    for iE = 1:(n_E_levels - 1)
        % mean
        avg_nonDefaultChoice.perElevel.([task_nm,'_',num2str(iE)]) = mean(nonDefaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)]),2,'omitnan');
        avg_conf.perElevel.([task_nm,'_',num2str(iE)]) = mean(conf_perSub.perElevel.([task_nm,'_',num2str(iE)]),2,'omitnan');
        avg_RT.perElevel.([task_nm,'_',num2str(iE)]) = mean(RT_perSub.perElevel.([task_nm,'_',num2str(iE)]),2,'omitnan');
        avg_Eperformance.perElevel.([task_nm,'_',num2str(iE)]) = mean(Eperformance_perSub.perElevel.([task_nm,'_',num2str(iE)]),2,'omitnan');
        % SEM
        sem_nonDefaultChoice.perElevel.([task_nm,'_',num2str(iE)]) = sem(nonDefaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)]),2);
        sem_conf.perElevel.([task_nm,'_',num2str(iE)]) = sem(conf_perSub.perElevel.([task_nm,'_',num2str(iE)]),2);
        sem_RT.perElevel.([task_nm,'_',num2str(iE)]) = sem(RT_perSub.perElevel.([task_nm,'_',num2str(iE)]),2);
        sem_Eperformance.perElevel.([task_nm,'_',num2str(iE)]) = sem(Eperformance_perSub.perElevel.([task_nm,'_',num2str(iE)]),2);
    end
    
    % per task per money level
    for iAbsMoney = 1:(n_R_levels - 1)
        % mean default choice proportion
        avg_nonDefaultChoice.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = mean(nonDefaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2,'omitnan');
        avg_nonDefaultChoice.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = mean(nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2,'omitnan');
        avg_nonDefaultChoice.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = mean(nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2,'omitnan');
        % mean confidence
        avg_conf.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = mean(conf_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2,'omitnan');
        avg_conf.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = mean(conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2,'omitnan');
        avg_conf.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = mean(conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2,'omitnan');
        % mean RT
        avg_RT.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = mean(RT_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2,'omitnan');
        avg_RT.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = mean(RT_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2,'omitnan');
        avg_RT.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = mean(RT_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2,'omitnan');
        % mean effort performance
        avg_Eperformance.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = mean(Eperformance_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2,'omitnan');
        avg_Eperformance.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = mean(Eperformance_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2,'omitnan');
        avg_Eperformance.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = mean(Eperformance_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2,'omitnan');
        
        % sem choice proportion
        sem_nonDefaultChoice.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = sem(nonDefaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        sem_nonDefaultChoice.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = sem(nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2);
        sem_nonDefaultChoice.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = sem(nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2);
        % sem confidence
        sem_conf.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = sem(conf_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        sem_conf.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = sem(conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2);
        sem_conf.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = sem(conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2);
        % sem RT
        sem_RT.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = sem(RT_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        sem_RT.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = sem(RT_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2);
        sem_RT.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = sem(RT_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2);
        % sem effort performance
        sem_Eperformance.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]) = sem(Eperformance_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]),2);
        sem_Eperformance.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]) = sem(Eperformance_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]),2);
        sem_Eperformance.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]) = sem(Eperformance_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]),2);
    end
    
end % loop through physical and mental effort

%% figures
% figure parameters
lWidth_50percentTrait = 2;
lWidth = 3;
pSize = 40;
bWidth = 0.4;
bDist = 0.2;
Em_col = [0 1 0];
Ep_col = [0 153/255 1];

% check choices = f(fatigue, R/P, level of R, level of E)
%% check choices = f(fatigue) displaying all subjects
fig;
% mark the 50% trait
plot(1:n_bins, 50*ones(1,n_bins),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
for iS = 1:NS
    plot(1:n_bins, nonDefaultChoice_perSub.Em_f_time(:,iS).*100,...
        'Color',[0 1-0.01*iS 0], 'LineWidth',1)
    plot(1:n_bins, nonDefaultChoice_perSub.Ep_f_time(:,iS).*100,...
        'Color',[0 153/255 1-0.01*iS], 'LineWidth',1)
end
% ylim([0 100]);
ylim([50 85]);
xlim([1 n_bins]);
xticks(1:n_bins);
xLabelNames = cell(1,n_bins);
nTrialsPerBin = nTrialsPerRun/n_bins;
for iBin = 1:n_bins
    xLabelNames{iBin} = [num2str( 1 + nTrialsPerBin*(iBin-1)),'-',num2str(nTrialsPerBin*iBin)];
end
xticklabels(xLabelNames);
xlabel('Trial');
ylabel('Choice high effort (%)');
ylabel('Choice (%)');
legend_size(pSize);
%% check choices = f(fatigue)
fig;
% mark the 50% trait
plot(1:n_bins, 50*ones(1,n_bins),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
% mental effort
bar_hdl.Em = jbfill(1:n_bins,...
    (avg_nonDefaultChoice.Em_f_time + sem_nonDefaultChoice.Em_f_time)'.*100,...
    (avg_nonDefaultChoice.Em_f_time - sem_nonDefaultChoice.Em_f_time)'.*100,...
    avg_nonDefaultChoice.Em_f_time'.*100,...
    Em_col);
% physical effort
bar_hdl.Ep = jbfill(1:n_bins,...
    (avg_nonDefaultChoice.Ep_f_time + sem_nonDefaultChoice.Ep_f_time)'.*100,...
    (avg_nonDefaultChoice.Ep_f_time - sem_nonDefaultChoice.Ep_f_time)'.*100,...
    avg_nonDefaultChoice.Ep_f_time'.*100,...
    Ep_col);
% ylim([0 100]);
ylim([40 100]);
xlim([1 n_bins]);
xticks(1:n_bins);
xLabelNames = cell(1,n_bins);
nTrialsPerBin = nTrialsPerRun/n_bins;
for iBin = 1:n_bins
    xLabelNames{iBin} = [num2str( 1 + nTrialsPerBin*(iBin-1)),'-',num2str(nTrialsPerBin*iBin)];
end
xticklabels(xLabelNames);
xlabel('trial number');
ylabel('Choice high effort option (%)');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend('Location','NorthWest');
legend_size(pSize);

%% check choices = f(R/P) per task
fig;
% mark the 50% trait
plot(0:5, 50*ones(1,6),...
    'LineWidth',lWidth_50percentTrait,...
    'Color','k','LineStyle',':');
hold on;
% R mental effort
bar_hdl.Em = bar(1, avg_nonDefaultChoice.Em_R,'FaceColor',Em_col);
errorbar(1, avg_nonDefaultChoice.Em_R, sem_nonDefaultChoice.Em_R,...
    'k','LineWidth',lWidth);
% P mental effort
bar(2, avg_nonDefaultChoice.Em_P,'FaceColor',Em_col);
errorbar(2, avg_nonDefaultChoice.Em_P, sem_nonDefaultChoice.Em_P,...
    'k','LineWidth',lWidth);
% R physical effort
bar_hdl.Ep = bar(3, avg_nonDefaultChoice.Ep_R,'FaceColor',Ep_col);
errorbar(3, avg_nonDefaultChoice.Ep_R, sem_nonDefaultChoice.Ep_R,...
    'k','LineWidth',lWidth);
% P physical effort
bar(4, avg_nonDefaultChoice.Ep_P,'FaceColor',Ep_col);
errorbar(4, avg_nonDefaultChoice.Ep_P, sem_nonDefaultChoice.Ep_P,...
    'k','LineWidth',lWidth);
xticks(1:4);
xticklabels({'R','P','R','P'});
ylim([0 100]);
xlim([0 5]);
ylabel('Choice high effort option (%)');
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
    bar_hdl.Em = bar(iAbsMoney-bDist,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),...
        'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iAbsMoney+bDist,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
        'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iAbsMoney-bDist,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),...
        sem_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),...
        'k','LineWidth',lWidth);
    errorbar(iAbsMoney+bDist,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
        sem_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
        'k','LineWidth',lWidth);
end
ylim([0 100]);
ylabel('Choice high effort option (%)');
ylabel('High effort choice (%)');
xticks(1:3);
xticklabels({'1','2','3'});
xlim([0 n_R_levels]);
xlabel('|Δ money| level');
legend_size(pSize);
ylim([0 100]);
ylabel('Choice high effort option (%)');
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
    bar_hdl.Em = bar(iMoney-bDist,...
        avg_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iMoney+bDist,...
        avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iMoney-bDist,...
        avg_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        sem_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        'k','LineWidth',lWidth);
    errorbar(iMoney+bDist,...
        avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        sem_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        'k','LineWidth',lWidth);
end
ylim([0 100]);
ylabel('Choice high effort option (%)');
xticks([-3:-(1), 1:3]);
xticklabels({'-3','-2','-1','1','2','3'});
xlim([-n_R_levels n_R_levels]);
xlabel('Money level');
legend_size(pSize);
ylim([0 100]);
ylabel('Choice high effort option (%)');
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
    bar_hdl.Em = bar(iE-bDist, avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]), 'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iE+bDist, avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]), 'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iE-bDist,...
        avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]),...
        sem_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]),...
        'k','LineWidth',lWidth);
    errorbar(iE+bDist,...
        avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]),...
        sem_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]),...
        'k','LineWidth',lWidth);
end
ylim([0 100]);
ylabel('Choice high effort option (%)');
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
    errorbar(iMoney-bDist,...
        avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        sem_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        'k','LineWidth',lWidth);
    errorbar(iMoney+bDist,...
        avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        sem_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        'k','LineWidth',lWidth);
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
    bar_hdl.Em = bar(iE-bDist,...
        avg_conf.perElevel.(['Em_',num2str(iE)]),...
        'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iE+bDist,...
        avg_conf.perElevel.(['Ep_',num2str(iE)]),...
        'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iE-bDist,...
        avg_conf.perElevel.(['Em_',num2str(iE)]),...
        sem_conf.perElevel.(['Em_',num2str(iE)]),...
        'k','LineWidth',lWidth);
    errorbar(iE+bDist,...
        avg_conf.perElevel.(['Ep_',num2str(iE)]),...
        sem_conf.perElevel.(['Ep_',num2str(iE)]),...
        'k','LineWidth',lWidth);
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
xticks(1:n_bins);
xticklabels(xLabelNames);
xlabel('Trial number');
ylabel('Confidence (%)');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check RT = f(money levels (splitting R and P trials))
fig;
hold on;
for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
    jMoney = abs(iMoney);
    if iMoney < 0
        taskCond = 'P';
    elseif iMoney > 0
        taskCond = 'R';
    end
    bar_hdl.Em = bar(iMoney-bDist,...
        avg_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iMoney+bDist,...
        avg_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iMoney-bDist,...
        avg_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        sem_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        'k','LineWidth',lWidth)
    errorbar(iMoney+bDist,...
        avg_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        sem_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        'k','LineWidth',lWidth)
end
xticks([-3:-(1), 1:3]);
xticklabels({'-3','-2','-1','1','2','3'});
xlim([-n_R_levels n_R_levels]);
ylim([1.7 2.4]);
ylabel('RT (s)');
xlabel('Money level');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check RT = f(E levels)
fig;
hold on;
for iE = 1:(n_E_levels - 1)
    bar_hdl.Em = bar(iE-bDist, avg_RT.perElevel.(['Em_',num2str(iE)]),...
        'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iE+bDist, avg_RT.perElevel.(['Ep_',num2str(iE)]),...
        'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iE-bDist,...
        avg_RT.perElevel.(['Em_',num2str(iE)]), sem_RT.perElevel.(['Em_',num2str(iE)]),...
        'k','LineWidth',lWidth);
    errorbar(iE+bDist,...
        avg_RT.perElevel.(['Ep_',num2str(iE)]), sem_RT.perElevel.(['Ep_',num2str(iE)]),...
        'k','LineWidth',lWidth);
end
ylim([0 5]);
ylabel('RT (s)');
xticks(1:3);
xticklabels({'1','2','3'});
xlim([0 n_E_levels]);
ylim([1.8 2.3]);
xlabel('Effort level');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check RT = f(time)
fig;
hold on;
% mental effort
bar_hdl.Em = jbfill(1:n_bins,...
    (avg_RT.Em_f_time + sem_RT.Em_f_time)',...
    (avg_RT.Em_f_time - sem_RT.Em_f_time)',...
    avg_RT.Em_f_time',...
    Em_col);
% physical effort
bar_hdl.Ep = jbfill(1:n_bins,...
    (avg_RT.Ep_f_time + sem_RT.Ep_f_time)',...
    (avg_RT.Ep_f_time - sem_RT.Ep_f_time)',...
    avg_RT.Ep_f_time',...
    Ep_col);
ylim([1.8 2.5]);
xlim([0 n_bins+1]);
xticks(1:n_bins);
xticklabels(xLabelNames);
xlabel('Trial number');
ylabel('RT (s)');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check effort performance = f(|delta money| level)
fig;
hold on;
for iAbsMoney = 1:(n_R_levels - 1)
    bar_hdl.Em = bar(iAbsMoney-bDist,...
        avg_Eperformance.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),...
        'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iAbsMoney+bDist,...
        avg_Eperformance.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
        'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iAbsMoney-bDist,...
        avg_Eperformance.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),...
        sem_Eperformance.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),...
        'k','LineWidth',lWidth);
    errorbar(iAbsMoney+bDist,...
        avg_Eperformance.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
        sem_Eperformance.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
        'k','LineWidth',lWidth);
end
ylim([90 100]);
ylabel('Effort performance (%)');
xticks(1:3);
xticklabels({'1','2','3'});
xlim([0 n_R_levels]);
xlabel('|Δ money| level');
legend_size(pSize);
ylim([90 100]);
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check effort performance = f(money levels (splitting R and P trials))
fig;
hold on;
for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
    jMoney = abs(iMoney);
    if iMoney < 0
        taskCond = 'P';
    elseif iMoney > 0
        taskCond = 'R';
    end
    bar_hdl.Em = bar(iMoney-bDist,...
        avg_Eperformance.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iMoney+bDist,...
        avg_Eperformance.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iMoney-bDist,...
        avg_Eperformance.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        sem_Eperformance.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
        'k','LineWidth',lWidth)
    errorbar(iMoney+bDist,...
        avg_Eperformance.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        sem_Eperformance.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
        'k','LineWidth',lWidth)
end
xticks([-3:-(1), 1:3]);
xticklabels({'-3','-2','-1','1','2','3'});
xlim([-n_R_levels n_R_levels]);
ylim([80 100]);
ylabel('Effort performance (%)');
xlabel('Money level');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check effort performance = f(E levels)
fig;
hold on;
for iE = 1:(n_E_levels - 1)
    bar_hdl.Em = bar(iE-bDist, avg_Eperformance.perElevel.(['Em_',num2str(iE)]),...
        'FaceColor',Em_col,'BarWidth',bWidth);
    bar_hdl.Ep = bar(iE+bDist, avg_Eperformance.perElevel.(['Ep_',num2str(iE)]),...
        'FaceColor',Ep_col,'BarWidth',bWidth);
    errorbar(iE-bDist,...
        avg_Eperformance.perElevel.(['Em_',num2str(iE)]), sem_Eperformance.perElevel.(['Em_',num2str(iE)]),...
        'k','LineWidth',lWidth);
    errorbar(iE+bDist,...
        avg_Eperformance.perElevel.(['Ep_',num2str(iE)]), sem_Eperformance.perElevel.(['Ep_',num2str(iE)]),...
        'k','LineWidth',lWidth);
end
ylabel('Effort Performance (%)');
xticks(1:3);
xticklabels({'1','2','3'});
xlim([0 n_E_levels]);
ylim([80 100]);
xlabel('Effort level');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

%% check effort performance = f(time)
fig;
hold on;
% mental effort
bar_hdl.Em = jbfill(1:n_bins,...
    (avg_Eperformance.Em_f_time + sem_Eperformance.Em_f_time)',...
    (avg_Eperformance.Em_f_time - sem_Eperformance.Em_f_time)',...
    avg_Eperformance.Em_f_time',...
    Em_col);
% physical effort
bar_hdl.Ep = jbfill(1:n_bins,...
    (avg_Eperformance.Ep_f_time + sem_Eperformance.Ep_f_time)',...
    (avg_Eperformance.Ep_f_time - sem_Eperformance.Ep_f_time)',...
    avg_Eperformance.Ep_f_time',...
    Ep_col);
ylim([90 100]);
xlim([0 n_bins+1]);
xticks(1:n_bins);
xticklabels(xLabelNames);
xlabel('Trial number');
ylabel('Effort Performance (%)');
legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
legend('boxoff');
legend_size(pSize);

end % function