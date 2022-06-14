function[avg_nonDefaultChoice, sem_nonDefaultChoice,...
    avg_conf, sem_conf,...
    avg_RT, sem_RT,...
    avg_Eperformance, sem_Eperformance] = checkChoiceNonDefaultProportion_group_f_metabolites(study_nm, n_bins)
% [avg_nonDefaultChoice, sem_nonDefaultChoice,...
%     avg_conf, sem_conf,...
%     avg_RT, sem_RT,...
%     avg_Eperformance, sem_Eperformance] = checkChoiceNonDefaultProportion_group_f_metabolites(study_nm, n_bins)
%% checkChoiceNonDefaultProportion_group_f_metabolites checks average proportion of choosing 
% the non-default option and the average level of confidence across
% participants depending on the levels of the selected metabolite.
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
if ~exist('study_nm','var') || isempty(study_nm)
%     study_names = {'fMRI_pilots','study1','study2'};
%     study_nm_idx = listdlg('ListString',study_names);
%     study_nm = study_names{study_nm_idx};
study_nm = 'study1';
end
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

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

%% define metabolite and ROI you want to focus on and corresponding subjects
[low_met_subs, high_met_subs, metabolite_nm] = medSplit_metabolites(subject_id);

%% extract data with loop through tasks
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
    
    %% extract the average across participants for each group
    groups = {'low','high'};
    nGrps = length(groups);
    for iGrp = 1:nGrps
        grp_nm = groups{iGrp};
        switch grp_nm
            case 'low'
                subsToInclude = low_met_subs;
            case 'high'
                subsToInclude = high_met_subs;
        end
        avg_nonDefaultChoice.(task_nm).(grp_nm) = mean(nonDefaultChoice_perSub.(task_nm)(:,subsToInclude),2,'omitnan');
        avg_nonDefaultChoice.([task_nm,'_R']).(grp_nm) = mean(nonDefaultChoice_perSub.([task_nm,'_R'])(:,subsToInclude),2,'omitnan');
        avg_nonDefaultChoice.([task_nm,'_P']).(grp_nm) = mean(nonDefaultChoice_perSub.([task_nm,'_P'])(:,subsToInclude),2,'omitnan');
        % SEM
        sem_nonDefaultChoice.(task_nm).(grp_nm) = sem(nonDefaultChoice_perSub.(task_nm)(:,subsToInclude),2);
        sem_nonDefaultChoice.([task_nm,'_R']).(grp_nm) = sem(nonDefaultChoice_perSub.([task_nm,'_R'])(:,subsToInclude),2);
        sem_nonDefaultChoice.([task_nm,'_P']).(grp_nm) = sem(nonDefaultChoice_perSub.([task_nm,'_P'])(:,subsToInclude),2);
        
        % fatigue
        % mean
        avg_nonDefaultChoice.([task_nm,'_f_time']).(grp_nm) = mean(nonDefaultChoice_perSub.([task_nm,'_f_time'])(:,subsToInclude),2,'omitnan');
        avg_conf.([task_nm,'_f_time']).(grp_nm) = mean(conf_perSub.([task_nm,'_f_time'])(:,subsToInclude),2,'omitnan');
        avg_RT.([task_nm,'_f_time']).(grp_nm) = mean(RT_perSub.([task_nm,'_f_time'])(:,subsToInclude),2,'omitnan');
        avg_Eperformance.([task_nm,'_f_time']).(grp_nm) = mean(Eperformance_perSub.([task_nm,'_f_time'])(:,subsToInclude),2,'omitnan');
        % sem
        sem_nonDefaultChoice.([task_nm,'_f_time']).(grp_nm) = sem(nonDefaultChoice_perSub.([task_nm,'_f_time'])(:,subsToInclude),2);
        sem_conf.([task_nm,'_f_time']).(grp_nm) = sem(conf_perSub.([task_nm,'_f_time'])(:,subsToInclude),2);
        sem_RT.([task_nm,'_f_time']).(grp_nm) = sem(RT_perSub.([task_nm,'_f_time'])(:,subsToInclude),2);
        sem_Eperformance.([task_nm,'_f_time']).(grp_nm) = sem(Eperformance_perSub.([task_nm,'_f_time'])(:,subsToInclude),2);
        
        % per task per effort level
        for iE = 1:(n_E_levels - 1)
            % mean
            avg_nonDefaultChoice.perElevel.([task_nm,'_',num2str(iE)]).(grp_nm) = mean(nonDefaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)])(:,subsToInclude),2,'omitnan');
            avg_conf.perElevel.([task_nm,'_',num2str(iE)]).(grp_nm) = mean(conf_perSub.perElevel.([task_nm,'_',num2str(iE)])(:,subsToInclude),2,'omitnan');
            avg_RT.perElevel.([task_nm,'_',num2str(iE)]).(grp_nm) = mean(RT_perSub.perElevel.([task_nm,'_',num2str(iE)])(:,subsToInclude),2,'omitnan');
            avg_Eperformance.perElevel.([task_nm,'_',num2str(iE)]).(grp_nm) = mean(Eperformance_perSub.perElevel.([task_nm,'_',num2str(iE)])(:,subsToInclude),2,'omitnan');
            % SEM
            sem_nonDefaultChoice.perElevel.([task_nm,'_',num2str(iE)]).(grp_nm) = sem(nonDefaultChoice_perSub.perElevel.([task_nm,'_',num2str(iE)])(:,subsToInclude),2);
            sem_conf.perElevel.([task_nm,'_',num2str(iE)]).(grp_nm) = sem(conf_perSub.perElevel.([task_nm,'_',num2str(iE)])(:,subsToInclude),2);
            sem_RT.perElevel.([task_nm,'_',num2str(iE)]).(grp_nm) = sem(RT_perSub.perElevel.([task_nm,'_',num2str(iE)])(:,subsToInclude),2);
            sem_Eperformance.perElevel.([task_nm,'_',num2str(iE)]).(grp_nm) = sem(Eperformance_perSub.perElevel.([task_nm,'_',num2str(iE)])(:,subsToInclude),2);
        end
        
        % per task per money level
        for iAbsMoney = 1:(n_R_levels - 1)
            % mean default choice proportion
            avg_nonDefaultChoice.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]).(grp_nm) = mean(nonDefaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            avg_nonDefaultChoice.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]).(grp_nm) = mean(nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            avg_nonDefaultChoice.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]).(grp_nm) = mean(nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            % mean confidence
            avg_conf.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]).(grp_nm) = mean(conf_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            avg_conf.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]).(grp_nm) = mean(conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            avg_conf.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]).(grp_nm) = mean(conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            % mean RT
            avg_RT.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]).(grp_nm) = mean(RT_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            avg_RT.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]).(grp_nm) = mean(RT_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            avg_RT.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]).(grp_nm) = mean(RT_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            % mean effort performance
            avg_Eperformance.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]).(grp_nm) = mean(Eperformance_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            avg_Eperformance.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]).(grp_nm) = mean(Eperformance_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            avg_Eperformance.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]).(grp_nm) = mean(Eperformance_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(:,subsToInclude),2,'omitnan');
            
            % sem choice proportion
            sem_nonDefaultChoice.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]).(grp_nm) = sem(nonDefaultChoice_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(:,subsToInclude),2);
            sem_nonDefaultChoice.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]).(grp_nm) = sem(nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(:,subsToInclude),2);
            sem_nonDefaultChoice.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]).(grp_nm) = sem(nonDefaultChoice_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(:,subsToInclude),2);
            % sem confidence
            sem_conf.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]).(grp_nm) = sem(conf_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(:,subsToInclude),2);
            sem_conf.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]).(grp_nm) = sem(conf_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(:,subsToInclude),2);
            sem_conf.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]).(grp_nm) = sem(conf_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(:,subsToInclude),2);
            % sem RT
            sem_RT.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]).(grp_nm) = sem(RT_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(:,subsToInclude),2);
            sem_RT.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]).(grp_nm) = sem(RT_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(:,subsToInclude),2);
            sem_RT.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]).(grp_nm) = sem(RT_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(:,subsToInclude),2);
            % sem effort performance
            sem_Eperformance.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)]).(grp_nm) = sem(Eperformance_perSub.perDeltaMoneylevel.([task_nm,'_',num2str(iAbsMoney)])(:,subsToInclude),2);
            sem_Eperformance.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)]).(grp_nm) = sem(Eperformance_perSub.perSignedMoneylevel.([task_nm,'_R_',num2str(iAbsMoney)])(:,subsToInclude),2);
            sem_Eperformance.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)]).(grp_nm) = sem(Eperformance_perSub.perSignedMoneylevel.([task_nm,'_P_',num2str(iAbsMoney)])(:,subsToInclude),2);
        end % absolute money level loo
    end % low/high metabolite group
    
end % loop through physical and mental effort

%% figures
% figure parameters
lWidth_50percentTrait = 2;
lWidth = 3;
pSize = 40;
bWidth = 0.4;
bDist = 0.2;
Em_col_low = [255 255 204]./255;
Em_col_high = [161 218 180]./255;
Ep_col_low = [65 182 196]./255;
Ep_col_high = [34 94 168]./255;

% check choices = f(fatigue, R/P, level of R, level of E)
%% check choices = f(fatigue)
fig;
% mark the 50% trait
plot(1:n_bins, 50*ones(1,n_bins),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
% mental effort
bar_hdl_Em_low = jbfill(1:n_bins,...
    (avg_nonDefaultChoice.Em_f_time.low + sem_nonDefaultChoice.Em_f_time.low)'.*100,...
    (avg_nonDefaultChoice.Em_f_time.low - sem_nonDefaultChoice.Em_f_time.low)'.*100,...
    avg_nonDefaultChoice.Em_f_time.low'.*100,...
    Em_col_low);
bar_hdl_Em_high = jbfill(1:n_bins,...
    (avg_nonDefaultChoice.Em_f_time.high + sem_nonDefaultChoice.Em_f_time.high)'.*100,...
    (avg_nonDefaultChoice.Em_f_time.high - sem_nonDefaultChoice.Em_f_time.high)'.*100,...
    avg_nonDefaultChoice.Em_f_time.high'.*100,...
    Em_col_high);
% physical effort
bar_hdl_Ep_low = jbfill(1:n_bins,...
    (avg_nonDefaultChoice.Ep_f_time.low + sem_nonDefaultChoice.Ep_f_time.low)'.*100,...
    (avg_nonDefaultChoice.Ep_f_time.low - sem_nonDefaultChoice.Ep_f_time.low)'.*100,...
    avg_nonDefaultChoice.Ep_f_time.low'.*100,...
    Ep_col_low);
bar_hdl_Ep_high = jbfill(1:n_bins,...
    (avg_nonDefaultChoice.Ep_f_time.high + sem_nonDefaultChoice.Ep_f_time.high)'.*100,...
    (avg_nonDefaultChoice.Ep_f_time.high - sem_nonDefaultChoice.Ep_f_time.high)'.*100,...
    avg_nonDefaultChoice.Ep_f_time.high'.*100,...
    Ep_col_high);
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
ylabel('Choice non-default option (%)');
legend([bar_hdl_Em_low, bar_hdl_Em_high,...
    bar_hdl_Ep_low, bar_hdl_Ep_high],...
    {'mental low','mental high','physical low','physical high'});
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
bar_hdl.Em_low = bar(1-bDist, avg_nonDefaultChoice.Em_R.low,...
    'FaceColor',Em_col_low,'BarWidth',bWidth);
errorbar(1-bDist, avg_nonDefaultChoice.Em_R.low, sem_nonDefaultChoice.Em_R.low,...
    'k','LineWidth',lWidth);
bar_hdl.Em_high = bar(1+bDist, avg_nonDefaultChoice.Em_R.high,...
    'FaceColor',Em_col_high,'BarWidth',bWidth);
errorbar(1+bDist, avg_nonDefaultChoice.Em_R.high, sem_nonDefaultChoice.Em_R.high,...
    'k','LineWidth',lWidth);
% P mental effort
bar(2-bDist, avg_nonDefaultChoice.Em_P.low,...
    'FaceColor',Em_col_low,'BarWidth',bWidth);
errorbar(2-bDist, avg_nonDefaultChoice.Em_P.low, sem_nonDefaultChoice.Em_P.low,...
    'k','LineWidth',lWidth);
bar(2+bDist, avg_nonDefaultChoice.Em_P.high,...
    'FaceColor',Em_col_high,'BarWidth',bWidth);
errorbar(2+bDist, avg_nonDefaultChoice.Em_P.high, sem_nonDefaultChoice.Em_P.high,...
    'k','LineWidth',lWidth);
% R physical effort
bar_hdl.Ep_low = bar(3-bDist, avg_nonDefaultChoice.Ep_R.low,...
    'FaceColor',Ep_col_low,'BarWidth',bWidth);
errorbar(3-bDist, avg_nonDefaultChoice.Ep_R.low, sem_nonDefaultChoice.Ep_R.low,...
    'k','LineWidth',lWidth);
bar_hdl.Ep_high = bar(3+bDist, avg_nonDefaultChoice.Ep_R.high,...
    'FaceColor',Ep_col_high,'BarWidth',bWidth);
errorbar(3+bDist, avg_nonDefaultChoice.Ep_R.high, sem_nonDefaultChoice.Ep_R.high,...
    'k','LineWidth',lWidth);
% P physical effort
bar(4-bDist, avg_nonDefaultChoice.Ep_P.low,...
    'FaceColor',Ep_col_low,'BarWidth',bWidth);
errorbar(4-bDist, avg_nonDefaultChoice.Ep_P.low, sem_nonDefaultChoice.Ep_P.low,...
    'k','LineWidth',lWidth);
bar(4+bDist, avg_nonDefaultChoice.Ep_P.high,...
    'FaceColor',Ep_col_high,'BarWidth',bWidth);
errorbar(4+bDist, avg_nonDefaultChoice.Ep_P.high, sem_nonDefaultChoice.Ep_P.high,...
    'k','LineWidth',lWidth);
xticks([1-bDist, 1+bDist, 2-bDist, 2+bDist, 3-bDist, 3+bDist, 4-bDist, 4+bDist]);
xticklabels({'R low','R high','P low','P high','R low','R high','P low','P high'});
ylim([0 100]);
xlim([0 5]);
ylabel('Choice non-default option (%)');
legend_size(pSize);

%% check choices = f(|delta money| level)
fig;
% mark the 50% trait
plot(0:n_bins, 50*ones(1,n_bins+1),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
for iAbsMoney = 1:(n_R_levels - 1)
    bar_hdl.Em_low = bar(iAbsMoney-bDist*3/2,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]).low,...
        'FaceColor',Em_col_low,'BarWidth',bWidth/2);
    bar_hdl.Em_high = bar(iAbsMoney-bDist/2,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]).high,...
        'FaceColor',Em_col_high,'BarWidth',bWidth/2);
    bar_hdl.Ep_low = bar(iAbsMoney+bDist/2,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]).low,...
        'FaceColor',Ep_col_low,'BarWidth',bWidth/2);
    bar_hdl.Ep_high = bar(iAbsMoney+bDist*3/2,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]).high,...
        'FaceColor',Ep_col_high,'BarWidth',bWidth/2);
    errorbar(iAbsMoney-bDist*3/2,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]).low,...
        sem_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]).low,...
        'k','LineWidth',lWidth);
    errorbar(iAbsMoney-bDist/2,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]).high,...
        sem_nonDefaultChoice.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]).high,...
        'k','LineWidth',lWidth);
    errorbar(iAbsMoney+bDist/2,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]).low,...
        sem_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]).low,...
        'k','LineWidth',lWidth);
    errorbar(iAbsMoney+bDist*3/2,...
        avg_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]).high,...
        sem_nonDefaultChoice.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]).high,...
        'k','LineWidth',lWidth);
end
ylim([0 100]);
xticks(1:3);
xticklabels({'1','2','3'});
xlim([0 n_R_levels]);
xlabel('|Δ money| level');
legend_size(pSize);
ylim([0 100]);
ylabel('Choice non-default option (%)');
legend([bar_hdl.Em_low, bar_hdl.Em_high, bar_hdl.Ep_low, bar_hdl.Ep_high],...
    {['mental l',metabolite_nm],['mental h',metabolite_nm],...
    ['physical l',metabolite_nm],['physical h',metabolite_nm],},...
    'Location','NorthWest');
legend('boxoff');
legend_size(pSize);

% %% check choices = f(money levels (splitting R and P trials))
% fig;
% % mark the 50% trait
% plot(-n_bins:n_bins, 50*ones(1,(n_bins*2)+1),...
%     'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
% hold on;
% for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
%     jMoney = abs(iMoney);
%     if iMoney < 0
%         taskCond = 'P';
%     elseif iMoney > 0
%         taskCond = 'R';
%     end
%     bar_hdl.Em = bar(iMoney-bDist,...
%         avg_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         'FaceColor',Em_col,'BarWidth',bWidth);
%     bar_hdl.Ep = bar(iMoney+bDist,...
%         avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         'FaceColor',Ep_col,'BarWidth',bWidth);
%     errorbar(iMoney-bDist,...
%         avg_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         sem_nonDefaultChoice.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         'k','LineWidth',lWidth);
%     errorbar(iMoney+bDist,...
%         avg_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         sem_nonDefaultChoice.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         'k','LineWidth',lWidth);
% end
% ylim([0 100]);
% ylabel('Choice non-default option (%)');
% xticks([-3:-(1), 1:3]);
% xticklabels({'-3','-2','-1','1','2','3'});
% xlim([-n_R_levels n_R_levels]);
% xlabel('Money level');
% legend_size(pSize);
% ylim([0 100]);
% ylabel('Choice non-default option (%)');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);

%% check choices = f(E level)
fig;
% mark the 50% trait
plot(0:n_bins, 50*ones(1,n_bins+1),...
    'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
hold on;
for iE = 1:(n_E_levels - 1)
    bar_hdl.Em_low = bar(iE-bDist*3/2, avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]).low,...
        'FaceColor',Em_col_low,'BarWidth',bWidth/2);
    bar_hdl.Em_high = bar(iE-bDist/2, avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]).high,...
        'FaceColor',Em_col_high,'BarWidth',bWidth/2);
    bar_hdl.Ep_low = bar(iE+bDist/2, avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]).low,...
        'FaceColor',Ep_col_low,'BarWidth',bWidth/2);
    bar_hdl.Ep_high = bar(iE+bDist*3/2, avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]).high,...
        'FaceColor',Ep_col_high,'BarWidth',bWidth/2);
    errorbar(iE-bDist*3/2,...
        avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]).low,...
        sem_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]).low,...
        'k','LineWidth',lWidth);
    errorbar(iE-bDist/2,...
        avg_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]).high,...
        sem_nonDefaultChoice.perElevel.(['Em_',num2str(iE)]).high,...
        'k','LineWidth',lWidth);
    errorbar(iE+bDist/2,...
        avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]).low,...
        sem_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]).low,...
        'k','LineWidth',lWidth);
    errorbar(iE+bDist*3/2,...
        avg_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]).high,...
        sem_nonDefaultChoice.perElevel.(['Ep_',num2str(iE)]).high,...
        'k','LineWidth',lWidth);
end
ylim([0 100]);
ylabel('Choice non-default option (%)');
xticks(1:3);
xticklabels({'1','2','3'});
xlim([0 n_E_levels]);
xlabel('Effort level');
legend([bar_hdl.Em_low, bar_hdl.Em_high,...
    bar_hdl.Ep_low, bar_hdl.Ep_high],...
    {['mental l',metabolite_nm],['mental h',metabolite_nm],...
    ['physical l',metabolite_nm],['physical h',metabolite_nm]});
legend('boxoff');
legend_size(pSize);

% %% check confidence = f(money levels (splitting R and P trials))
% fig;
% % mark the 50% trait
% plot(-n_bins:n_bins, 50*ones(1,(n_bins*2)+1),...
%     'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
% hold on;
% for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
%     jMoney = abs(iMoney);
%     if iMoney < 0
%         taskCond = 'P';
%     elseif iMoney > 0
%         taskCond = 'R';
%     end
%     bar_hdl.Em = bar(iMoney-bDist, avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]), 'FaceColor',Em_col,'BarWidth',bWidth);
%     bar_hdl.Ep = bar(iMoney+bDist, avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]), 'FaceColor',Ep_col,'BarWidth',bWidth);
%     errorbar(iMoney-bDist,...
%         avg_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         sem_conf.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         'k','LineWidth',lWidth);
%     errorbar(iMoney+bDist,...
%         avg_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         sem_conf.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         'k','LineWidth',lWidth);
% end
% xticks([-3:-(1), 1:3]);
% xticklabels({'-3','-2','-1','1','2','3'});
% xlim([-n_R_levels n_R_levels]);
% ylim([0 100]);
% ylabel('Level of confidence');
% xlabel('Money level');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);
% 
% %% check confidence = f(E levels)
% fig;
% % mark the 50% trait
% plot(0:n_bins, 50*ones(1,n_bins+1),...
%     'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
% hold on;
% for iE = 1:(n_E_levels - 1)
%     bar_hdl.Em = bar(iE-bDist,...
%         avg_conf.perElevel.(['Em_',num2str(iE)]),...
%         'FaceColor',Em_col,'BarWidth',bWidth);
%     bar_hdl.Ep = bar(iE+bDist,...
%         avg_conf.perElevel.(['Ep_',num2str(iE)]),...
%         'FaceColor',Ep_col,'BarWidth',bWidth);
%     errorbar(iE-bDist,...
%         avg_conf.perElevel.(['Em_',num2str(iE)]),...
%         sem_conf.perElevel.(['Em_',num2str(iE)]),...
%         'k','LineWidth',lWidth);
%     errorbar(iE+bDist,...
%         avg_conf.perElevel.(['Ep_',num2str(iE)]),...
%         sem_conf.perElevel.(['Ep_',num2str(iE)]),...
%         'k','LineWidth',lWidth);
% end
% ylim([0 100]);
% ylabel('Confidence');
% xticks(1:3);
% xticklabels({'1','2','3'});
% xlim([0 n_E_levels]);
% xlabel('Effort level');
% legend_size(pSize);
% ylabel('Level of confidence');
% xlabel('Effort level');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);
% 
% %% check confidence = f(time)
% fig;
% % mark the 50% trait
% plot(1:n_bins, 50*ones(1,n_bins),...
%     'LineWidth',lWidth_50percentTrait,'Color','k','LineStyle',':');
% hold on;
% % mental effort
% bar_hdl.Em = jbfill(1:n_bins,...
%     (avg_conf.Em_f_time + sem_conf.Em_f_time)'.*100,...
%     (avg_conf.Em_f_time - sem_conf.Em_f_time)'.*100,...
%     avg_conf.Em_f_time'.*100,...
%     Em_col);
% % physical effort
% bar_hdl.Ep = jbfill(1:n_bins,...
%     (avg_conf.Ep_f_time + sem_conf.Ep_f_time)'.*100,...
%     (avg_conf.Ep_f_time - sem_conf.Ep_f_time)'.*100,...
%     avg_conf.Ep_f_time'.*100,...
%     Ep_col);
% ylim([0 100]);
% xlim([0 n_bins+1]);
% xlabel('trial bins');
% ylabel('Confidence (%)');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);
% 
% %% check RT = f(money levels (splitting R and P trials))
% fig;
% hold on;
% for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
%     jMoney = abs(iMoney);
%     if iMoney < 0
%         taskCond = 'P';
%     elseif iMoney > 0
%         taskCond = 'R';
%     end
%     bar_hdl.Em = bar(iMoney-bDist,...
%         avg_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         'FaceColor',Em_col,'BarWidth',bWidth);
%     bar_hdl.Ep = bar(iMoney+bDist,...
%         avg_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         'FaceColor',Ep_col,'BarWidth',bWidth);
%     errorbar(iMoney-bDist,...
%         avg_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         sem_RT.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         'k','LineWidth',lWidth)
%     errorbar(iMoney+bDist,...
%         avg_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         sem_RT.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         'k','LineWidth',lWidth)
% end
% xticks([-3:-(1), 1:3]);
% xticklabels({'-3','-2','-1','1','2','3'});
% xlim([-n_R_levels n_R_levels]);
% ylim([1.5 3]);
% ylabel('RT (s)');
% xlabel('Money level');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);
% 
% %% check RT = f(E levels)
% fig;
% hold on;
% for iE = 1:(n_E_levels - 1)
%     bar_hdl.Em = bar(iE-bDist, avg_RT.perElevel.(['Em_',num2str(iE)]),...
%         'FaceColor',Em_col,'BarWidth',bWidth);
%     bar_hdl.Ep = bar(iE+bDist, avg_RT.perElevel.(['Ep_',num2str(iE)]),...
%         'FaceColor',Ep_col,'BarWidth',bWidth);
%     errorbar(iE-bDist,...
%         avg_RT.perElevel.(['Em_',num2str(iE)]), sem_RT.perElevel.(['Em_',num2str(iE)]),...
%         'k','LineWidth',lWidth);
%     errorbar(iE+bDist,...
%         avg_RT.perElevel.(['Ep_',num2str(iE)]), sem_RT.perElevel.(['Ep_',num2str(iE)]),...
%         'k','LineWidth',lWidth);
% end
% ylim([0 5]);
% ylabel('RT (s)');
% xticks(1:3);
% xticklabels({'1','2','3'});
% xlim([0 n_E_levels]);
% ylim([1.5 3]);
% xlabel('Effort level');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);
% 
% %% check RT = f(time)
% fig;
% hold on;
% % mental effort
% bar_hdl.Em = jbfill(1:n_bins,...
%     (avg_RT.Em_f_time + sem_RT.Em_f_time)',...
%     (avg_RT.Em_f_time - sem_RT.Em_f_time)',...
%     avg_RT.Em_f_time',...
%     Em_col);
% % physical effort
% bar_hdl.Ep = jbfill(1:n_bins,...
%     (avg_RT.Ep_f_time + sem_RT.Ep_f_time)',...
%     (avg_RT.Ep_f_time - sem_RT.Ep_f_time)',...
%     avg_RT.Ep_f_time',...
%     Ep_col);
% ylim([1.5 3]);
% xlim([0 n_bins+1]);
% xlabel('trial bins');
% ylabel('RT (s)');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);
% 
% %% check effort performance = f(|delta money| level)
% fig;
% hold on;
% for iAbsMoney = 1:(n_R_levels - 1)
%     bar_hdl.Em = bar(iAbsMoney-bDist,...
%         avg_Eperformance.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),...
%         'FaceColor',Em_col,'BarWidth',bWidth);
%     bar_hdl.Ep = bar(iAbsMoney+bDist,...
%         avg_Eperformance.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
%         'FaceColor',Ep_col,'BarWidth',bWidth);
%     errorbar(iAbsMoney-bDist,...
%         avg_Eperformance.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),...
%         sem_Eperformance.perDeltaMoneylevel.(['Em_',num2str(iAbsMoney)]),...
%         'k','LineWidth',lWidth);
%     errorbar(iAbsMoney+bDist,...
%         avg_Eperformance.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
%         sem_Eperformance.perDeltaMoneylevel.(['Ep_',num2str(iAbsMoney)]),...
%         'k','LineWidth',lWidth);
% end
% ylim([90 100]);
% ylabel('Effort performance (%)');
% xticks(1:3);
% xticklabels({'1','2','3'});
% xlim([0 n_R_levels]);
% xlabel('|Δ money| level');
% legend_size(pSize);
% ylim([90 100]);
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);
% 
% %% check effort performance = f(money levels (splitting R and P trials))
% fig;
% hold on;
% for iMoney = [-(n_R_levels-1):(-1), 1:(n_R_levels - 1)]
%     jMoney = abs(iMoney);
%     if iMoney < 0
%         taskCond = 'P';
%     elseif iMoney > 0
%         taskCond = 'R';
%     end
%     bar_hdl.Em = bar(iMoney-bDist,...
%         avg_Eperformance.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         'FaceColor',Em_col,'BarWidth',bWidth);
%     bar_hdl.Ep = bar(iMoney+bDist,...
%         avg_Eperformance.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         'FaceColor',Ep_col,'BarWidth',bWidth);
%     errorbar(iMoney-bDist,...
%         avg_Eperformance.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         sem_Eperformance.perSignedMoneylevel.(['Em_',taskCond','_',num2str(jMoney)]),...
%         'k','LineWidth',lWidth)
%     errorbar(iMoney+bDist,...
%         avg_Eperformance.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         sem_Eperformance.perSignedMoneylevel.(['Ep_',taskCond','_',num2str(jMoney)]),...
%         'k','LineWidth',lWidth)
% end
% xticks([-3:-(1), 1:3]);
% xticklabels({'-3','-2','-1','1','2','3'});
% xlim([-n_R_levels n_R_levels]);
% ylim([80 100]);
% ylabel('Effort performance (%)');
% xlabel('Money level');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);
% 
% %% check effort performance = f(E levels)
% fig;
% hold on;
% for iE = 1:(n_E_levels - 1)
%     bar_hdl.Em = bar(iE-bDist, avg_Eperformance.perElevel.(['Em_',num2str(iE)]),...
%         'FaceColor',Em_col,'BarWidth',bWidth);
%     bar_hdl.Ep = bar(iE+bDist, avg_Eperformance.perElevel.(['Ep_',num2str(iE)]),...
%         'FaceColor',Ep_col,'BarWidth',bWidth);
%     errorbar(iE-bDist,...
%         avg_Eperformance.perElevel.(['Em_',num2str(iE)]), sem_Eperformance.perElevel.(['Em_',num2str(iE)]),...
%         'k','LineWidth',lWidth);
%     errorbar(iE+bDist,...
%         avg_Eperformance.perElevel.(['Ep_',num2str(iE)]), sem_Eperformance.perElevel.(['Ep_',num2str(iE)]),...
%         'k','LineWidth',lWidth);
% end
% ylabel('Effort Performance (%)');
% xticks(1:3);
% xticklabels({'1','2','3'});
% xlim([0 n_E_levels]);
% ylim([80 100]);
% xlabel('Effort level');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);
% 
% %% check effort performance = f(time)
% fig;
% hold on;
% % mental effort
% bar_hdl.Em = jbfill(1:n_bins,...
%     (avg_Eperformance.Em_f_time + sem_Eperformance.Em_f_time)',...
%     (avg_Eperformance.Em_f_time - sem_Eperformance.Em_f_time)',...
%     avg_Eperformance.Em_f_time',...
%     Em_col);
% % physical effort
% bar_hdl.Ep = jbfill(1:n_bins,...
%     (avg_Eperformance.Ep_f_time + sem_Eperformance.Ep_f_time)',...
%     (avg_Eperformance.Ep_f_time - sem_Eperformance.Ep_f_time)',...
%     avg_Eperformance.Ep_f_time',...
%     Ep_col);
% ylim([80 100]);
% xlim([0 n_bins+1]);
% xlabel('trial bins');
% ylabel('Effort Performance (%)');
% legend([bar_hdl.Em, bar_hdl.Ep],'mental','physical');
% legend('boxoff');
% legend_size(pSize);

end % function