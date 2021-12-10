function[betas, pvalues] = logitfit_choices_group(computerRoot, study_nm, figDispGroup, figDispIndiv, n_NV_bins)
% [betas, pvalues] = logitfit_choices_group(computerRoot, study_nm, figDisp, figDispIndiv, n_NV_bins)
%
% INPUTS
% computerRoot: pathway where data is
%
% study_nm: study name
%
% figDisp: display group figures (1) or not (0)
%
% figDispIndiv: display individual figures (1) or not (0)
%
% n_NV_bins: number of bins for net value
%
% OUTPUT
% betas: structure with betas
%
% pvalues: structure with p.values

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% study names
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% by default, display group figures
if ~exist('figDispGroup','var') || isempty(figDispGroup)
    figDispGroup = 1;
    disp(['figDispGroup was not defined in the inputs so that by default ',...
        'group figures are displayed.']);
end

%% by default, do not display individual figures
if ~exist('figDispIndiv','var') || isempty(figDispIndiv)
    figDispIndiv = 0;
    disp(['figDispGroup was not defined in the inputs so that by default ',...
        'individual figures are not displayed.']);
end

%% if not defined in the inputs, define by default some number of bins for the graphs
if ~exist('n_NV_bins','var') || isempty(n_NV_bins)
    n_NV_bins = 6;
end
%% define R/P/E levels
P_levels = [-3, -2, -1];
nPlevels = length(P_levels);
R_levels = [1, 2, 3];
nRlevels = length(R_levels);
money_levels = [-3, -2, -1, 1, 2, 3];
nMoneyLevels = length(money_levels);
E_levels = [1, 2, 3];
nELevels = length(E_levels);

%% working directories
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];
resultFolder_a = [studyBehaviorFolder,'results',filesep];
if ~exist(resultFolder_a,'dir')
    mkdir(resultFolder_a);
end
resultFolder = [resultFolder_a,'figures',filesep];
if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end

%% subject selection
[subject_id, NS] = LGCM_subject_selection(study_nm);

%% initialize variables of interest
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
        case 2
            task_id = 'Em';
    end
    [betas.(task_id).Mdl1.kb0,...
        betas.(task_id).Mdl1.kMoney,...
        betas.(task_id).Mdl1.kEffort] = deal(NaN(1,NS));
    
    [choiceNonDef.perMoneyLevel.(task_id),...
        choiceFitNonDef.perMoneyLevel.Mdl1.(task_id),...
        choiceFitNonDef.perMoneyLevel.Mdl2.(task_id),...
        choiceFitNonDef.perMoneyLevel.Mdl3.(task_id),...
        choiceFitNonDef.perMoneyLevel.Mdl4.(task_id)] = deal(NaN(nMoneyLevels,NS));
    [choiceNonDef.perEffortLevel.(task_id),...
        choiceFitNonDef.perEffortLevel.Mdl1.(task_id),...
        choiceFitNonDef.perEffortLevel.Mdl2.(task_id),...
        choiceFitNonDef.perEffortLevel.Mdl3.(task_id),...
        choiceFitNonDef.perEffortLevel.Mdl4.(task_id)] = deal(NaN(nELevels,NS));
    
    [choiceNonDef.perNVLevel.Mdl1.(task_id),...
        choiceFitNonDef.perNVLevel.Mdl1.(task_id),...
        NV_bins.Mdl1.(task_id)] = deal(NaN(n_NV_bins, NS));
end % physical/mental loop

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    % load individual data
    [betas_tmp, choices_tmp] = logitfit_choices(computerRoot, study_nm, sub_nm, figDispIndiv, n_NV_bins);
    
    % pool data across subjects
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
            case 2
                task_id = 'Em';
        end
        % extract betas
        betas.(task_id).Mdl1.kb0(iS) = betas_tmp.(task_id).Mdl1.kb0;
        betas.(task_id).Mdl1.kMoney(iS) = betas_tmp.(task_id).Mdl1.kMoney;
        betas.(task_id).Mdl1.kEffort(iS) = betas_tmp.(task_id).Mdl1.kEffort;
        % extract choices
        choiceNonDef.perMoneyLevel.(task_id)(:,iS) = choices_tmp.choiceNonDef.perMoneyLevel.(task_id);
        choiceNonDef.perEffortLevel.(task_id)(:,iS) = choices_tmp.choiceNonDef.perEffortLevel.(task_id);
        % extract fit
        choiceNonDef.perNVLevel.Mdl1.(task_id)(:,iS) = choices_tmp.choiceNonDef.perNVLevel.Mdl1.(task_id);
        choiceFitNonDef.perNVLevel.Mdl1.(task_id)(:,iS) = choices_tmp.choiceFitNonDef.perNVLevel.Mdl1.(task_id);
        NV_bins.Mdl1.(task_id)(:,iS) = choices_tmp.NV_bins.Mdl1.(task_id);
        choiceFitNonDef.perMoneyLevel.Mdl1.(task_id)(:,iS) = choices_tmp.choiceFitNonDef.perMoneyLevel.Mdl1.(task_id);
        choiceFitNonDef.perEffortLevel.Mdl1.(task_id)(:,iS) = choices_tmp.choiceFitNonDef.perEffortLevel.Mdl1.(task_id);
    end % physical/mental loop
end % subject loop

%% average data
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
            task_fullName = 'physical';
        case 2
            task_id = 'Em';
            task_fullName = 'mental';
    end
    
    [m_choiceNonDef.perMoneyLevel.(task_id),...
        sem_choiceNonDef.perMoneyLevel.(task_id)] = mean_sem_sd(choiceNonDef.perMoneyLevel.(task_id), 2);
    [m_choiceFitNonDef.perMoneyLevel.Mdl1.(task_id),...
        sem_choiceFitNonDef.perMoneyLevel.Mdl1.(task_id)] = mean_sem_sd(choiceFitNonDef.perMoneyLevel.Mdl1.(task_id), 2);
    [m_choiceNonDef.perEffortLevel.(task_id),...
        sem_choiceNonDef.perEffortLevel.(task_id)] = mean_sem_sd(choiceNonDef.perEffortLevel.(task_id), 2);
    [m_choiceFitNonDef.perEffortLevel.Mdl1.(task_id),...
        sem_choiceFitNonDef.perEffortLevel.Mdl1.(task_id)] = mean_sem_sd(choiceFitNonDef.perEffortLevel.Mdl1.(task_id), 2);
    [m_NV_bins.Mdl1.(task_id),...
        sem_NV_bins.Mdl1.(task_id)] = mean_sem_sd(NV_bins.Mdl1.(task_id), 2);
    [m_choiceNonDef.perNVLevel.Mdl1.(task_id),...
        sem_choiceNonDef.perNVLevel.Mdl1.(task_id)] = mean_sem_sd(choiceNonDef.perNVLevel.Mdl1.(task_id), 2);
    [m_choiceFitNonDef.perNVLevel.Mdl1.(task_id),...
        sem_choiceFitNonDef.perNVLevel.Mdl1.(task_id)] = mean_sem_sd(choiceFitNonDef.perNVLevel.Mdl1.(task_id), 2);
    % average betas
    [betas.mean.(task_id).Mdl1.kb0,...
        betas.sem.(task_id).Mdl1.kb0,...
        betas.sd.(task_id).Mdl1.kb0] = mean_sem_sd(betas.(task_id).Mdl1.kb0,2);
    [betas.mean.(task_id).Mdl1.kMoney,...
        betas.sem.(task_id).Mdl1.kMoney,...
        betas.sd.(task_id).Mdl1.kMoney] = mean_sem_sd(betas.(task_id).Mdl1.kMoney,2);
    [betas.mean.(task_id).Mdl1.kEffort,...
        betas.sem.(task_id).Mdl1.kEffort,...
        betas.sd.(task_id).Mdl1.kEffort] = mean_sem_sd(betas.(task_id).Mdl1.kEffort,2);
    
    %% test how significant betas are
    [~,pvalues.(task_id).Mdl1.kb0] = ttest(betas.(task_id).Mdl1.kb0);
    [~,pvalues.(task_id).Mdl1.kMoney] = ttest(betas.(task_id).Mdl1.kMoney);
    [~,pvalues.(task_id).Mdl1.kEffort] = ttest(betas.(task_id).Mdl1.kEffort);
    
    %% display average data
    if figDispGroup == 1
        pSize = 30;
        lWidth = 3;
        lWidth_borders = 1;
        %% model 1
        %% display choice = f(net value)
        fig;
        %         pointMdl1 = errorbar(m_NV_bins.Mdl1.(task_id),...
        %             m_choiceNonDef.perNVLevel.Mdl1.(task_id),...
        %             sem_choiceNonDef.perNVLevel.Mdl1.(task_id),...
        %             sem_choiceNonDef.perNVLevel.Mdl1.(task_id),...
        %             sem_NV_bins.Mdl1.(task_id),...
        %             sem_NV_bins.Mdl1.(task_id));
        pointMdl1 = errorbar(m_NV_bins.Mdl1.(task_id),...
            m_choiceNonDef.perNVLevel.Mdl1.(task_id),...
            sem_choiceNonDef.perNVLevel.Mdl1.(task_id));
        pointMdl1.Color = [0 0 0];
        pointMdl1.Marker = 'o';
        pointMdl1.LineStyle = 'none';
        pointMdl1.LineWidth = lWidth;
        hold on;
        line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
        line(xlim(),[1 1],'LineWidth',lWidth_borders,'Color',[0 0 0]);
        lHdlMdl1 = plot(m_NV_bins.Mdl1.(task_id), m_choiceFitNonDef.perNVLevel.Mdl1.(task_id));
        lHdlMdl1.LineStyle = '--';
        lHdlMdl1.LineWidth = lWidth;
        lHdlMdl1.Color = [143 0 0]./255;
        ylim([-0.2 1.2]);
        xlabel(['Net value ',task_fullName]);
        ylabel('Choice non-default option (%)');
        legend_size(pSize);
        %         saveas(resultFolder)
        
        %% choice non-default = f(money levels)
        fig;
        pointMdl1 = errorbar(money_levels,...
            m_choiceNonDef.perMoneyLevel.(task_id),...
            sem_choiceNonDef.perMoneyLevel.(task_id));
        pointMdl1.Color = [0 0 0];
        pointMdl1.Marker = 'o';
        pointMdl1.LineStyle = 'none';
        pointMdl1.LineWidth = lWidth;
        hold on;
        line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
        line(xlim(),[1 1],'LineWidth',lWidth_borders,'Color',[0 0 0]);
        lHdlMdl1 = plot(money_levels, m_choiceFitNonDef.perMoneyLevel.Mdl1.(task_id));
        lHdlMdl1.LineStyle = '--';
        lHdlMdl1.LineWidth = lWidth;
        lHdlMdl1.Color = [143 0 0]./255;
        ylim([-0.2 1.2]);
        xlabel([task_fullName,' Money level']);
        ylabel('Choice non-default option (%)');
        legend_size(pSize);
%         saveas(resultFolder)
        
        %% choice non-default = f(effort levels)
        fig;
        pointMdl1 = errorbar(E_levels,...
            m_choiceNonDef.perEffortLevel.(task_id),...
            sem_choiceNonDef.perEffortLevel.(task_id));
        pointMdl1.Color = [0 0 0];
        pointMdl1.Marker = 'o';
        pointMdl1.LineStyle = 'none';
        pointMdl1.LineWidth = lWidth;
        hold on;
        line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
        line(xlim(),[1 1],'LineWidth',lWidth_borders,'Color',[0 0 0]);
        lHdlMdl1 = plot(E_levels,...
            m_choiceFitNonDef.perEffortLevel.Mdl1.(task_id));
        lHdlMdl1.LineStyle = '--';
        lHdlMdl1.LineWidth = lWidth;
        lHdlMdl1.Color = [143 0 0]./255;
        ylim([-0.2 1.2]);
        xlabel([task_fullName,' effort level']);
        ylabel('Choice non-default option (%)');
        legend_size(pSize);
%         saveas(resultFolder)
        
    end % display figure
    
end % physical/mental loop

end % function