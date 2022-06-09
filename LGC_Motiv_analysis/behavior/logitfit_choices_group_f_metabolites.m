function [betas, pvalues] = logitfit_choices_group_f_metabolites()
% [betas, pvalues] = logitfit_choices_group_f_metabolites()
%
% INPUTS
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
%% by default, display monetary levels instead of actual monetary amounts
if ~exist('dispMoneyOrLevels','var') || isempty(dispMoneyOrLevels)
    dispMoneyOrLevels = 'levels';
    %     dispMoneyOrLevels = 'money';
end
%% if not defined in the inputs, define by default some number of bins for the graphs
if ~exist('n_NV_bins','var') || isempty(n_NV_bins)
    n_NV_bins = 6;
end
%% if not defined in the inputs, define by default some number of bins for the graphs
if ~exist('n_trialN_bins','var') || isempty(n_trialN_bins)
    n_trialN_bins = 6;
end
%% define R/P/E levels
money_levels = [-3, -2, -1, 1, 2, 3];
nMoneyLevels = length(money_levels);
E_levels = [1, 2, 3];
nELevels = length(E_levels);
nMdl = 4;

%% define metabolite and ROI you want to focus on
% ROI
ROIs = {'dmPFC','aIns'};
nROIs = length(ROIs);
ROI_idx = spm_input('Metabolites in which brain area?',1,'m',...
    ROIs,1:nROIs,0);
ROI_nm = ROIs{ROI_idx};
% select metabolite of interest
metabolites = {'Mac','Ala','Asp','PCho','Cr','PCr','GABA',...
    'Gln','Glu','GSH','Gly','Ins','Lac','NAA','Scyllo','Tau',...
    'Asc','Glc','NAAG','GPC','PE','Ser',...
    'NAA_NAAG','Glu_Gln','GPC_PCho','Cr_PCr','Gly_Ins','Gln_div_Glu'};
n_met = length(metabolites);
metabolite_idx = spm_input('Which metabolite to focus on?',1,'m',...
    metabolites,1:n_met,0);
metabolite_nm = metabolites{metabolite_idx};

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
[condition] = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);
% store subject list to know which beta corresponds to which subject
betas.subList = subject_id;

%% extract all metabolites
[metabolites] = metabolite_load(subject_id);
% focus on metabolite and brain area selected
metabolite_allSubs = metabolites.(ROI_nm).(metabolite_nm);
% perform a median split based on the metabolite selected in the ROI
% selected
med_metabolite_allSubs = median(metabolite_allSubs,'omitnan');
% extract index of participants with low or high level of metabolites
low_met_subs = metabolite_allSubs <= med_metabolite_allSubs;
high_met_subs = metabolite_allSubs > med_metabolite_allSubs;

%% initialize variables of interest
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
        case 2
            task_id = 'Em';
    end
    [betas.(task_id).mdl_1.kMoney,...
        betas.(task_id).mdl_1.kEffort,...
        betas.(task_id).mdl_2.kR,...
        betas.(task_id).mdl_2.kP,...
        betas.(task_id).mdl_2.kEffort,...
        betas.(task_id).mdl_3.kMoney,...
        betas.(task_id).mdl_3.kEffort,...
        betas.(task_id).mdl_3.kFatigue,...
        betas.(task_id).mdl_4.kR,...
        betas.(task_id).mdl_4.kP,...
        betas.(task_id).mdl_4.kEffort,...
        betas.(task_id).mdl_4.kFatigue] = deal(NaN(1,NS));

    [choiceNonDef.perMoneyLevel.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl_1.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl_2.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl_3.(task_id),...
        choiceFitNonDef.perMoneyLevel.mdl_4.(task_id)] = deal(NaN(nMoneyLevels,NS));
    [choiceNonDef.perEffortLevel.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl_1.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl_2.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl_3.(task_id),...
        choiceFitNonDef.perEffortLevel.mdl_4.(task_id)] = deal(NaN(nELevels,NS));
    [choiceNonDef.perTrialN.(task_id),...
        choiceFitNonDef.perTrialN.mdl_1.(task_id),...
        choiceFitNonDef.perTrialN.mdl_2.(task_id),...
        choiceFitNonDef.perTrialN.mdl_3.(task_id),...
        choiceFitNonDef.perTrialN.mdl_4.(task_id)] = deal(NaN(n_trialN_bins,NS));

    [choiceNonDef.perNVLevel.mdl_1.(task_id),...
        choiceFitNonDef.perNVLevel.mdl_1.(task_id),...
        NV_bins.mdl_1.(task_id),...
        choiceNonDef.perNVLevel.mdl_2.(task_id),...
        choiceFitNonDef.perNVLevel.mdl_2.(task_id),...
        NV_bins.mdl_2.(task_id),...
        choiceNonDef.perNVLevel.mdl_3.(task_id),...
        choiceFitNonDef.perNVLevel.mdl_3.(task_id),...
        NV_bins.mdl_3.(task_id),...
        choiceNonDef.perNVLevel.mdl_4.(task_id),...
        choiceFitNonDef.perNVLevel.mdl_4.(task_id),...
        NV_bins.mdl_4.(task_id)] = deal(NaN(n_NV_bins, NS));
    actualMoney_values.(task_id) = NaN(nMoneyLevels, NS);
end % physical/mental loop

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    % load individual data
    [betas_tmp, choices_tmp] = logitfit_choices(computerRoot, study_nm, sub_nm,...
        figDispIndiv, dispMoneyOrLevels, n_NV_bins, n_trialN_bins);
    trialN_levels = choices_tmp.trialN_bins.Ep;

    % pool data across subjects
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
            case 2
                task_id = 'Em';
        end
        % extract betas
        % model 1
        betas.(task_id).mdl_1.kMoney(iS) = betas_tmp.(task_id).mdl_1.kMoney;
        betas.(task_id).mdl_1.kEffort(iS) = betas_tmp.(task_id).mdl_1.kEffort;
        % model 2
        betas.(task_id).mdl_2.kR(iS) = betas_tmp.(task_id).mdl_2.kR;
        betas.(task_id).mdl_2.kP(iS) = betas_tmp.(task_id).mdl_2.kP;
        betas.(task_id).mdl_2.kEffort(iS) = betas_tmp.(task_id).mdl_2.kEffort;
        % model 3
        betas.(task_id).mdl_3.kMoney(iS) = betas_tmp.(task_id).mdl_3.kMoney;
        betas.(task_id).mdl_3.kEffort(iS) = betas_tmp.(task_id).mdl_3.kEffort;
        betas.(task_id).mdl_3.kFatigue(iS) = betas_tmp.(task_id).mdl_3.kFatigue;
        % model 4
        betas.(task_id).mdl_4.kR(iS) = betas_tmp.(task_id).mdl_4.kR;
        betas.(task_id).mdl_4.kP(iS) = betas_tmp.(task_id).mdl_4.kP;
        betas.(task_id).mdl_4.kEffort(iS) = betas_tmp.(task_id).mdl_4.kEffort;
        betas.(task_id).mdl_4.kFatigue(iS) = betas_tmp.(task_id).mdl_4.kFatigue;

        % extract choices
        choiceNonDef.perMoneyLevel.(task_id)(:,iS) = choices_tmp.choiceNonDef.perMoneyLevel.(task_id);
        choiceNonDef.perEffortLevel.(task_id)(:,iS) = choices_tmp.choiceNonDef.perEffortLevel.(task_id);
        choiceNonDef.perTrialN.(task_id)(:,iS) = choices_tmp.choiceNonDef.perTrialN.(task_id);
        % extract fit
        for iMdl = 1:nMdl
            mdl_nm = ['mdl_',num2str(iMdl)];
            choiceNonDef.perNVLevel.(mdl_nm).(task_id)(:,iS) = choices_tmp.choiceNonDef.perNVLevel.(mdl_nm).(task_id);
            choiceFitNonDef.perNVLevel.(mdl_nm).(task_id)(:,iS) = choices_tmp.choiceFitNonDef.perNVLevel.(mdl_nm).(task_id);
            NV_bins.(mdl_nm).(task_id)(:,iS) = choices_tmp.NV_bins.(mdl_nm).(task_id);
            choiceFitNonDef.perMoneyLevel.(mdl_nm).(task_id)(:,iS) = choices_tmp.choiceFitNonDef.perMoneyLevel.(mdl_nm).(task_id);
            choiceFitNonDef.perEffortLevel.(mdl_nm).(task_id)(:,iS) = choices_tmp.choiceFitNonDef.perEffortLevel.(mdl_nm).(task_id);
            choiceFitNonDef.perTrialN.(mdl_nm).(task_id)(:,iS) = choices_tmp.choiceFitNonDef.perTrialN.(mdl_nm).(task_id);
        end

        % extract actual money levels (linked to IP measurement)
        actualMoney_values.(task_id)(:, iS) = choices_tmp.actualMoney_values.(task_id);
    end % physical/mental loop
end % subject loop

%% average data depending on metabolite levels
group_names = {['low_',metabolite_nm],['high_',metabolite_nm]};
n_grps = length(group_names);
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
            task_fullName = 'physical';
        case 2
            task_id = 'Em';
            task_fullName = 'mental';
    end

    % loop through 2 groups (low vs high metabolites)
    for iGrp = 1:n_grps
        grp_nm = group_names{iGrp};
        switch grp_nm
            case {['low_',metabolite_nm]}
                subjectsToInclude = low_met_subs;
            case {['high_',metabolite_nm]}
                subjectsToInclude = high_met_subs;
        end
        [m_choiceNonDef.perMoneyLevel.(task_id).(grp_nm),...
            sem_choiceNonDef.perMoneyLevel.(task_id).(grp_nm)] = mean_sem_sd(choiceNonDef.perMoneyLevel.(task_id)(:,subjectsToInclude), 2);
        [m_choiceNonDef.perEffortLevel.(task_id).(grp_nm),...
            sem_choiceNonDef.perEffortLevel.(task_id).(grp_nm)] = mean_sem_sd(choiceNonDef.perEffortLevel.(task_id)(:,subjectsToInclude), 2);
        [m_choiceNonDef.perTrialN.(task_id).(grp_nm),...
            sem_choiceNonDef.perTrialN.(task_id).(grp_nm)] = mean_sem_sd(choiceNonDef.perTrialN.(task_id)(:,subjectsToInclude), 2);
        for iMdl = 1:nMdl
            mdl_nm = ['mdl_',num2str(iMdl)];
            [m_choiceFitNonDef.perMoneyLevel.(mdl_nm).(task_id).(grp_nm),...
                sem_choiceFitNonDef.perMoneyLevel.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(choiceFitNonDef.perMoneyLevel.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            [m_choiceFitNonDef.perEffortLevel.(mdl_nm).(task_id).(grp_nm),...
                sem_choiceFitNonDef.perEffortLevel.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(choiceFitNonDef.perEffortLevel.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            [m_NV_bins.(mdl_nm).(task_id).(grp_nm),...
                sem_NV_bins.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(NV_bins.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            [m_choiceNonDef.perNVLevel.(mdl_nm).(task_id).(grp_nm),...
                sem_choiceNonDef.perNVLevel.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(choiceNonDef.perNVLevel.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            [m_choiceFitNonDef.perNVLevel.(mdl_nm).(task_id).(grp_nm),...
                sem_choiceFitNonDef.perNVLevel.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(choiceFitNonDef.perNVLevel.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            [m_choiceFitNonDef.perTrialN.(mdl_nm).(task_id).(grp_nm),...
                sem_choiceFitNonDef.perTrialN.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(choiceFitNonDef.perTrialN.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            % average betas
            if ismember(iMdl,[1,3])
                [betas.mean.(task_id).(mdl_nm).kMoney.(grp_nm),...
                    betas.sem.(task_id).(mdl_nm).kMoney.(grp_nm),...
                    betas.sd.(task_id).(mdl_nm).kMoney.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).kMoney(:,subjectsToInclude),2);
            elseif ismember(iMdl,[2,4])
                [betas.mean.(task_id).(mdl_nm).kR.(grp_nm),...
                    betas.sem.(task_id).(mdl_nm).kR.(grp_nm),...
                    betas.sd.(task_id).(mdl_nm).kR.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).kR(:,subjectsToInclude),2);
                [betas.mean.(task_id).(mdl_nm).kP.(grp_nm),...
                    betas.sem.(task_id).(mdl_nm).kP.(grp_nm),...
                    betas.sd.(task_id).(mdl_nm).kP.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).kP(:,subjectsToInclude),2);
            end
            [betas.mean.(task_id).(mdl_nm).kEffort.(grp_nm),...
                betas.sem.(task_id).(mdl_nm).kEffort.(grp_nm),...
                betas.sd.(task_id).(mdl_nm).kEffort.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).kEffort(:,subjectsToInclude),2);
            if ismember(iMdl,[3,4])
                [betas.mean.(task_id).(mdl_nm).kFatigue.(grp_nm),...
                    betas.sem.(task_id).(mdl_nm).kFatigue.(grp_nm),...
                    betas.sd.(task_id).(mdl_nm).kFatigue.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).kFatigue(:,subjectsToInclude),2);
            end

            %% test how significant betas are
            if ismember(iMdl,[1,3])
                [~,pvalues.(task_id).(mdl_nm).kMoney.(grp_nm)] = ttest(betas.(task_id).(mdl_nm).kMoney(:,subjectsToInclude));
            elseif ismember(iMdl,[2,4])
                [~,pvalues.(task_id).(mdl_nm).kR.(grp_nm)] = ttest(betas.(task_id).(mdl_nm).kR(:,subjectsToInclude));
                [~,pvalues.(task_id).(mdl_nm).kP.(grp_nm)] = ttest(betas.(task_id).(mdl_nm).kP(:,subjectsToInclude));
            end
            [~,pvalues.(task_id).(mdl_nm).kEffort.(grp_nm)] = ttest(betas.(task_id).(mdl_nm).kEffort(:,subjectsToInclude));
            if ismember(iMdl,[3,4])
                [~,pvalues.(task_id).(mdl_nm).kFatigue.(grp_nm)] = ttest(betas.(task_id).(mdl_nm).kFatigue(:,subjectsToInclude));
            end

            %% average money levels
            [m_actualMoney_values.(task_id),...
                sem_actualMoney_values.(task_id)] = mean_sem_sd(actualMoney_values.(task_id), 2);
        end % model loop

        % compare the betas of high group
        for iMdl = 1:nMdl
            mdl_nm = ['mdl_',num2str(iMdl)];
            if ismember(iMdl,[1,3])
                [~,pvalues.(task_id).(mdl_nm).kMoney.(grp_nm)] = ttest2(betas.(task_id).(mdl_nm).kMoney(:,low_met_subs),...
                    betas.(task_id).(mdl_nm).kMoney(:,high_met_subs));
            elseif ismember(iMdl,[2,4])
                [~,pvalues.(task_id).(mdl_nm).kR.(grp_nm)] = ttest2(betas.(task_id).(mdl_nm).kR(:,low_met_subs),...
                    betas.(task_id).(mdl_nm).kR(:,high_met_subs));
                [~,pvalues.(task_id).(mdl_nm).kP.(grp_nm)] = ttest2(betas.(task_id).(mdl_nm).kP(:,low_met_subs),...
                    betas.(task_id).(mdl_nm).kP(:,high_met_subs));
            end
            [~,pvalues.(task_id).(mdl_nm).kEffort.(grp_nm)] = ttest2(betas.(task_id).(mdl_nm).kEffort(:,low_met_subs),...
                    betas.(task_id).(mdl_nm).kEffort(:,high_met_subs));
            if ismember(iMdl,[3,4])
                [~,pvalues.(task_id).(mdl_nm).kFatigue.(grp_nm)] = ttest2(betas.(task_id).(mdl_nm).kFatigue(:,low_met_subs),...
                    betas.(task_id).(mdl_nm).kFatigue(:,high_met_subs));
            end
        end
    end % loop through low/high metabolite levels

    %% display average data
    if figDispGroup == 1
        pSize = 30;
        lWidth = 3;
        %         lWidth_borders = 1;
        %% loop through models
        for iMdl = 4:nMdl
            mdl_nm = ['mdl_',num2str(iMdl)];

            %% display choice = f(net value)
            fig;
            for iGrp = 1:n_grps
                grp_nm = group_names{iGrp};
                %         pointMdl1 = errorbar(m_NV_bins.mdl_1.(task_id).(grp_nm),...
                %             m_choiceNonDef.perNVLevel.mdl_1.(task_id).(grp_nm),...
                %             sem_choiceNonDef.perNVLevel.mdl_1.(task_id).(grp_nm),...
                %             sem_choiceNonDef.perNVLevel.mdl_1.(task_id).(grp_nm),...
                %             sem_NV_bins.mdl_1.(task_id).(grp_nm),...
                %             sem_NV_bins.mdl_1.(task_id).(grp_nm));
                pointMdl = errorbar(m_NV_bins.(mdl_nm).(task_id).(grp_nm),...
                    m_choiceNonDef.perNVLevel.(mdl_nm).(task_id).(grp_nm).*100,...
                    sem_choiceNonDef.perNVLevel.(mdl_nm).(task_id).(grp_nm).*100);
                pointMdl.Color = [0 0 0];
                switch iGrp
                    case 1
                        pointMdl.Marker = 'o';
                    case 2
                        pointMdl.Marker = '+';
                end
                pointMdl.LineStyle = 'none';
                pointMdl.LineWidth = lWidth;
                hold on;
                %             line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                %             line(xlim(),[100 100],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                %             line(xlim(),[50 50],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                %             line([0 0],[0 100],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                lHdlMdl = plot(m_NV_bins.(mdl_nm).(task_id).(grp_nm),...
                    m_choiceFitNonDef.perNVLevel.(mdl_nm).(task_id).(grp_nm).*100);
                switch iGrp
                    case 1
                        lHdlMdl.LineStyle = '--';
                    case 2
                        lHdlMdl.LineStyle = ':';
                end
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = [143 0 0]./255;
            end % low/high metabolites loop
            ylim([0 100]);
            %             xlabel(['Net value ',task_fullName,' - model ',num2str(iMdl)]);
            xlabel('Net value non-default - default option');
            ylabel('Choice non-default option (%)');
            legend_size(pSize);
            %         saveas(resultFolder)

            %% choice non-default = f(money levels)
            fig;
            for iGrp = 1:n_grps
                grp_nm = group_names{iGrp};
                switch dispMoneyOrLevels
                    case 'money' % show also error bar in the X dimension
                        money_or_levels = m_actualMoney_values.(task_id);
                        pointMdl = errorbar(m_actualMoney_values.(task_id),...
                            m_choiceNonDef.perMoneyLevel.(task_id).(grp_nm).*100,...
                            m_choiceNonDef.perMoneyLevel.(task_id).(grp_nm)-sem_choiceNonDef.perMoneyLevel.(task_id).(grp_nm).*100,...
                            m_choiceNonDef.perMoneyLevel.(task_id).(grp_nm)+sem_choiceNonDef.perMoneyLevel.(task_id).(grp_nm).*100,...
                            m_actualMoney_values.(task_id).(grp_nm)-sem_actualMoney_values.(task_id).(grp_nm),...
                            m_actualMoney_values.(task_id).(grp_nm)+sem_actualMoney_values.(task_id).(grp_nm));
                    case 'levels'
                        money_or_levels = money_levels;
                        pointMdl = errorbar(money_levels,...
                            m_choiceNonDef.perMoneyLevel.(task_id).(grp_nm).*100,...
                            sem_choiceNonDef.perMoneyLevel.(task_id).(grp_nm).*100);
                end
                pointMdl.Color = [0 0 0];
                switch iGrp
                    case 1
                        pointMdl.Marker = 'o';
                    case 2
                        pointMdl.Marker = '+';
                end
                pointMdl.LineStyle = 'none';
                pointMdl.LineWidth = lWidth;
                hold on;
                %             line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                %             line(xlim(),[100 100],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                lHdlMdl = plot(money_or_levels,...
                    m_choiceFitNonDef.perMoneyLevel.(mdl_nm).(task_id).(grp_nm).*100);
                switch iGrp
                    case 1
                        lHdlMdl.LineStyle = '--';
                    case 2
                        lHdlMdl.LineStyle = ':';
                end
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = [143 0 0]./255;
            end % low/high metabolites loop
            %             ylim([-0.2 1.2]);
            ylim([0 100]);
            switch dispMoneyOrLevels
                case 'money'
                    xlabel([task_fullName,' Money (€) - model ',num2str(iMdl)]);
                case 'levels'
                    xlabel([task_fullName,' Money level - model ',num2str(iMdl)]);
            end
            ylabel('Choice non-default option (%)');
            legend_size(pSize);
            %         saveas(resultFolder)

            %% choice non-default = f(effort levels)
            fig;
            for iGrp = 1:n_grps
                grp_nm = group_names{iGrp};
                pointMdl = errorbar(E_levels,...
                    m_choiceNonDef.perEffortLevel.(task_id).(grp_nm).*100,...
                    sem_choiceNonDef.perEffortLevel.(task_id).(grp_nm).*100);
                pointMdl.Color = [0 0 0];
                switch iGrp
                    case 1
                        pointMdl.Marker = 'o';
                    case 2
                        pointMdl.Marker = '+';
                end
                pointMdl.LineStyle = 'none';
                pointMdl.LineWidth = lWidth;
                hold on;
                %             line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                %             line(xlim(),[100 100],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                lHdlMdl = plot(E_levels,...
                    m_choiceFitNonDef.perEffortLevel.(mdl_nm).(task_id).(grp_nm).*100);
                switch iGrp
                    case 1
                        lHdlMdl.LineStyle = '--';
                    case 2
                        lHdlMdl.LineStyle = ':';
                end
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = [143 0 0]./255;
            end % low/high metabolites loop
            %             ylim([-0.2 1.2]);
            ylim([0 100]);
            xlabel([task_fullName,' effort level - model ',num2str(iMdl)]);
            ylabel('Choice non-default option (%)');
            legend_size(pSize);
            %         saveas(resultFolder)

            %% choice non-default = f(trial number)
            fig;
            for iGrp = 1:n_grps
                grp_nm = group_names{iGrp};
                pointMdl = errorbar(trialN_levels,...
                    m_choiceNonDef.perTrialN.(task_id).(grp_nm).*100,...
                    sem_choiceNonDef.perTrialN.(task_id).(grp_nm).*100);
                pointMdl.Color = [0 0 0];
                switch iGrp
                    case 1
                        pointMdl.Marker = 'o';
                    case 2
                        pointMdl.Marker = '+';
                end
                pointMdl.LineStyle = 'none';
                pointMdl.LineWidth = lWidth;
                hold on;
                %             line(xlim(),[0 0],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                %             line(xlim(),[100 100],'LineWidth',lWidth_borders,'Color',[0 0 0]);
                lHdlMdl = plot(trialN_levels,...
                    m_choiceFitNonDef.perTrialN.(mdl_nm).(task_id).(grp_nm).*100);
                switch iGrp
                    case 1
                        lHdlMdl.LineStyle = '--';
                    case 2
                        lHdlMdl.LineStyle = ':';
                end
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = [143 0 0]./255;
            end % low/high metabolites loop
            %             ylim([-0.2 1.2]);
            ylim([0 100]);
            xlabel([task_fullName,' trial number - model ',num2str(iMdl)]);
            ylabel('Choice non-default option (%)');
            legend_size(pSize);
            %         saveas(resultFolder)
        end % model loop
    end % display figure

end % physical/mental loop


end % function