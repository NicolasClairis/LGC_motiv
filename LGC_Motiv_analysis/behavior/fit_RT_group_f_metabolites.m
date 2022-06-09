function[betas, pval] = fit_RT_group_f_metabolites(dispMoneyOrLevels, RT_type)
% [betas, pval] = fit_RT_group(dispMoneyOrLevels, RT_type)
%
% INPUTS
% dispMoneyOrLevels: display actual money ('money') or reward levels
% ('levels')
%
% RT_type:
% 'raw': raw reaction times (in seconds)
% 'zscored': zscore reaction times per run
% 'log': take log(RT)
%
% OUTPUT
% betas: structure with betas
%
% pval: structure with p.values

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
nTrialsPerRun = 54;

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
condition = subject_condition();
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
    [betas.(task_id).mdl_0.b0,...
        betas.(task_id).mdl_0.bR,...
        betas.(task_id).mdl_0.bP,...
        betas.(task_id).mdl_0.bRP,...
        betas.(task_id).mdl_0.bE,...
        betas.(task_id).mdl_0.bConf,...
        betas.(task_id).mdl_1.b0,...
        betas.(task_id).mdl_1.bNV,...
        betas.(task_id).mdl_1.bConf,...
        betas.(task_id).mdl_2.b0,...
        betas.(task_id).mdl_2.bNV,...
        betas.(task_id).mdl_2.bConf,...
        betas.(task_id).mdl_3.b0,...
        betas.(task_id).mdl_3.bNV,...
        betas.(task_id).mdl_3.bConf,...
        betas.(task_id).mdl_4.b0,...
        betas.(task_id).mdl_4.bNV,...
        betas.(task_id).mdl_4.bConf] = deal(NaN(1,NS));

    [RT.perMoneyLevel.(task_id),...
        RTfit.perMoneyLevel.mdl_0.(task_id),...
        RTfit.perMoneyLevel.mdl_1.(task_id),...
        RTfit.perMoneyLevel.mdl_2.(task_id),...
        RTfit.perMoneyLevel.mdl_3.(task_id),...
        RTfit.perMoneyLevel.mdl_4.(task_id)] = deal(NaN(nMoneyLevels,NS));
    [RT.perEffortLevel.(task_id),...
        RTfit.perEffortLevel.mdl_0.(task_id),...
        RTfit.perEffortLevel.mdl_1.(task_id),...
        RTfit.perEffortLevel.mdl_2.(task_id),...
        RTfit.perEffortLevel.mdl_3.(task_id),...
        RTfit.perEffortLevel.mdl_4.(task_id)] = deal(NaN(nELevels,NS));
    [RT.perTrialN.(task_id),...
        RTfit.perTrialN.mdl_0.(task_id),...
        RTfit.perTrialN.mdl_1.(task_id),...
        RTfit.perTrialN.mdl_2.(task_id),...
        RTfit.perTrialN.mdl_3.(task_id),...
        RTfit.perTrialN.mdl_4.(task_id)] = deal(NaN(n_trialN_bins,NS));
    [RT.perConfidenceLevel.(task_id),...
        RTfit.perConfidenceLevel.mdl_0.(task_id),...
        RTfit.perConfidenceLevel.mdl_1.(task_id),...
        RTfit.perConfidenceLevel.mdl_2.(task_id),...
        RTfit.perConfidenceLevel.mdl_3.(task_id),...
        RTfit.perConfidenceLevel.mdl_4.(task_id)] = deal(NaN(2,NS));
    [RT.RP.(task_id),...
        RTfit.RP.mdl_0.(task_id),...
        RTfit.RP.mdl_1.(task_id),...
        RTfit.RP.mdl_2.(task_id),...
        RTfit.RP.mdl_3.(task_id),...
        RTfit.RP.mdl_4.(task_id)] = deal(NaN(2,NS));

    [RT.perNVlevel.mdl_1.(task_id),...
        RTfit.perNVlevel.mdl_1.(task_id),...
        NV_bins.mdl_1.(task_id),...
        RT.perNVlevel.mdl_2.(task_id),...
        RTfit.perNVlevel.mdl_2.(task_id),...
        NV_bins.mdl_2.(task_id),...
        RT.perNVlevel.mdl_3.(task_id),...
        RTfit.perNVlevel.mdl_3.(task_id),...
        NV_bins.mdl_3.(task_id),...
        RT.perNVlevel.mdl_4.(task_id),...
        RTfit.perNVlevel.mdl_4.(task_id),...
        NV_bins.mdl_4.(task_id)] = deal(NaN(n_NV_bins, NS));
    actualMoney_values.(task_id) = NaN(nMoneyLevels, NS);
end % physical/mental loop

%% extract individual data
for iS = 1:NS
    sub_nm = subject_id{iS};

    % get the data
    [betas_RT_perSub, RTstruct_perSub, pval_perSub] = fit_RT(computerRoot, study_nm, sub_nm,...
        figDispIndiv, dispMoneyOrLevels, n_NV_bins, n_trialN_bins, RT_type);
    trialN_levels = RTstruct_perSub.trialN_bins.Ep;

    % pool data across subjects
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
            case 2
                task_id = 'Em';
        end

        % extract bins RT
        RT.perMoneyLevel.(task_id)(:, iS)       = RTstruct_perSub.RT_bins.perMoneyLevel.(task_id);
        RT.perEffortLevel.(task_id)(:, iS)      = RTstruct_perSub.RT_bins.perEffortLevel.(task_id);
        RT.perTrialN.(task_id)(:, iS)           = RTstruct_perSub.RT_bins.perTrialN.(task_id);
        RT.RP.(task_id)(:, iS)                  = RTstruct_perSub.RT_bins.RP.(task_id);
        RT.perConfidenceLevel.(task_id)(:, iS)  = RTstruct_perSub.RT_bins.perConfidenceLevel.(task_id);

        % extract betas and bins
        for iMdl = 0:nMdl
            mdl_nm = ['mdl_',num2str(iMdl)];

            % extract betas
            betas.(task_id).(mdl_nm).b0(iS) = betas_RT_perSub.(task_id).(mdl_nm).b0;
            betas.(task_id).(mdl_nm).bConf(iS) = betas_RT_perSub.(task_id).(mdl_nm).bConf;
            switch mdl_nm
                case 'mdl_0'
                    betas.(task_id).(mdl_nm).bR(iS)     = betas_RT_perSub.(task_id).(mdl_nm).bR;
                    betas.(task_id).(mdl_nm).bP(iS)     = betas_RT_perSub.(task_id).(mdl_nm).bP;
                    betas.(task_id).(mdl_nm).bRP(iS)    = betas_RT_perSub.(task_id).(mdl_nm).bRP;
                    betas.(task_id).(mdl_nm).bE(iS)     = betas_RT_perSub.(task_id).(mdl_nm).bE;
                case {'mdl_1','mdl_2','mdl_3','mdl_4'}
                    betas.(task_id).(mdl_nm).bNV(iS) = betas_RT_perSub.(task_id).(mdl_nm).bNV;
            end

            % extract bins for fit
            RTfit.perMoneyLevel.(mdl_nm).(task_id)(:, iS)       = RTstruct_perSub.RTfit_bins.perMoneyLevel.(mdl_nm).(task_id);
            RTfit.perEffortLevel.(mdl_nm).(task_id)(:, iS)      = RTstruct_perSub.RTfit_bins.perEffortLevel.(mdl_nm).(task_id);
            RTfit.perTrialN.(mdl_nm).(task_id)(:, iS)           = RTstruct_perSub.RTfit_bins.perTrialN.(mdl_nm).(task_id);
            RTfit.RP.(mdl_nm).(task_id)(:, iS)                  = RTstruct_perSub.RTfit_bins.RP.(mdl_nm).(task_id);
            RTfit.perConfidenceLevel.(mdl_nm).(task_id)(:, iS)  = RTstruct_perSub.RTfit_bins.perConfidenceLevel.(mdl_nm).(task_id);
            if iMdl > 0
                RT.perNVlevel.(mdl_nm).(task_id)(:, iS)     = RTstruct_perSub.RT_bins.perNVlevel.(mdl_nm).(task_id);
                RTfit.perNVlevel.(mdl_nm).(task_id)(:, iS)  = RTstruct_perSub.RTfit_bins.perNVlevel.(mdl_nm).(task_id);
                NV_bins.(mdl_nm).(task_id)(:,iS)            = RTstruct_perSub.NV_bins.perNVlevel.(mdl_nm).(task_id);
            end
        end % model loop

    end % physical/mental loop
end % subject loop

%% average data depending on metabolite levels
group_names = {['low_',metabolite_nm],['high_',metabolite_nm]};
n_grps = length(group_names);
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
        case 2
            task_id = 'Em';
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
        % average RT
        [m_RT.perMoneyLevel.(task_id).(grp_nm),...
            sem_RT.perMoneyLevel.(task_id).(grp_nm)] = mean_sem_sd(RT.perMoneyLevel.(task_id)(:,subjectsToInclude), 2);
        [m_RT.perEffortLevel.(task_id).(grp_nm),...
            sem_RT.perEffortLevel.(task_id).(grp_nm)] = mean_sem_sd(RT.perEffortLevel.(task_id)(:,subjectsToInclude), 2);
        [m_RT.perTrialN.(task_id).(grp_nm),...
            sem_RT.perTrialN.(task_id).(grp_nm)] = mean_sem_sd(RT.perTrialN.(task_id)(:,subjectsToInclude), 2);
        [m_RT.RP.(task_id).(grp_nm),...
            sem_RT.RP.(task_id).(grp_nm)] = mean_sem_sd(RT.RP.(task_id)(:,subjectsToInclude), 2);
        [m_RT.perConfidenceLevel.(task_id).(grp_nm),...
            sem_RT.perConfidenceLevel.(task_id).(grp_nm)] = mean_sem_sd(RT.perConfidenceLevel.(task_id)(:,subjectsToInclude), 2);

        for iMdl = 0:nMdl
            mdl_nm = ['mdl_',num2str(iMdl)];

            if iMdl > 0
                [m_RT.perNVlevel.(mdl_nm).(task_id).(grp_nm),...
                    sem_RT.perNVlevel.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(RT.perNVlevel.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
                [m_RTfit.perNVlevel.(mdl_nm).(task_id).(grp_nm),...
                    sem_RTfit.perNVlevel.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(RTfit.perNVlevel.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
                m_NV_bins.(mdl_nm).(task_id).(grp_nm) = mean_sem_sd(NV_bins.(mdl_nm).(task_id)(:,subjectsToInclude),2);
            end

            % extract betas
            % beta 0
            [betas.mean.(task_id).(mdl_nm).b0.(grp_nm),...
                betas.sem.(task_id).(mdl_nm).b0.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).b0(:,subjectsToInclude),2);
            [~,pval_tmp] = ttest(betas.(task_id).(mdl_nm).b0(:,subjectsToInclude));
            pval.(task_id).(mdl_nm).b0.(grp_nm) = pval_tmp;
            % beta confidence
            [betas.mean.(task_id).(mdl_nm).bConf.(grp_nm),...
                betas.sem.(task_id).(mdl_nm).bConf.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).bConf(:,subjectsToInclude),2);
            [~,pval_tmp] = ttest(betas.(task_id).(mdl_nm).bConf(:,subjectsToInclude));
            pval.(task_id).(mdl_nm).bConf.(grp_nm) = pval_tmp;
            switch mdl_nm
                case 'mdl_0'
                    % beta net value
                    [betas.mean.(task_id).(mdl_nm).bR.(grp_nm),...
                        betas.sem.(task_id).(mdl_nm).bR.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).bR(:,subjectsToInclude),2);
                    [~,pval_tmp] = ttest(betas.(task_id).(mdl_nm).bR(:,subjectsToInclude));
                    pval.(task_id).(mdl_nm).bR.(grp_nm) = pval_tmp;
                    % beta punishment
                    [betas.mean.(task_id).(mdl_nm).bP.(grp_nm),...
                        betas.sem.(task_id).(mdl_nm).bP.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).bP(:,subjectsToInclude),2);
                    [~,pval_tmp] = ttest(betas.(task_id).(mdl_nm).bP(:,subjectsToInclude));
                    pval.(task_id).(mdl_nm).bP.(grp_nm) = pval_tmp;
                    % beta R vs P
                    [betas.mean.(task_id).(mdl_nm).bRP.(grp_nm),...
                        betas.sem.(task_id).(mdl_nm).bRP.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).bRP(:,subjectsToInclude),2);
                    [~,pval_tmp] = ttest(betas.(task_id).(mdl_nm).bRP(:,subjectsToInclude));
                    pval.(task_id).(mdl_nm).bRP.(grp_nm) = pval_tmp;
                    % beta effort
                    [betas.mean.(task_id).(mdl_nm).bE.(grp_nm),...
                        betas.sem.(task_id).(mdl_nm).bE.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).bE(:,subjectsToInclude),2);
                    [~,pval_tmp] = ttest(betas.(task_id).(mdl_nm).bE(:,subjectsToInclude));
                    pval.(task_id).(mdl_nm).bE.(grp_nm) = pval_tmp;
                case {'mdl_1','mdl_2','mdl_3','mdl_4'}
                    % beta net value
                    [betas.mean.(task_id).(mdl_nm).bNV.(grp_nm),...
                        betas.sem.(task_id).(mdl_nm).bNV.(grp_nm)] = mean_sem_sd(betas.(task_id).(mdl_nm).bNV(:,subjectsToInclude),2);
                    [~,pval_tmp] = ttest(betas.(task_id).(mdl_nm).bNV(:,subjectsToInclude));
                    pval.(task_id).(mdl_nm).bNV.(grp_nm) = pval_tmp;
            end

            % average RT fit
            [m_RTfit.perMoneyLevel.(mdl_nm).(task_id).(grp_nm),...
                sem_RTfit.perMoneyLevel.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(RTfit.perMoneyLevel.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            [m_RTfit.perEffortLevel.(mdl_nm).(task_id).(grp_nm),...
                sem_RTfit.perEffortLevel.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(RTfit.perEffortLevel.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            [m_RTfit.perTrialN.(mdl_nm).(task_id).(grp_nm),...
                sem_RTfit.perTrialN.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(RTfit.perTrialN.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            [m_RTfit.RP.(mdl_nm).(task_id).(grp_nm),...
                sem_RTfit.RP.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(RTfit.RP.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
            [m_RTfit.perConfidenceLevel.(mdl_nm).(task_id).(grp_nm),...
                sem_RTfit.perConfidenceLevel.(mdl_nm).(task_id).(grp_nm)] = mean_sem_sd(RTfit.perConfidenceLevel.(mdl_nm).(task_id)(:,subjectsToInclude), 2);
        end % model loop

    end % loop through low/high metabolite levels
end % physical/mental loop

%% display the data
if figDispGroup == 1
    pSize = 30;
    lWidth = 3;
    lWidth_borders = 1;

    %% loop through models (for the fit)
    for iMdl = 0:nMdl
        mdl_nm = ['mdl_',num2str(iMdl)];

        %% loop through tasks
        for iPM = 1:2
            switch iPM
                case 1
                    task_id = 'Ep';
                    task_fullName = 'physical';
                    fitColor = [0 153/255 1];
                case 2
                    task_id = 'Em';
                    task_fullName = 'mental';
                    fitColor = [0 1 0];
            end

            %% RT = f(monetary amounts)
            switch dispMoneyOrLevels
                case 'money'
                    money_or_levels = actualMoney_values.(task_id);
                case 'levels'
                    money_or_levels = money_levels;
            end
            fig;
            for iGrp = 1:n_grps
                grp_nm = group_names{iGrp};
                pointMdl = errorbar(money_or_levels,...
                    m_RT.perMoneyLevel.(task_id).(grp_nm),...
                    sem_RT.perMoneyLevel.(task_id).(grp_nm));
                pointMdl.Color = [0 0 0];
                pointMdl.Marker = 'o';
                pointMdl.MarkerSize = 10;
                switch iGrp
                    case 1
                        pointMdl.MarkerFaceColor = [0 0 0];
                    case 2
                        pointMdl.MarkerFaceColor = [1 1 1];
                end
                pointMdl.LineStyle = 'none';
                pointMdl.LineWidth = lWidth;
                hold on;
                lHdlMdl = plot(money_or_levels,...
                    m_RTfit.perMoneyLevel.(['mdl_',num2str(iMdl)]).(task_id).(grp_nm));
                switch iGrp
                    case 1
                        lHdlMdl.LineStyle = '--';
                    case 2
                        lHdlMdl.LineStyle = ':';
                end
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = fitColor;
            end % low/high metabolites loop
            switch dispMoneyOrLevels
                case 'money'
                    xlabel([task_fullName,' Money (€) - model ',num2str(iMdl)]);
                case 'levels'
                    xlabel([task_fullName,' Money level - model ',num2str(iMdl)]);
            end
            y_legend(RT_type);
            legend_size(pSize);

            %% RT = f(R/P)
            fig;
            for iGrp = 1:n_grps
                grp_nm = group_names{iGrp};
                pointMdl = errorbar(1:2,...
                    m_RT.RP.(task_id).(grp_nm),...
                    sem_RT.RP.(task_id).(grp_nm));
                pointMdl.Color = [0 0 0];
                pointMdl.Marker = 'o';
                pointMdl.MarkerSize = 10;
                switch iGrp
                    case 1
                        pointMdl.MarkerFaceColor = [0 0 0];
                    case 2
                        pointMdl.MarkerFaceColor = [1 1 1];
                end
                pointMdl.LineStyle = 'none';
                pointMdl.LineWidth = lWidth;
                hold on;
                lHdlMdl = plot(1:2,...
                    m_RTfit.RP.(['mdl_',num2str(iMdl)]).(task_id).(grp_nm));
                switch iGrp
                    case 1
                        lHdlMdl.LineStyle = '--';
                    case 2
                        lHdlMdl.LineStyle = ':';
                end
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = fitColor;
            end % low/high metabolites loop
            xticks(1:2);
            xticklabels({'R','P'});
            xlim([0.5 2.5]);
            xlabel(['Model ',num2str(iMdl)]);
            y_legend(RT_type);
            legend_size(pSize);

            %% RT = f(E levels)
            fig;
            for iGrp = 1:n_grps
                grp_nm = group_names{iGrp};
                pointMdl = errorbar(E_levels,...
                    m_RT.perEffortLevel.(task_id).(grp_nm),...
                    sem_RT.perEffortLevel.(task_id).(grp_nm));
                pointMdl.Color = [0 0 0];
                pointMdl.Marker = 'o';
                pointMdl.MarkerSize = 10;
                switch iGrp
                    case 1
                        pointMdl.MarkerFaceColor = [0 0 0];
                    case 2
                        pointMdl.MarkerFaceColor = [1 1 1];
                end
                pointMdl.LineStyle = 'none';
                pointMdl.LineWidth = lWidth;
                hold on;
                lHdlMdl = plot(E_levels,...
                    m_RTfit.perEffortLevel.(['mdl_',num2str(iMdl)]).(task_id).(grp_nm));
                switch iGrp
                    case 1
                        lHdlMdl.LineStyle = '--';
                    case 2
                        lHdlMdl.LineStyle = ':';
                end
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = fitColor;
            end % low/high metabolites loop
            xlabel([task_fullName,' effort level - model ',num2str(iMdl)']);
            y_legend(RT_type);
            legend_size(pSize);

            %% RT = f(NV)
            if iMdl > 0
                fig;
                for iGrp = 1:n_grps
                    grp_nm = group_names{iGrp};
                    pointMdl = errorbar(m_NV_bins.(mdl_nm).(task_id).(grp_nm),...
                        m_RT.perNVlevel.(mdl_nm).(task_id).(grp_nm),...
                        sem_RT.perNVlevel.(mdl_nm).(task_id).(grp_nm));
                    pointMdl.Color = [0 0 0];
                    pointMdl.Marker = 'o';
                    pointMdl.MarkerSize = 10;
                    switch iGrp
                        case 1
                            pointMdl.MarkerFaceColor = [0 0 0];
                        case 2
                            pointMdl.MarkerFaceColor = [1 1 1];
                    end
                    pointMdl.LineStyle = 'none';
                    pointMdl.LineWidth = lWidth;
                    hold on;
                    lHdlMdl = plot(m_NV_bins.(mdl_nm).(task_id).(grp_nm),...
                        m_RTfit.perNVlevel.(['mdl_',num2str(iMdl)]).(task_id).(grp_nm));
                    switch iGrp
                        case 1
                            lHdlMdl.LineStyle = '--';
                        case 2
                            lHdlMdl.LineStyle = ':';
                    end
                    lHdlMdl.LineWidth = lWidth;
                    lHdlMdl.Color = fitColor;
                end % low/high metabolites loop
                xlabel([task_fullName,' net value - model ',num2str(iMdl)']);
                y_legend(RT_type);
                legend_size(pSize);
            end

            %% RT = f(confidence)
            fig;
            for iGrp = 1:n_grps
                grp_nm = group_names{iGrp};
                pointMdl = errorbar(1:2,...
                    m_RT.perConfidenceLevel.(task_id).(grp_nm),...
                    sem_RT.perConfidenceLevel.(task_id).(grp_nm));
                pointMdl.Color = [0 0 0];
                pointMdl.Marker = 'o';
                pointMdl.MarkerSize = 10;
                switch iGrp
                    case 1
                        pointMdl.MarkerFaceColor = [0 0 0];
                    case 2
                        pointMdl.MarkerFaceColor = [1 1 1];
                end
                pointMdl.LineStyle = 'none';
                pointMdl.LineWidth = lWidth;
                hold on;
                lHdlMdl = plot(1:2,...
                    m_RTfit.perConfidenceLevel.(['mdl_',num2str(iMdl)]).(task_id).(grp_nm));
                switch iGrp
                    case 1
                        lHdlMdl.LineStyle = '--';
                    case 2
                        lHdlMdl.LineStyle = ':';
                end
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = fitColor;
            end % low/high metabolites loop
            xticks(1:2);
            xticklabels({'low','high'});
            xlabel([task_fullName,' Confidence - model ',num2str(iMdl)']);
            xlim([0.5 2.5]);
            y_legend(RT_type);
            legend_size(pSize);

            %% RT = f(trial number)
            fig;
            for iGrp = 1:n_grps
                grp_nm = group_names{iGrp};
                pointMdl = errorbar(trialN_levels,...
                    m_RT.perTrialN.(task_id).(grp_nm),...
                    sem_RT.perTrialN.(task_id).(grp_nm));
                pointMdl.Color = [0 0 0];
                pointMdl.Marker = 'o';
                pointMdl.MarkerSize = 10;
                switch iGrp
                    case 1
                        pointMdl.MarkerFaceColor = [0 0 0];
                    case 2
                        pointMdl.MarkerFaceColor = [1 1 1];
                end
                pointMdl.LineStyle = 'none';
                pointMdl.LineWidth = lWidth;
                hold on;
                lHdlMdl = plot(trialN_levels,...
                    m_RTfit.perTrialN.(['mdl_',num2str(iMdl)]).(task_id).(grp_nm));
                switch iGrp
                    case 1
                        lHdlMdl.LineStyle = '--';
                    case 2
                        lHdlMdl.LineStyle = ':';
                end
                lHdlMdl.LineWidth = lWidth;
                lHdlMdl.Color = fitColor;
            end % low/high metabolites loop
            xlabel([task_fullName,' trial number - model ',num2str(iMdl)']);
            xlim([0 nTrialsPerRun]);
            y_legend(RT_type);
            legend_size(pSize);
        end % physical/mental task
    end % model loop
end % display figures

end % function


function[] = y_legend(RT_type)
switch RT_type
    case 'raw'
        ylabel('RT (s)');
    case 'zscored'
        ylabel('z(RT) (a.u.)');
    case 'log'
        ylabel('log(RT) (a.u.)');
end
end