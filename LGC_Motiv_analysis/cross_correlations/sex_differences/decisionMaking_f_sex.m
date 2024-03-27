function[choices, parameters, RT, deltaIP] = decisionMaking_f_sex(fig_disp, rmv_outliers_yn)
% [choices, parameters, RT, deltaIP] = decisionMaking_f_sex(fig_disp, rmv_outliers_yn)
% decisionMaking_f_sex will compare all the variables related to the
% decision-making process of the experiment between males and females.
% Variables include the total proportion of choices (HE), the proportion of
% physical (HPE) or mental (HME) choices, and all the parameters derived
% from our computational model (kR, kP, kEp, kEm, kFp, kLm, bias).
%
% INPUTS
% fig_disp: display figures? (1) yes (0) no
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% choices: structure comparing choice-related variables between male and female
%
% parameters: structure comparing behavioral parameters between male and female
%
% RT: structure comparing RT between male and female
%
% deltaIP: difference for indifference point between male and female

%% define inputs by default
% display figure by default
if ~exist('fig_disp','var') || isempty(fig_disp) || ~ismember(fig_disp,[0,1])
    fig_disp = 1;
end
% remove outliers by default
if ~exist('rmv_outliers_yn','var') || isempty(rmv_outliers_yn) || ~ismember(rmv_outliers_yn,[0,1])
    rmv_outliers_yn = 1;
end

%% subject selection
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

%% load variables of interest
% load proportion of high-effort choices
fig_disp0 = 0;
[choice_hE_males] = choice_hE_proportion(study_nm, condition, male_CIDS, fig_disp0);
[choice_hE_females] = choice_hE_proportion(study_nm, condition, female_CIDS, fig_disp0);

% load behavioral parameters
prm_names = {'kR','kP','kEp','kEm',...
    'kBiasM','kFp','kLm'};
nPrm = length(prm_names);
[prm_males, mdlType, mdlN] = prm_extraction(study_nm, male_CIDS);
[prm_females, ~, ~] = prm_extraction(study_nm, female_CIDS);

% load RT
[RT_summary_males] = RT_range(male_CIDS, condition, fig_disp0);
[RT_summary_females] = RT_range(female_CIDS, condition, fig_disp0);

% load indifference point
[deltaIP_Ep_males,...
    deltaIP_Em_males] = extract_IP(male_CIDS, male_NS);
[deltaIP_Ep_females,...
    deltaIP_Em_females] = extract_IP(female_CIDS, female_NS);

%% remove outliers
if rmv_outliers_yn == 1
    % filter choices
    % for males
    [~,~,choice_hE_males.EpEm] = rmv_outliers_3sd(choice_hE_males.EpEm);
    [~,~,choice_hE_males.Ep] = rmv_outliers_3sd(choice_hE_males.Ep);
    [~,~,choice_hE_males.Em] = rmv_outliers_3sd(choice_hE_males.Em);
    % same for females
    [~,~,choice_hE_females.EpEm] = rmv_outliers_3sd(choice_hE_females.EpEm);
    [~,~,choice_hE_females.Ep] = rmv_outliers_3sd(choice_hE_females.Ep);
    [~,~,choice_hE_females.Em] = rmv_outliers_3sd(choice_hE_females.Em);
    
    % filter parameters
    for iPrm = 1:nPrm
        prm_nm = prm_names{iPrm};
        [~,~,prm_males.(prm_nm)] = rmv_outliers_3sd(prm_males.(prm_nm));
        [~,~,prm_females.(prm_nm)] = rmv_outliers_3sd(prm_females.(prm_nm));
    end
    
    % filter RT
    % males
    [~,~,RT_summary_males.EpEm.mean_RT.allSubs.choice] = rmv_outliers_3sd(RT_summary_males.EpEm.mean_RT.allSubs.choice);
    [~,~,RT_summary_males.Ep.mean_RT.allSubs.choice] = rmv_outliers_3sd(RT_summary_males.Ep.mean_RT.allSubs.choice);
    [~,~,RT_summary_males.Em.mean_RT.allSubs.choice] = rmv_outliers_3sd(RT_summary_males.Em.mean_RT.allSubs.choice);
    % females
    [~,~,RT_summary_females.EpEm.mean_RT.allSubs.choice] = rmv_outliers_3sd(RT_summary_females.EpEm.mean_RT.allSubs.choice);
    [~,~,RT_summary_females.Ep.mean_RT.allSubs.choice] = rmv_outliers_3sd(RT_summary_females.Ep.mean_RT.allSubs.choice);
    [~,~,RT_summary_females.Em.mean_RT.allSubs.choice] = rmv_outliers_3sd(RT_summary_females.Em.mean_RT.allSubs.choice);
    
    % indifference point
    % males
    [~,~,deltaIP_Ep_males] = rmv_outliers_3sd(deltaIP_Ep_males);
    [~,~,deltaIP_Em_males] = rmv_outliers_3sd(deltaIP_Em_males);
    % females
    [~,~,deltaIP_Ep_females] = rmv_outliers_3sd(deltaIP_Ep_females);
    [~,~,deltaIP_Em_females] = rmv_outliers_3sd(deltaIP_Em_females);
end % remove outliers

%% compare data
% HE choices
[~,choices.HE.pval] = ttest2(choice_hE_males.EpEm, choice_hE_females.EpEm);
[choices.HE.m_males,...
    choices.HE.sem_males] = mean_sem_sd(choice_hE_males.EpEm, 2);
[choices.HE.m_females,...
    choices.HE.sem_females] = mean_sem_sd(choice_hE_females.EpEm, 2);
% HPE choices
[~,choices.HPE.pval] = ttest2(choice_hE_males.Ep, choice_hE_females.Ep);
[choices.HPE.m_males,...
    choices.HPE.sem_males] = mean_sem_sd(choice_hE_males.Ep, 2);
[choices.HPE.m_females,...
    choices.HPE.sem_females] = mean_sem_sd(choice_hE_females.Ep, 2);
% HME choices
[~,choices.HME.pval] = ttest2(choice_hE_males.Em, choice_hE_females.Em);
[choices.HME.m_males,...
    choices.HME.sem_males] = mean_sem_sd(choice_hE_males.Em, 2);
[choices.HME.m_females,...
    choices.HME.sem_females] = mean_sem_sd(choice_hE_females.Em, 2);


% behavioral parameters
for iP = 1:nPrm
    prm_nm = prm_names{iP};
    [~,parameters.(prm_nm).pval] = ttest2(prm_males.(prm_nm),...
        prm_females.(prm_nm));
    [parameters.(prm_nm).m_males,...
        parameters.(prm_nm).sem_males] = mean_sem_sd(prm_males.(prm_nm), 2);
    [parameters.(prm_nm).m_females,...
        parameters.(prm_nm).sem_females] = mean_sem_sd(prm_females.(prm_nm), 2);
    
    % store significant
    if parameters.(prm_nm).pval < 0.05
        parameters.signif_pval.(prm_nm) = parameters.(prm_nm).pval;
    end
end % parameter loop


% RT choices per task
% HE choices
[~,RT.EpEm.pval] = ttest2(RT_summary_males.EpEm.mean_RT.allSubs.choice,...
    RT_summary_females.EpEm.mean_RT.allSubs.choice);
[RT.EpEm.m_males,...
    RT.EpEm.sem_males] = mean_sem_sd(RT_summary_males.EpEm.mean_RT.allSubs.choice, 2);
[RT.EpEm.m_females,...
    RT.EpEm.sem_females] = mean_sem_sd(RT_summary_females.EpEm.mean_RT.allSubs.choice, 2);
% HPE choices
[~,RT.Ep.pval] = ttest2(RT_summary_males.Ep.mean_RT.allSubs.choice,...
    RT_summary_females.Ep.mean_RT.allSubs.choice);
[RT.Ep.m_males,...
    RT.Ep.sem_males] = mean_sem_sd(RT_summary_males.Ep.mean_RT.allSubs.choice, 2);
[RT.Ep.m_females,...
    RT.Ep.sem_females] = mean_sem_sd(RT_summary_females.Ep.mean_RT.allSubs.choice, 2);
% HME choices
[~,RT.Em.pval] = ttest2(RT_summary_males.Em.mean_RT.allSubs.choice,...
    RT_summary_females.Em.mean_RT.allSubs.choice);
[RT.Em.m_males,...
    RT.Em.sem_males] = mean_sem_sd(RT_summary_males.Em.mean_RT.allSubs.choice, 2);
[RT.Em.m_females,...
    RT.Em.sem_females] = mean_sem_sd(RT_summary_females.Em.mean_RT.allSubs.choice, 2);


% indifference point
% physical effort
[~,deltaIP.Ep.pval] = ttest2(deltaIP_Ep_males, deltaIP_Ep_females);
[deltaIP.Ep.m_males,...
    deltaIP.Ep.sem_males] = mean_sem_sd(deltaIP_Ep_males, 2);
[deltaIP.Ep.m_females,...
    deltaIP.Ep.sem_females] = mean_sem_sd(deltaIP_Ep_females, 2);
% mental effort
[~,deltaIP.Em.pval] = ttest2(deltaIP_Em_males, deltaIP_Em_females);
[deltaIP.Em.m_males,...
    deltaIP.Em.sem_males] = mean_sem_sd(deltaIP_Em_males, 2);
[deltaIP.Em.m_females,...
    deltaIP.Em.sem_females] = mean_sem_sd(deltaIP_Em_females, 2);

%% figures
if fig_disp == 1
    %% general parameters for figures
    [pSize, lW, col, mSize] = general_fig_prm;
    female_col = col.red;
    male_col = col.blue_dark;
    tasks = {'EpEm','Ep','Em'};
    nTasks = length(tasks);

    %% choices
    fig;
    for iT = 1:nTasks
        task_nm = tasks{iT};
        % x coordinates
        jPos_male = 1 + 2*(iT - 1);
        jPos_female = 2 + 2*(iT - 1);

        % show male vs female data
        ok_males = ~isnan(choice_hE_males.(task_nm));
        male_violin = Violin({choice_hE_males.(task_nm)(ok_males)},jPos_male,...
            'ViolinColor',{male_col});
        ok_females = ~isnan(choice_hE_females.(task_nm));
        female_violin = Violin({choice_hE_females.(task_nm)(ok_females)},jPos_female,...
            'ViolinColor',{female_col});

        % add p.value indication if difference is significant
        switch task_nm
            case 'EpEm'
                pval_tmp = choices.HE.pval;
            case 'Ep'
                pval_tmp = choices.HPE.pval;
            case 'Em'
                pval_tmp = choices.HME.pval;
        end
        [l_hdl, star_hdl] = add_pval_comparison(choice_hE_males.(task_nm),...
            choice_hE_females.(task_nm),...
            pval_tmp, jPos_male, jPos_female, 'NS');
        ylim([0 120]);
        ylabel('Choices (%)');
        xticks(1.5:2:nTasks*2);
        xticklabels({'HE','HPE','HME'});
        legend_size(pSize);
    end % task loop


    %% behavioral parameters
    fig;
    for iP = 1:nPrm
        prm_nm = prm_names{iP};
        % x coordinates
        jPos_male = 1 + 2*(iP - 1);
        jPos_female = 2 + 2*(iP - 1);

        % show male vs female data
        ok_males = ~isnan(prm_males.(prm_nm));
        male_violin = Violin({prm_males.(prm_nm)(ok_males)},jPos_male,...
            'ViolinColor',{male_col});
        ok_females = ~isnan(prm_females.(prm_nm));
        female_violin = Violin({prm_females.(prm_nm)(ok_females)},jPos_female,...
            'ViolinColor',{female_col});

        % add p.value indication if difference is significant
        [l_hdl, star_hdl] = add_pval_comparison(prm_males.(prm_nm),...
            prm_females.(prm_nm),...
            parameters.(prm_nm).pval, jPos_male, jPos_female, 'NS');
        ylabel('Parameters');
        xticks(1.5:2:nPrm*2);
        xticklabels(prm_names);
        legend_size(pSize);
    end % task loop

    %% RT
    fig;
    for iT = 1:nTasks
        task_nm = tasks{iT};
        % x coordinates
        jPos_male = 1 + 2*(iT - 1);
        jPos_female = 2 + 2*(iT - 1);

        % show male vs female data
        ok_males = ~isnan(RT_summary_males.(task_nm).mean_RT.allSubs.choice);
        male_violin = Violin({RT_summary_males.(task_nm).mean_RT.allSubs.choice(ok_males)},jPos_male,...
            'ViolinColor',{male_col});
        ok_females = ~isnan(RT_summary_females.(task_nm).mean_RT.allSubs.choice);
        female_violin = Violin({RT_summary_females.(task_nm).mean_RT.allSubs.choice(ok_females)},jPos_female,...
            'ViolinColor',{female_col});

        % add p.value indication if difference is significant
        [l_hdl, star_hdl] = add_pval_comparison(RT_summary_males.(task_nm).mean_RT.allSubs.choice,...
            RT_summary_females.(task_nm).mean_RT.allSubs.choice,...
            RT.(task_nm).pval, jPos_male, jPos_female, 'NS');
        ylabel('RT (s)');
        xticks(1.5:2:nTasks*2);
        xticklabels({'HE','HPE','HME'});
        legend_size(pSize);
    end % task loop

    %% IP
    fig;

    % show male vs female data
    ok_Ep_males = ~isnan(deltaIP_Ep_males);
    male_violin = Violin({deltaIP_Ep_males(ok_Ep_males)},1,...
        'ViolinColor',{male_col});
    ok_Ep_females = ~isnan(deltaIP_Ep_females);
    female_violin = Violin({deltaIP_Ep_females(ok_Ep_females)},2,...
        'ViolinColor',{female_col});
    ok_Em_males = ~isnan(deltaIP_Em_males);
    male_violin = Violin({deltaIP_Em_males(ok_Em_males)},3,...
        'ViolinColor',{male_col});
    ok_Em_females = ~isnan(deltaIP_Em_females);
    female_violin = Violin({deltaIP_Em_females(ok_Em_females)},4,...
        'ViolinColor',{female_col});

    % add p.value indication if difference is significant
    [l_hdl, star_hdl] = add_pval_comparison(deltaIP_Ep_males,...
        deltaIP_Ep_females,...
        deltaIP.Ep.pval, 1, 2, 'NS');
    [l_hdl, star_hdl] = add_pval_comparison(deltaIP_Em_males,...
        deltaIP_Em_females,...
        deltaIP.Em.pval, 3, 4, 'NS');
    ylabel('δIP (CHF)');
    xticks([1.5,3.5]);
    xticklabels({'Ep','Em'});
    legend_size(pSize);
end % figure

end % function