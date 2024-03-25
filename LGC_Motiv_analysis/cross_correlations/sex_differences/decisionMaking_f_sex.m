function[choices, parameters, RT] = decisionMaking_f_sex()
% [choices, parameters, RT] = decisionMaking_f_sex()
% decisionMaking_f_sex will compare all the variables related to the
% decision-making process of the experiment between males and females.
% Variables include the total proportion of choices (HE), the proportion of
% physical (HPE) or mental (HME) choices, and all the parameters derived
% from our computational model (kR, kP, kEp, kEm, kFp, kLm, bias).
%
% OUTPUTS
% choices: structure comparing choice-related variables between male and female
%
% parameters: structure comparing behavioral parameters between male and female
%
% RT: structure comparing RT between male and female

%% subject selection
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

%% load variables of interest
% load proportion of high-effort choices
fig_disp = 0;
[choice_hE_males] = choice_hE_proportion(study_nm, condition, male_CIDS, fig_disp);
[choice_hE_females] = choice_hE_proportion(study_nm, condition, female_CIDS, fig_disp);

% load behavioral parameters
[prm_males, mdlType, mdlN] = prm_extraction(study_nm, male_CIDS);
[prm_females, ~, ~] = prm_extraction(study_nm, female_CIDS);

% load RT
[RT_summary_males] = RT_range(male_CIDS, condition, fig_disp);
[RT_summary_females] = RT_range(female_CIDS, condition, fig_disp);

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
prm_names = {'kR','kP','kEp','kEm',...
    'kBiasM','kFp','kLm','kR_div_kP'};
nPrm = length(prm_names);
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

end % function