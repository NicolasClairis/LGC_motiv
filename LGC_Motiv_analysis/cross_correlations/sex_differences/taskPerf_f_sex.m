% function[calib] = taskPerf_f_sex(fig_disp)
%[calib] = taskPerf_f_sex(fig_disp)
% perf_f_sex will compare different variables related to performance
% (maximal performance during calibration, average RT, difference in
% accuracy, overshoot, etc.) between males and females
%
% INPUTS
% fig_disp: display figure (1) or not (0)
%
% OUTPUTS
% calib: structure with calibration variables for physical (MVC + PCSA) and
% mental (NMP) tasks
%
%

%% define inputs by default
% display figure by default
if ~exist('fig_disp','var') || isempty(fig_disp) || ~ismember(fig_disp,[0,1])
    fig_disp = 1;
end

%% subject selection
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

%% extract performance variables

% calibration
% extract male data
[male_MVC, male_PCSA] = extract_grip_MVC_PCSA(study_nm, male_CIDS, male_NS);
[male_NMP] = extract_Nback_NMP(study_nm, male_CIDS, male_NS);
% extract female data
[female_MVC, female_PCSA] = extract_grip_MVC_PCSA(study_nm, female_CIDS, female_NS);
[female_NMP] = extract_Nback_NMP(study_nm, female_CIDS, female_NS);

% max performance before/after each run
nMaxPerf = 4;
maxPerf_males = maxPerfEvolutionAcrossRuns_group([], 0, 0, study_nm, condition, male_CIDS, male_NS);
maxPerf_females = maxPerfEvolutionAcrossRuns_group([], 0, 0, study_nm, condition, female_CIDS, female_NS);

% load RT (performance)
fig_disp0 = 0;
[RT_summary_males] = RT_range(male_CIDS, condition, fig_disp0);
[RT_summary_females] = RT_range(female_CIDS, condition, fig_disp0);

% accuracy
% physical: overshoot + time% below threshold + final perf
% mental: n.correct/n.errors (% correct ~= perf) and n.errors

%% perform comparisons

% calibration
% max performance in physical task (in Newtons)
[~,calib.MVC.pval] = ttest2(male_MVC, female_MVC);
[calib.MVC.m_males, calib.MVC.sem_males] = mean_sem_sd(male_MVC,2);
[calib.MVC.m_females, calib.MVC.sem_females] = mean_sem_sd(female_MVC,2);

% theoretical Fmax (based on forearm values)
[~,calib.PCSA.pval] = ttest2(male_PCSA, female_PCSA);
[calib.PCSA.m_males, calib.PCSA.sem_males] = mean_sem_sd(male_PCSA,2);
[calib.PCSA.m_females, calib.PCSA.sem_females] = mean_sem_sd(female_PCSA,2);

% max performance in mental effort task
[~,calib.NMP.pval] = ttest2(male_NMP, female_NMP);
[calib.NMP.m_males, calib.NMP.sem_males] = mean_sem_sd(male_NMP,2);
[calib.NMP.m_females, calib.NMP.sem_females] = mean_sem_sd(female_NMP,2);


% max perf before/after each run
% perform a repeated measures ANOVA to compare male/females and see if
% there is any effect of time that would be different between the two
% groups
sex = [repmat({'m'},male_NS,1); repmat({'f'},female_NS,1)];
% test physical max
maxP_Ep = [maxPerf_males.Ep.allData'; maxPerf_females.Ep.allData'];
t_maxP_Ep = table(sex, maxP_Ep(:,1), maxP_Ep(:,2),maxP_Ep(:,3),maxP_Ep(:,4),...
    'VariableNames',{'sex','preR1','postR1','preR2','postR2'});
timePoints = table([1 2 3 4]','VariableNames',{'Time'});
rm_Ep = fitrm(t_maxP_Ep,'preR1-postR2~sex','WithinDesign',timePoints);
ranovatbl_Ep = ranova(rm_Ep);
% post-hoc test comparing males and females for each timepoint
comp_Ep = multcompare(rm_Ep,'sex','By','Time');
pval_maxP_Ep = comp_Ep.pValue(1:2:end);
% margmean(rm_Ep,'Time'); % if you want to look at direction of the effect
% you can extract the marginal mean

% test mental max
maxP_Em = [maxPerf_males.Em.allData'; maxPerf_females.Em.allData'];
t_maxP_Em = table(sex, maxP_Em(:,1), maxP_Em(:,2),maxP_Em(:,3),maxP_Em(:,4),...
    'VariableNames',{'sex','preR1','postR1','preR2','postR2'});
timePoints = table([1 2 3 4]','VariableNames',{'Time'});
rm_Em = fitrm(t_maxP_Em,'preR1-postR2~sex','WithinDesign',timePoints);
ranovatbl_Em = ranova(rm_Em);
% post-hoc test comparing males and females for each timepoint
comp_Em = multcompare(rm_Em,'sex','By','Time');
pval_maxP_Em = comp_Em.pValue(1:2:end);
% margmean(rm_Em,'Time'); % if you want to look at direction of the effect
% you can extract the marginal mean


% RT
% RT choices per task
% across all efforts
[~,RT.EpEm.pval] = ttest2(RT_summary_males.EpEm.mean_RT.allSubs.perf,...
    RT_summary_females.EpEm.mean_RT.allSubs.perf);
[RT.EpEm.m_males,...
    RT.EpEm.sem_males] = mean_sem_sd(RT_summary_males.EpEm.mean_RT.allSubs.perf, 2);
[RT.EpEm.m_females,...
    RT.EpEm.sem_females] = mean_sem_sd(RT_summary_females.EpEm.mean_RT.allSubs.perf, 2);
% physical efforts
[~,RT.Ep.pval] = ttest2(RT_summary_males.Ep.mean_RT.allSubs.perf,...
    RT_summary_females.Ep.mean_RT.allSubs.perf);
[RT.Ep.m_males,...
    RT.Ep.sem_males] = mean_sem_sd(RT_summary_males.Ep.mean_RT.allSubs.perf, 2);
[RT.Ep.m_females,...
    RT.Ep.sem_females] = mean_sem_sd(RT_summary_females.Ep.mean_RT.allSubs.perf, 2);
% mental efforts
[~,RT.Em.pval] = ttest2(RT_summary_males.Em.mean_RT.allSubs.perf,...
    RT_summary_females.Em.mean_RT.allSubs.perf);
[RT.Em.m_males,...
    RT.Em.sem_males] = mean_sem_sd(RT_summary_males.Em.mean_RT.allSubs.perf, 2);
[RT.Em.m_females,...
    RT.Em.sem_females] = mean_sem_sd(RT_summary_females.Em.mean_RT.allSubs.perf, 2);

%% figure display
if fig_disp == 1
    %% general parameters for figures
    [pSize, lW, col, mSize] = general_fig_prm;
    male_col = col.blue_dark;
    female_col = col.red;
    tasks = {'EpEm','Ep','Em'};
    nTasks = length(tasks);
    
    %% calibration
    fig;
    
    % physical effort
    subplot(1,2,1); hold on;
    % show male vs female data for MVC
    ok_MVC_males = ~isnan(male_MVC);
    male_violin = Violin({male_MVC(ok_MVC_males)},1,...
        'ViolinColor',{male_col});
    ok_MVC_females = ~isnan(female_MVC);
    female_violin = Violin({female_MVC(ok_MVC_females)},2,...
        'ViolinColor',{female_col});
    
    % show male vs female data for MVC
    ok_PCSA_males = ~isnan(male_PCSA);
    male_violin = Violin({male_PCSA(ok_PCSA_males)},3,...
        'ViolinColor',{male_col});
    ok_PCSA_females = ~isnan(female_PCSA);
    female_violin = Violin({female_PCSA(ok_PCSA_females)},4,...
        'ViolinColor',{female_col});
    
    % add p.value indication if difference is significant
    [l_hdl, star_hdl] = add_pval_comparison(male_MVC,...
        female_MVC,...
        calib.MVC.pval, 1, 2, 'NS');
    [l_hdl, star_hdl] = add_pval_comparison(male_PCSA,...
        female_PCSA,...
        calib.PCSA.pval, 3, 4, 'NS');
    
    ylabel('Force calibration (N)');
    xticks([1.5, 3.5]);
    xticklabels({'MVC','PCSA'});
    legend_size(pSize);
    
    % mental effort
    subplot(1,2,2);
    % show male vs female data for NMP
    ok_NMP_males = ~isnan(male_NMP);
    male_violin = Violin({male_NMP(ok_NMP_males)},1,...
        'ViolinColor',{male_col});
    ok_NMP_females = ~isnan(female_NMP);
    female_violin = Violin({female_NMP(ok_NMP_females)},2,...
        'ViolinColor',{female_col});
    
    % add p.value indication if difference is significant
    [l_hdl, star_hdl] = add_pval_comparison(male_NMP,...
        female_NMP,...
        calib.NMP.pval, 1, 2, 'NS');
    
    ylabel('2-back max calibration');
    legend_size(pSize);
    
    %% max performance
    fig;
    
    % physical
    subplot(1,2,1); hold on;
    for iP = 1:nMaxPerf
        jR_male = 1 + 2*(iP - 1);
        jR_female = 2 + 2*(iP - 1);
        
        % show male vs female data for max perf
        ok_maxP_Ep_males = ~isnan(maxPerf_males.Ep.allData(iP,:));
        male_violin = Violin({maxPerf_males.Ep.allData(iP, ok_maxP_Ep_males)},jR_male,...
            'ViolinColor',{male_col});
        ok_maxP_Ep_females = ~isnan(maxPerf_females.Ep.allData(iP,:));
        female_violin = Violin({maxPerf_females.Ep.allData(iP,ok_maxP_Ep_females)},jR_female,...
            'ViolinColor',{female_col});
        
        % add p.value indication if difference is significant
        [l_hdl, star_hdl] = add_pval_comparison(maxPerf_males.Ep.allData(iP,:),...
            maxPerf_females.Ep.allData(iP,:),...
            pval_maxP_Ep(iP), jR_male, jR_female, 'NS');
    end % run loop
    ylabel('Performance (%)');
    xticks(1.5:2:(nMaxPerf*2));
    xticklabels({'pre-r1','post-r1','pre-r2','post-r2'});
    legend_size(pSize);
    
    % mental
    subplot(1,2,2); hold on;
    for iP = 1:nMaxPerf
        jR_male = 1 + 2*(iP - 1);
        jR_female = 2 + 2*(iP - 1);
        
        % show male vs female data for max perf
        ok_maxP_Em_males = ~isnan(maxPerf_males.Em.allData(iP,:));
        male_violin = Violin({maxPerf_males.Em.allData(iP, ok_maxP_Em_males)},jR_male,...
            'ViolinColor',{male_col});
        ok_maxP_Em_females = ~isnan(maxPerf_females.Em.allData(iP,:));
        female_violin = Violin({maxPerf_females.Em.allData(iP,ok_maxP_Em_females)},jR_female,...
            'ViolinColor',{female_col});
        
        % add p.value indication if difference is significant
        [l_hdl, star_hdl] = add_pval_comparison(maxPerf_males.Em.allData(iP,:),...
            maxPerf_females.Em.allData(iP,:),...
            pval_maxP_Em(iP), jR_male, jR_female, 'NS');
    end % run loop
    ylabel('Performance (%)');
    xticks(1.5:2:(nMaxPerf*2));
    xticklabels({'pre-r1','post-r1','pre-r2','post-r2'});
    legend_size(pSize);
    
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
        xticklabels({'E','Ep','Em'});
        legend_size(pSize);
    end % task loop
end % figure display

% end % function