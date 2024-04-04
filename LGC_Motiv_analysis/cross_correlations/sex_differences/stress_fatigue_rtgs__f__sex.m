function[stress_rtgs, fatigue_rtgs] = stress_fatigue_rtgs__f__sex(fig_disp, rmv_outliers_yn)
% [stress_rtgs, fatigue_rtgs] = stress_fatigue_rtgs__f__sex(fig_disp, rmv_outliers_yn)
% stress_fatigue_rtgs__f__sex will look at stress and fatigue ratings
% across the experiment duration in function of sex.
%
% INPUTS
% fig_disp: display figures? (1) yes (0) no
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% stress_rtgs: structure comparing stress ratings between male and female
%
% fatigue_rtgs: structure comparing fatigue ratings between male and female

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

%% load stress and fatigue ratings on the day of the experiment
n_TimeSamples = 4;
[stress_males, fatigue_males] = deal(NaN(n_TimeSamples, male_NS));
[stress_females, fatigue_females] = deal(NaN(n_TimeSamples, female_NS));
[~,...
    stress_males(1,:), stress_males(2,:),...
    stress_males(3,:), stress_males(4,:)] = extract_subjective_stress_ratings(study_nm, male_CIDS, male_NS);
[~,...
    stress_females(1,:), stress_females(2,:),...
    stress_females(3,:), stress_females(4,:)] = extract_subjective_stress_ratings(study_nm, female_CIDS, female_NS);
[~,...
    fatigue_males(1,:), fatigue_males(2,:),...
    fatigue_males(3,:), fatigue_males(4,:)] = extract_subjective_fatigue_ratings(study_nm, male_CIDS, male_NS);
[~,...
    fatigue_females(1,:), fatigue_females(2,:),...
    fatigue_females(3,:), fatigue_females(4,:)] = extract_subjective_fatigue_ratings(study_nm, female_CIDS, female_NS);

%% extract mean and SEM per sex and perform the comparison
[stress_rtgs.m_males, stress_rtgs.sem_males] = mean_sem_sd(stress_males, 2);
[stress_rtgs.m_females, stress_rtgs.sem_females] = mean_sem_sd(stress_females, 2);
[fatigue_rtgs.m_males, fatigue_rtgs.sem_males] = mean_sem_sd(fatigue_males, 2);
[fatigue_rtgs.m_females, fatigue_rtgs.sem_females] = mean_sem_sd(fatigue_females, 2);

% perform a repeated measures ANOVA to compare male/females and see if
% there is any effect of time that would be different between the two
% groups
sex = [repmat({'m'},male_NS,1); repmat({'f'},female_NS,1)];

% test stress ratings
stress_ratings_pooled = [stress_males'; stress_females'];
t_stress_rtg = table(sex, stress_ratings_pooled(:,1), stress_ratings_pooled(:,2),stress_ratings_pooled(:,3),stress_ratings_pooled(:,4),...
    'VariableNames',{'sex','preR1','postR1','preR2','postR2'});
timePoints = table([1 2 3 4]','VariableNames',{'Time'});
rm_stress = fitrm(t_stress_rtg,'preR1-postR2~sex','WithinDesign',timePoints);
ranovatbl_stress = ranova(rm_stress);
% post-hoc test comparing males and females for each timepoint
comp_stress = multcompare(rm_stress,'sex','By','Time');
pval_stress = comp_stress.pValue(1:2:end);
stress_rtgs.pval = pval_stress;
stress_rtgs.anova_table = ranovatbl_stress;
% margmean(rm_stress,'Time'); % if you want to look at direction of the effect
% you can extract the marginal mean

% test fatigue ratings
fatigue_ratings_pooled = [fatigue_males'; fatigue_females'];
t_fatigue_rtg = table(sex, fatigue_ratings_pooled(:,1), fatigue_ratings_pooled(:,2),fatigue_ratings_pooled(:,3),fatigue_ratings_pooled(:,4),...
    'VariableNames',{'sex','preR1','postR1','preR2','postR2'});
timePoints = table([1 2 3 4]','VariableNames',{'Time'});
rm_fatigue = fitrm(t_fatigue_rtg,'preR1-postR2~sex','WithinDesign',timePoints);
ranovatbl_fatigue = ranova(rm_fatigue);
% post-hoc test comparing males and females for each timepoint
comp_fatigue = multcompare(rm_fatigue,'sex','By','Time');
pval_fatigue = comp_fatigue.pValue(1:2:end);
fatigue_rtgs.anova_table = ranovatbl_fatigue;
% margmean(rm_fatigue,'Time'); % if you want to look at direction of the effect
% you can extract the marginal mean
fatigue_rtgs.pval = pval_fatigue;

%% figure display
if fig_disp == 1
    %% general parameters for figures
    [pSize, lW, col, mSize] = general_fig_prm;
    male_col = col.blue_dark;
    female_col = col.red;
    
    %% stress ratings
    fig;
    
    for iT = 1:n_TimeSamples
        jR_male = 1 + 2*(iT - 1);
        jR_female = 2 + 2*(iT - 1);
        
        % show male vs female data for max perf
        ok_males = ~isnan(stress_males(iT,:));
        male_violin = Violin({stress_males(iT, ok_males)},jR_male,...
            'ViolinColor',{male_col});
        ok_females = ~isnan(stress_females(iT,:));
        female_violin = Violin({stress_females(iT,ok_females)},jR_female,...
            'ViolinColor',{female_col});
        
        % add p.value indication if difference is significant
        [l_hdl, star_hdl] = add_pval_comparison(stress_males(iT,:),...
            stress_females(iT,:),...
            pval_stress(iT), jR_male, jR_female, 'NS');
    end % run loop
    ylabel('Stress rating');
    xticks(1.5:2:(n_TimeSamples*2));
    xticklabels({'pre-MRS','post-MRS','pre-MRI','post-MRI'});
    legend_size(pSize);
    
    %% fatigue ratings
    fig;
    
    for iT = 1:n_TimeSamples
        jR_male = 1 + 2*(iT - 1);
        jR_female = 2 + 2*(iT - 1);
        
        % show male vs female data for max perf
        ok_males = ~isnan(fatigue_males(iT,:));
        male_violin = Violin({fatigue_males(iT, ok_males)},jR_male,...
            'ViolinColor',{male_col});
        ok_females = ~isnan(fatigue_females(iT,:));
        female_violin = Violin({fatigue_females(iT,ok_females)},jR_female,...
            'ViolinColor',{female_col});
        
        % add p.value indication if difference is significant
        [l_hdl, star_hdl] = add_pval_comparison(fatigue_males(iT,:),...
            fatigue_females(iT,:),...
            pval_fatigue(iT), jR_male, jR_female, 'NS');
    end % run loop
    ylabel('Fatigue rating');
    xticks(1.5:2:(n_TimeSamples*2));
    xticklabels({'pre-MRS','post-MRS','pre-MRI','post-MRI'});
    legend_size(pSize);
end % figure display
end % function